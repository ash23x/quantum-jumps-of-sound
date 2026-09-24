// ---- Quantum-jumps-of-sound toy model: physics core (plain JS, no deps) ----
function mulberry32(a){return function(){a|=0;a=a+0x6D2B79F5|0;let t=Math.imul(a^a>>>15,1|a);t=t+Math.imul(t^t>>>7,61|t)^t;return((t^t>>>14)>>>0)/4294967296;};}
const H_PLANCK=6.62607015e-34, Q_E=1.602176634e-19;

// derived measurement parameters from the knobs
function derive(p){
  const tsel = 1/p.chi2;                 // number-selective pulse: just long enough to pick out one phonon line
  const tc = tsel + p.tro;               // one reading = selective pulse + readout and reset
  const pdec = 1-Math.exp(-(0.5*tsel+0.5*p.tro)/p.T1q); // qubit relaxes before it is read
  const wq = 1/(Math.PI*p.T1q);          // qubit FWHM, Hz (T2 taken equal to T1)
  const w1 = wq + 1/(2*Math.PI*p.T1m);   // FWHM of the one-phonon line
  // pulse window ~ chi2 wide centred on the one-phonon line (Lorentzian overlap model)
  const eta1 = (2/Math.PI)*Math.atan(p.chi2/w1);                                   // n=1 line caught by the window
  const x0 = (Math.atan(3*p.chi2/wq)-Math.atan(p.chi2/wq))/Math.PI;               // n=0 line leaking into it
  const e1 = 1-(1-p.eps)*(1-pdec)*eta1;  // P(no click | one phonon)
  const e0 = p.eps+(1-p.eps)*(1-pdec)*x0; // P(click | no phonon): leakage clicks decay too
  const perLife = p.T1m/tc;
  const resolved = p.chi2 >= 2*w1;
  const D = (1-e1)*Math.log((1-e1)/e0) + e1*Math.log(e1/(1-e0)); // nats per reading
  const nReq = D>0 ? Math.log(99)/D : Infinity;                  // readings for 99:1 odds
  const visible = perLife >= 5*nReq;
  const T1eff = p.T1m/(1+3*p.nth);       // a warm one-phonon state can also leave upwards
  return {tsel,tc,pdec,e1,e0,wq,w1,eta1,x0,perLife,resolved,D,nReq,visible,T1eff};
}

function generator(kappa,nth,N){
  const Q=[];for(let i=0;i<N;i++)Q.push(new Float64Array(N));
  for(let n=0;n<N;n++){
    const down=n>0?kappa*n*(nth+1):0, up=n<N-1?kappa*(n+1)*nth:0;
    if(n>0)Q[n][n-1]=down; if(n<N-1)Q[n][n+1]=up; Q[n][n]=-(down+up);
  }
  return Q;
}
function matmul(A,B){const N=A.length,C=[];for(let i=0;i<N;i++){const r=new Float64Array(N);for(let k=0;k<N;k++){const a=A[i][k];if(a===0)continue;for(let j=0;j<N;j++)r[j]+=a*B[k][j];}C.push(r);}return C;}
function eye(N){const I=[];for(let i=0;i<N;i++){const r=new Float64Array(N);r[i]=1;I.push(r);}return I;}
function expm(Q,t){
  const N=Q.length;let norm=0;
  for(let i=0;i<N;i++){let s=0;for(let j=0;j<N;j++)s+=Math.abs(Q[i][j]*t);norm=Math.max(norm,s);}
  const s=norm>0.25?Math.ceil(Math.log2(norm/0.25)):0, f=t/Math.pow(2,s);
  const A=Q.map(r=>r.map(v=>v*f));
  let term=eye(N),res=eye(N);
  for(let k=1;k<=20;k++){term=matmul(term,A);for(let i=0;i<N;i++)for(let j=0;j<N;j++){term[i][j]/=k;res[i][j]+=term[i][j];}}
  for(let i=0;i<s;i++)res=matmul(res,res);
  for(let i=0;i<N;i++){let z=0;for(let j=0;j<N;j++){if(res[i][j]<0)res[i][j]=0;z+=res[i][j];}for(let j=0;j<N;j++)res[i][j]/=z;}
  return res;
}

// one run: truth (Gillespie on Fock states), clicks, forward (real-time) + smoothed posteriors
function oneRun(p,d,K,N,P0,P,rng){
  const kappa=1/p.T1m, nth=p.nth;
  const truth=new Uint8Array(K), clicks=new Uint8Array(K);
  let n = rng()<p.pload?1:0, t=0;
  for(let k=0;k<K;k++){
    const tk=d.tsel/2+k*d.tc;
    for(;;){
      const down=n>0?kappa*n*(nth+1):0, up=n<N-1?kappa*(n+1)*nth:0, R=down+up;
      if(R<=0){t=tk;break;}
      const dt=-Math.log(1-rng())/R;
      if(t+dt>tk){t=tk;break;}
      t+=dt; n=(rng()*R<down)?n-1:n+1;
    }
    truth[k]=n;
    clicks[k]= rng() < (n===1?1-d.e1:d.e0) ? 1:0;
  }
  // forward
  const F=new Float64Array(K*N), fwd1=new Float64Array(K);
  let a=new Float64Array(N); a[0]=1-p.pload; a[1]=p.pload;
  for(let k=0;k<K;k++){
    const M=k===0?P0:P, b=new Float64Array(N);
    for(let i=0;i<N;i++){const ai=a[i];if(ai===0)continue;const r=M[i];for(let j=0;j<N;j++)b[j]+=ai*r[j];}
    const c=clicks[k];let z=0;
    for(let j=0;j<N;j++){const L=j===1?(c?1-d.e1:d.e1):(c?d.e0:1-d.e0);b[j]*=L;z+=b[j];}
    if(!(z>0)){b.fill(1/N);z=1;} else for(let j=0;j<N;j++)b[j]/=z;
    a=b; F.set(a,k*N); fwd1[k]=a[1];
  }
  // backward -> smoothed P(n=1)
  const sm1=new Float64Array(K); let be=new Float64Array(N).fill(1);
  for(let k=K-1;k>=0;k--){
    let z=0,s1=0;for(let j=0;j<N;j++){const g=F[k*N+j]*be[j];z+=g;if(j===1)s1=g;}
    sm1[k]=z>0?s1/z:0;
    if(k>0){
      const c=clicks[k],nb=new Float64Array(N);let zz=0;
      const Lb=new Float64Array(N);for(let j=0;j<N;j++){const L=j===1?(c?1-d.e1:d.e1):(c?d.e0:1-d.e0);Lb[j]=L*be[j];}
      for(let i=0;i<N;i++){let s=0;const r=P[i];for(let j=0;j<N;j++)s+=r[j]*Lb[j];nb[i]=s;zz+=s;}
      for(let i=0;i<N;i++)nb[i]/=(zz>0?zz:1); be=nb;
    }
  }
  // herald = the real-time belief reaches 99:1 odds
  let kh=-1;
  for(let k=0;k<K;k++)if(fwd1[k]>=0.99){kh=k;break;}
  return {truth,clicks,fwd1,sm1,kh};
}

function setup(p){
  const d=derive(p), N=p.nth<=0.1?4:6;
  const Q=generator(1/p.T1m,p.nth,N);
  const P=expm(Q,d.tc), P0=expm(Q,d.tsel/2);
  return {d,N,P,P0};
}
function recordLength(p,d,Kmax){ return Math.max(40,Math.min(Kmax,Math.ceil(6*p.T1m/d.tc))); }

// first stay in the one-phonon state after the herald, read off the smoothed posterior.
// A drop only counts once it has lasted W readings, so noise dips are bridged and later
// thermal re-excitations are not added to the first stay.
function firstStay(o,K,d,T1eff){
  if(o.kh<0) return null;
  const W=Math.max(3,Math.ceil(0.1*T1eff/d.tc));
  let kend=-1;
  for(let k=o.kh+1;k<K;k++){
    if(o.sm1[k]<0.5){let ok=true;for(let j=k;j<Math.min(K,k+W);j++){if(o.sm1[j]>=0.5){ok=false;k=j;break;}}if(ok){kend=k;break;}}
  }
  const stop = Math.min(K, kend<0?K:kend+W);
  let dw=0; for(let k=o.kh;k<stop;k++) dw+=o.sm1[k];
  return {stay:Math.max(0,dw*d.tc-d.tc/2), ended:kend>=0, kend};
}
function manyRuns(p,seed){
  const {d,N,P,P0}=setup(p);
  const K=recordLength(p,d,20000);
  const runs=Math.max(60,Math.min(400,Math.floor(3e7/(2*N*N*K))));
  const rng=mulberry32(seed);
  const waits=[];let nh=0,nj=0,sumW=0,good=0;
  for(let r=0;r<runs;r++){
    const o=oneRun(p,d,K,N,P0,P,rng);
    const f=firstStay(o,K,d,d.T1eff); if(!f) continue; nh++; if(o.truth[o.kh]===1) good++;
    sumW+=f.stay; if(f.ended){nj++;waits.push(f.stay);}   // censored stays still add time
  }
  const T1est = nj>0 ? sumW/nj : NaN;
  return {runs,nh,nj,waits,T1est,K,d,heraldFid: nh>0?good/nh:NaN};
}
if(typeof module!=='undefined')module.exports={firstStay,mulberry32,derive,generator,expm,oneRun,setup,recordLength,manyRuns,H_PLANCK,Q_E};
