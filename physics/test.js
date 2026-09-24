// Regression check for the physics core: can the lifetime be read back from simulated clicks alone?
// Run: node physics/test.js   (Node 18+)
const C=require('./core.js');
const base={T1m:2.1e-3,chi2:328e3,nth:0.02,T1q:30e-6,eps:0.02,tro:2e-6,pload:0.9};
const cases=[['Stanford 2026 values',{}],['slow qubit, 3 us',{T1q:3e-6}],['short lifetime, 50 us',{T1m:50e-6}],['2010 FBAR lifetime, 6.1 ns',{T1m:6.1e-9}]];
let fail=0;
for(const [name,over] of cases){
  const p={...base,...over}, d=C.derive(p); const ests=[]; let nh=0,runs=0;
  for(let s=1;s<=10;s++){const m=C.manyRuns(p,s*7919);if(isFinite(m.T1est))ests.push(m.T1est);nh+=m.nh;runs+=m.runs;}
  const m=ests.length?ests.reduce((a,b)=>a+b,0)/ests.length:NaN;
  const se=ests.length>1?Math.sqrt(ests.reduce((a,b)=>a+(b-m)**2,0)/(ests.length-1))/Math.sqrt(ests.length):NaN;
  let verdict;
  if(!d.visible&&nh===0) verdict='ok (never heralded, as expected)';
  else { const z=Math.abs(m-d.T1eff)/Math.max(se,1e-12); verdict = z<3||Math.abs(m/d.T1eff-1)<0.05 ? 'ok' : 'FAIL'; if(verdict==='FAIL') fail++; }
  console.log(`${name.padEnd(28)} expected ${(d.T1eff*1e3).toFixed(4)} ms  recovered ${isFinite(m)?(m*1e3).toFixed(4):'   -  '} ms  heralded ${nh}/${runs}  ${verdict}`);
}
process.exit(fail?1:0);
