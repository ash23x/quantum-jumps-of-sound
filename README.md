# Quantum jumps of sound: a working model

An interactive model of the Stanford experiment that watched a single phonon leave a vibrating resonator in real time: T. Makihara, E. Szakiel *et al.*, "Quantum jumps of sound", *Science* (2026), doi:10.1126/science.aeh7535.

## What it does

- Simulates the resonator's phonon number as a jump process on number states, sampled by repeated quantum non-demolition readings from a dispersively coupled qubit.
- Runs a Bayes filter over the noisy readings (the live estimate), heralds a phonon at 99:1 odds, and reads the phonon lifetime back from the clicks alone.
- Maps where devices sit on phonon lifetime against qubit shift per phonon, including the 2010 FBAR (the first mechanical object put in its quantum ground state) and a bench FBAR at Q 1000.
- Turns a Butterworth–Van Dyke fit (C0, Cm) into the qubit coupling g, the qubit nonlinearity α and the shift per phonon 2χ, and shows why a 50 Ω-matched FBAR's own C0 swamps a qubit.

## Use it

Open `index.html` in a browser. There is no build step and nothing to install.

To check the physics: `node physics/test.js` (Node 18 or later). It confirms the lifetime read back from simulated clicks matches the input, and that a 6.1 ns resonator is never heralded.

## What is from the paper and what is assumed

From the paper: a 2.1 ms phonon lifetime, a 328 kHz qubit shift per phonon and 85% herald fidelity. The qubit lifetime, readout error, timings, temperature and the coupling examples are assumptions and are labelled as such on the page. Line shapes and pulse selectivity are Lorentzian estimates; there is no pulse shaping and no measurement back-action on the phonon.

## Sources

1. T. Makihara *et al.*, "Quantum jumps of sound", *Science* (2026), doi:10.1126/science.aeh7535.
2. A. D. O'Connell *et al.*, "Quantum ground state and single-phonon control of a mechanical resonator", *Nature* 464, 697–703 (2010).
3. R. G. Gruenke-Freudenstein *et al.*, "Surface and bulk two-level system losses in lithium niobate acoustic resonators", arXiv:2501.08291 (2025).
4. Y. Chu *et al.*, "Quantum acoustics with superconducting qubits", *Science* 358, 199–202 (2017).
5. J. Koch *et al.*, "Charge-insensitive qubit design derived from the Cooper pair box", *Phys. Rev. A* 76, 042319 (2007).

Built with Claude (Anthropic).
