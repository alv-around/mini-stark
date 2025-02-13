# mini-stark

In this repo I implement a stark prover for learning purposes using arkworks.
**Do not use this in a production environment.**

🚧 This repo is WIP 🚧 I intend to add more features with time.

Issues/comments/critics all very welcome.

## Roadmap

Planned features(in no specific order):

- [ ] add Fiat-Shamir
- [ ] commit to multiple polynomials
- [ ] add support for AIR arithmetization
- [x] extend `merkle.rs` to get trees with 4/8 nodes
- [ ] Add `rayon` for parallelization
- [ ] Add `criterion` and integration tests to benchmark performance
- [ ] perf optimizations such as:
  - [ ] Sort $\Omega$ for less merkle paths
  - [ ] griding
  - [ ] STIR ??
  - [ ] WHIR ??
  - [ ] ..
