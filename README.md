# mini-stark

In this repo I implement a stark prover for learning purposes using arkworks.
**Do not use this in a production environment.**

🚧 This repo is WIP 🚧 I intend to add more features with time.

Issues/comments/critics all very welcome.

## Roadmap

Planned features(in no specific order):

- [ ] add Fiat-Shamir
- [x] commit to multiple polynomials
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

## Reference

Here is a list of some resources, I used:

- The IOP and FRI commitment scheme are mostly take from [anatomy of a STARK](https://aszepieniec.github.io/stark-anatomy/stark) and [Chp.8 of zk-learning.org](https://rdi.berkeley.edu/zk-learning/assets/lecture8.pdf).
- For the AIR arithmatization, I took the implementation from [OpenZK](https://www.youtube.com/watch?v=H3AKu03AwYc). [Notes on air arithmetization](https://cronokirby.com/posts/2022/09/notes-on-stark-arithmetization/) is a very good primer on air.
- FRI-batching from risc0 knowledge database
...

The list is not complete, so if you see recognize something who needs attribution, hit me up :)
