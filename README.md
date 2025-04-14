# mini-stark

In this repo I implement a stark prover for learning purposes using arkworks.
**Do not use this in a production environment.**

🚧 This repo is WIP 🚧 I intend to add more features with time.

currently system status:

- Complete ✅
- Soundness :x:
- Non-interactive ✅
- Succint ✅
- Zero-Knowledge^* ❓

Issues/comments/critics all very welcome.

## Roadmap

Planned features(in no specific order):

- [x] add Fiat-Shamir
- [x] commit to multiple polynomials
- [x] add support for AIR arithmetization
- [x] extend `merkle.rs` to get trees with 4/8 nodes
- [x] add zk: random trace padding and domain coset
- [x] Stark PIOP
  - [ ] add DEEP-ALI to IOPP 🚧
- [ ] Add `rayon` for parallelization
- [ ] Add `criterion` and integration tests to benchmark performance
- [ ] perf optimizations such as:
  - [ ] griding
  - [ ] STIR ??
  - [ ] WHIR ??
  - [ ] ..

## Comments

### Zero-Knowledge

  The zk property was implemented following [stark-by-hand](https://dev.risczero.com/proof-system/stark-by-hand).
  It was shown in the [following paper](https://eprint.iacr.org/2024/1037.pdf#cite.FRISummary) do not meet the requirements to assert the specific property of zero knowledge.
But this is still a topic of discussion.

## Reference

Here is a list of some resources, I used:

- The IOPP has been mostly taken from [stark-by-hand by risc0](https://dev.risczero.com/proof-system/stark-by-hand) and to a lesser degree from [anatomy of a STARK](https://aszepieniec.github.io/stark-anatomy/stark).
- FRI is based [anatomy of a STARK](https://aszepieniec.github.io/stark-anatomy/stark) and [Chp.8 of zk-learning.org](https://rdi.berkeley.edu/zk-learning/assets/lecture8.pdf) and the following [FRI Summary](https://eprint.iacr.org/2022/1216.pdf) in no specific order.
- For the AIR arithmetization, I took great inspiration from [OpenZK](https://www.youtube.com/watch?v=H3AKu03AwYc). [Notes on air arithmetization](https://cronokirby.com/posts/2022/09/notes-on-stark-arithmetization/) is a very good primer on air.
- For the zk-property this implementation follows [stark-by-hand](https://dev.risczero.com/proof-system/stark-by-hand#lesson-5-zk-commitments-of-the-trace-data) and I am aware of it limitations discussed [here](https://eprint.iacr.org/2024/1037.pdf#cite.FRISummary).
- I took FRI-batching from risc0 knowledge database
...

The list is not complete, so if you see recognize something who needs attribution, hit me up :)
