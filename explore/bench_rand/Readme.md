# Readme

Various benchmarks for optimization purposes.

To run the samples, install `libabsl-dev` and `libtestu01-0-dev` on Ubuntu.

TODO:

- Using VMT19937 commit [fa93111](https://github.com/fabiocannizzo/VMT19937/commit/fa93111bfc6f56f25c315990430a3487cdff9935). However, this random generator requires bundling of file sized ~48 MiB INTO THE BINARY.
- File `vigna.h` is from [xoshiro](https://github.com/nessan/xoshiro/) under tag [v1.0.0](https://github.com/nessan/xoshiro/releases/tag/v1.1.0) and commit [176fa19](https://github.com/nessan/xoshiro/commit/176fa191c8493e4c5cb06a44bc083010664fe39b).

- File from <https://www.pcg-random.org/posts/some-prng-implementations.html>:
  - [`jsf.hpp`](https://gist.github.com/imneme/85cff47d4bad8de6bdeb671f9c76c814)
  - [`gjrand.hpp`](https://gist.github.com/imneme/7a783e20f71259cc13e219829bcea4ac)
  - [`sfc.hpp`](https://gist.github.com/imneme/f1f7821f07cf76504a97f6537c818083)
  - [`lehmer.hpp`](https://gist.github.com/imneme/aeae7628565f15fb3fef54be8533e39c)
  - [`splitmix.hpp`](https://gist.github.com/imneme/6179748664e88ef3c34860f44309fc71)
  - [`arc4.hpp`](https://gist.github.com/imneme/4f2bf4b4f3a221ef051cf108d6b64d5a)
