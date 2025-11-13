# ECQ

`ECQ` is a library based on [FLINT](https://github.com/flintlib/flint/) to do arithmetic with elliptic curves over $\mathbb Q$. In particular, computing rational points on the curve.

## $2$-descent

Complete $2$-descent is work in progress and is the first goal of this project.

### Basic setup

- [x] Basic data structure
- [x] Computation of the possible value for $b_i$

### Ternary quadratic forms

- [x] Reduction of TQF
- [ ] Finding primitive solution with LLL
- [ ] Parameterization of the solution space

### Hyper-elliptic curve of genus $2, 3$

- [ ] Solubility criterion
- [ ] Basic enumeration
- [ ] Sieving for rational points
- [ ] Higher order descent

## Installation

For now, `ECQ` is dependent on `FLINT`, specifically on this [PR](https://github.com/flintlib/flint/pull/2486), I might just add the full content of the PR to `factor_addition` is it ends up not being merged.

To make the current demo :

```bash
make
```

and run it with :

```bash
./ecq
```
