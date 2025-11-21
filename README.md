# ECQ

`ECQ` is a library based on [FLINT](https://github.com/flintlib/flint/) to do arithmetic with elliptic curves over $\mathbb Q$. In particular, computing rational points on the curve.

## $2$-descent

Complete $2$-descent is work in progress and is the first goal of this project.

### Basic setup

- [x] Basic data structure
- [ ] Minimal model (expect at $p = 2$ to keep $a_1 = a_3 = 0$)
- [x] Computation of the possible value for $b_i$
- [x] The $2$-descent itself for full $2$-torsion

### Ternary quadratic forms

- [x] Reduction of TQF
- [x] Finding primitive solution with LLL
- [x] Parameterization of the solution space

### Quartics

- [ ] Equivalence of quartics
- [ ] Minimal model and reduction

### Hyper-elliptic curve of genus $2, 3$

- [x] Local solubility criterion
- [x] Basic enumeration
- [ ] Sieving for rational points
- [ ] Higher order descent

### General $2$-descent

- [x] Finding quartics in full $2$-torsion case
- [ ] Finding quartics in $2$-isogeny case
- [ ] Finding quartics without $2$-torsion

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

## Performance

For performance we consider the following curve :

$$y^2 = (x - 7265)(x - 649)(x + 7557)$$

which as of last build find the point `(-6591312805886080952551/905772661979601025 : 154861252810623357814342127836368/862042768525758359255917375)` in :

```bash
real    2m32,583s
user    2m32,411s
sys     0m0,020s
```
