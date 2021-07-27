# L-BFGS-B

## Software for Large-scale Bound-constrained Optimization
[`L-BFGS-B`](http://users.iems.northwestern.edu/~nocedal/lbfgsb.html) is a limited-memory quasi-Newton code for bound-constrained optimization,
i.e., for problems where the only constraints are of the form `l <= x <= u`.
It is intended for problems in which information on the Hessian matrix is difficult to obtain,
or for large dense problems.
`L-BFGS-B` can also be used for unconstrained problems, and in this case performs similarly to its predecessor,
algorithm [`L-BFGS`](http://users.iems.northwestern.edu/~nocedal/lbfgs.html) (Harwell routine VA15). The algorithm is implemented in Fortran 77.

## Authors

* [Ciyou Zhu](http://web.archive.org/web/19990129014554/http://www.ece.nwu.edu/%7Eciyou/)
* Richard Byrd
* [Jorge Nocedal](http://www.ece.northwestern.edu/~nocedal)
* [Jose Luis Morales](http://web.archive.org/web/20090903033420/http://www.ece.northwestern.edu:80/~morales)

## Related Publications

* R. H. Byrd, P. Lu, J. Nocedal and C. Zhu. [A Limited Memory Algorithm for Bound Constrained Optimization](https://doi.org/10.1137/0916069) (1995), SIAM Journal on Scientific and Statistical Computing, Vol. 16, Num. 5, pp. 1190-1208
* C. Zhu, R. H. Byrd and J. Nocedal. [L-BFGS-B: Algorithm 778: L-BFGS-B, FORTRAN routines for large scale bound constrained optimization](https://doi.org/10.1145/279232.279236) (1997), ACM Transactions on Mathematical Software, Vol. 23, Num. 4, pp. 550-560
* J.L. Morales and J. Nocedal. [L-BFGS-B: Remark on Algorithm 778: L-BFGS-B, FORTRAN routines for large scale bound constrained optimization](https://doi.org/10.1145/2049662.2049669) (2011), ACM Transactions on Mathematical Software, Vol. 38, Num. 1
* R. H. Byrd, J. Nocedal and R. B. Schnabel. [Representations of quasi-Newton matrices and their use in limited memory methods](https://doi.org/10.1007/BF01582063) (1994), Mathematical Programming, Vol. 63, pp. 129-156

Note that the subspace minimization in the [LBFGSpp](https://github.com/yixuan/LBFGSpp) implementation
is an exact minimization subject to the bounds, based on the BOXCQP algorithm:
* C. Voglis and I. E. Lagaris, [BOXCQP: An Algorithm for Bound Constrained Convex Quadratic Problems](http://www.cs.uoi.gr/~voglis/boxcqp.pdf) (2004), 1st International Conference "From Scientific Computing to Computational Engineering", Athens, Greece

For an eagle-eye overview of `L-BFGS-B` and the genealogy `BFGS`->`L-BFGS`->`L-BFGS-B`,
see [Henao's Master's thesis](https://cs.nyu.edu/overton/mstheses/henao/msthesis.pdf).

## Related Software

* [wilmerhenao/L-BFGS-B-NS](https://github.com/wilmerhenao/L-BFGS-B-NS): An L-BFGS-B-NS Optimizer for Non-Smooth Functions
* [pcarbo/lbfgsb-matlab](https://github.com/pcarbo/lbfgsb-matlab): A MATLAB interface for L-BFGS-B
* [bgranzow/L-BFGS-B](https://github.com/bgranzow/L-BFGS-B): A pure Matlab implementation of L-BFGS-B (LBFGSB)
* [constantino-garcia/lbfgsb_cpp_wrapper](https://github.com/constantino-garcia/lbfgsb_cpp_wrapper): A simple C++ wrapper around the original Fortran L-BGSG-B routine
* [yixuan/LBFGSpp](https://github.com/yixuan/LBFGSpp): A header-only C++ library for L-BFGS and L-BFGS-B algorithms
* [chokkan/liblbfgs](https://github.com/chokkan/liblbfgs): libLBFGS: a library of Limited-memory Broyden-Fletcher-Goldfarb-Shanno (L-BFGS)
* [mkobos/lbfgsb_wrapper](https://github.com/mkobos/lbfgsb_wrapper): Java wrapper for the Fortran L-BFGS-B algorithm
* [yuhonglin/Lbfgsb.jl](https://github.com/yuhonglin/Lbfgsb.jl): A Julia wrapper of the l-bfgs-b fortran library
* [Gnimuc/LBFGSB.jl](https://github.com/Gnimuc/LBFGSB.jl): Julia wrapper for L-BFGS-B Nonlinear Optimization Code
* [afbarnard/go-lbfgsb](https://github.com/afbarnard/go-lbfgsb): L-BFGS-B optimization for Go, C, Fortran 2003
* [nepluno/lbfgsb-gpu](https://github.com/nepluno/lbfgsb-gpu): An open source library for the GPU-implementation of L-BFGS-B algorithm
* [Chris00/L-BFGS-ocaml](https://github.com/Chris00/L-BFGS-ocaml):  OCaml bindings for L-BFGS
* [dwicke/L-BFGS-B-Lua](https://github.com/dwicke/L-BFGS-B-Lua): L-BFGS-B lua wrapper around a L-BFGS-B C implementation
* [avieira/python_lbfgsb](https://github.com/avieira/python_lbfgsb): Pure Python-based L-BFGS-B implementation
* [ybyygu/rust-lbfgsb](https://github.com/ybyygu/rust-lbfgsb): Ergonomic bindings to L-BFGS-B code for Rust
* [rforge/lbfgsb3c](https://rdrr.io/rforge/lbfgsb3c): Limited Memory BFGS Minimizer with Bounds on Parameters with optim() 'C' Interface for R
* [florafauna/optimParallel-python](https://github.com/florafauna/optimParallel-python): A parallel version of `scipy.optimize.minimize(method='L-BFGS-B')`

## Notes on this repository

[I (J. Schilling)](https://github.com/jonathanschilling/) took the freedom to

* put the L-BFGS-B code obtained from [the original website](http://users.iems.northwestern.edu/~nocedal/Software/Lbfgsb.3.0.tar.gz) up in this repository,
* divide the subroutines and functions into separate files,
* convert parts of the documentation into a format understandable to [doxygen](https://www.doxygen.nl/index.html) and
* replace the included BLAS/LINPACK routines with calls to user-provided BLAS/LAPACK routines.
  To be precise, the calls to LINPACK's `dtrsl` were replace with calls to LAPACK's `dtrsm`
  and the calls to LINPACK's `dpofa` were replaces with calls to LAPACK's `dpotrf`.

## Building

A CMake setup is provided for L-BFGS-B in this repository.
External modules for `BLAS` and `LAPACK` have to be installed on your system.
Then, building the shared library `liblbgsb.so` and the examples `driver*.f` and `driver*.f90` works as follows:

```bash
> mkdir build
> cd build
> cmake ..
> make -j
```

The resulting shared library and the driver executables can be found in the `build` directory.

## Concluding remarks

The current release is version 3.0. The distribution file was last changed on 02/08/11.

This work was in no way intending to infringe any copyrights or take credit for others' work.
Feel free to contact me at any time in case you noticed something against the rules.
Above documentation is obtained from the archived version of the [original manual](http://web.archive.org/web/19991005125105/http://www.ece.nwu.edu:80/%7Eciyou/pp9/pp9.html).

A PDF version of the documentation is available here: [L-BFGS-B.pdf](https://github.com/jonathanschilling/L-BFGS-B/blob/gh-pages/L-BFGS-B.pdf).
