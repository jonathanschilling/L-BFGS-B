# L-BFGS-B

## Software for Large-scale Bound-constrained Optimization
[`L-BFGS-B`](http://users.iems.northwestern.edu/~nocedal/lbfgsb.html) is a limited-memory quasi-Newton code for bound-constrained optimization,
i.e., for problems where the only constraints are of the form `l <= x <= u`.
It is intended for problems in which information on the Hessian matrix is difficult to obtain,
or for large dense problems.
`L-BFGS-B` can also be used for unconstrained problems, and in this case performs similarly to its predecessor,
algorithm [`L-BFGS`](http://users.iems.northwestern.edu/~nocedal/lbfgs.html) (Harwell routine VA15). The algorithm is implemented in Fortran 77.

The current release is version 3.0. The distribution file was last changed on 02/08/11.

## Authors

* [Ciyou Zhu](http://web.archive.org/web/19990129014554/http://www.ece.nwu.edu/%7Eciyou/)
* Richard Byrd
* [Jorge Nocedal](http://www.ece.northwestern.edu/~nocedal)
* [Jose Luis Morales](http://web.archive.org/web/20080509084403/http://www.ece.northwestern.edu:80/~morales/)

## Relevant Publications

* R. H. Byrd, P. Lu, J. Nocedal and C. Zhu. [A Limited Memory Algorithm for Bound Constrained Optimization](https://doi.org/10.1137/0916069), (1995), SIAM Journal on Scientific and Statistical Computing , 16, 5, pp. 1190-1208
* C. Zhu, R. H. Byrd and J. Nocedal. [L-BFGS-B: Algorithm 778: L-BFGS-B, FORTRAN routines for large scale bound constrained optimization](https://doi.org/10.1145/279232.279236) (1997), ACM Transactions on Mathematical Software, Vol 23, Num. 4, pp. 550 - 560
* J.L. Morales and J. Nocedal. [L-BFGS-B: Remark on Algorithm 778: L-BFGS-B, FORTRAN routines for large scale bound constrained optimization](https://doi.org/10.1145/2049662.2049669) (2011), ACM Transactions on Mathematical Software, Vol. 38, Num. 1

## Notes on this repository

[I](https://github.com/jonathanschilling/) took the freedom to

* put the L-BFGS-B code obtained from [the original website](http://users.iems.northwestern.edu/~nocedal/Software/Lbfgsb.3.0.tar.gz) up in this repository,
* divide the subroutines and functions into separate files,
* convert parts of the documentation into a format understandable to [doxygen](https://www.doxygen.nl/index.html) and
* adjust the `Makefile` to accomodate the separate files and additionally generate a statically linked `liblbfgsb.a` library.

This work was in no way intending to infringe any copyrights or take credit for others' work.
Feel free to contact me at any time in case you noticed something against the rules.
