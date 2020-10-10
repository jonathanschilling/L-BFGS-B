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
The following documentation is obtained from the archived version of the [original manual](http://web.archive.org/web/19991005125105/http://www.ece.nwu.edu:80/%7Eciyou/pp9/pp9.html).

## Introduction

The purpose of algorithm L-BFGS-B is to minimize a nonlinear function of `n` variables,

```
min f(x)
```

subject to the simple bounds

```
l <= x <= u ,
```

where the vectors `l` and `u` represent lower and upper bounds on the variables.
Not all the variables need to have bounds;
in fact the algorithm is also appropriate and efficient for solving unconstrained problems.
The user must supply the gradient `g`, but knowledge about the Hessian matrix of `f` is not required.
For this reason the algorithm can be useful for solving large problems
in which the Hessian matrix is difficult to compute or is dense.

The algorithm is described in detail in [8], and proceeds roughly as follows.
At each iteration a limited memory BFGS approximation to the Hessian is updated.
This limited memory matrix is used to define a quadratic model of the objective function `f`.
A search direction is then computed using a two-stage approach:
first, the gradient projection method [15],[3],[18],[9] is used to identify a set of active variables,
i.e., variables that will be held at their bounds;
then the quadratic model is approximately minimized with respect to the free variables.
The search direction is defined to be the vector leading from the current iterate to this approximate minimizer.
Finally a line search is performed along the search direction using the subroutine described in [17].
A novel feature of the algorithm is that the limited memory BFGS matrices
are represented in a compact form [7] that is efficient for bound constrained problems.

The user can control the amount of storage required by L-BFGS-B
by selecting a parameter `m` that determines the number of BFGS corrections saved.
The algorithm requires roughly `(12+2m)n` storage locations,
and since small values of `m` (say `3 <= m <= 20`) are recommended,
it can be used to solve very large problems.
The computational cost of one iteration of the algorithm is modest,
ranging from `(4 m + 1) n` multiplications when no bounds are active,
to approximately `m^2 n` multiplications when all variables are at their bounds.

If no bounds are active at the solution,
it is appropriate to stop the iteration when the norm of the gradient `g` is sufficiently small.
The corresponding quantity for the case when some bounds are active
is the norm of the projected gradient, which we denote by `||proj g||`,
and which is defined, for example, in [8].
Both the output of L-BFGS-B and its documentation, make reference to the projected gradient.

L-BFGS-B is an extension of the limited memory algorithm (L-BFGS) for unconstrained optimization described in [16] and implemented as Harwell routine VA15 [12]. The main improvement is the ability of L-BFGS-B to deal with bounds on the variables. Even though this requirement makes the new algorithm far more complex than its predecessor, the two codes perform similarly on unconstrained problems. Therefore L-BFGS-B could be considered to supersede L-BFGS - except for one fact that can be important in some applications: L-BFGS-B requires 8 more -vectors of storage.

L-BFGS-B is, at present, the only limited memory quasi-Newton algorithm capable of handling bounds on the variables; other published codes [5], [6], [13], [16] are only able to solve unconstrained problems. We note also that the nonlinear conjugate gradient method [14], which is used for solving many large unconstrained problems, has not been adequately extended to handle bounds on the variables, and L-BFGS-B can be used in its place.

The advantages of L-BFGS-B are: (i) the code is easy to use, and the user need not supply information about the Hessian matrix or about the structure of the objective function; (ii) the storage requirements of the algorithm are modest and can be controlled by the user; (iii) the cost of the iteration is low, and is independent of the properties of the objective function. Due to this, L-BFGS-B is recommended for solving large problems in which the Hessian matrix is not sparse or is difficult to compute.

However L-BFGS-B suffers from the following drawbacks: (i) it is not rapidly convergent, and on difficult problems can take a large number of function evaluations to converge; (ii) on highly ill-conditioned problems the algorithm may fail to obtain high accuracy in the solution; (iii) the algorithm cannot make use of knowledge about the structure of the problem to accelerate convergence. 
