### To-Do

* add Model with nonlinear dynamics
* remove error codes - use enums?
* Allow warmstart - change interface
* Look for other ways to make this faster
* python script for plotting output

_Second priority_

* integrate cost for change in control input
* add docs, short explanation of algorithm?
* move to using Eigen::Ref
* Apparently std::function has some overhead - test alternatives
* Test for back-pass, based on linear system?

### Algorithm

_Pros/cons vs other trajectory optimization algorithms_

* iLQR/DDP searches in space of control trajectories; it solves an $m$-dimensional problem $N$ times instead of solving a single problem of size $mN$ (as a direct method would). So complexity is $O(Nm^3)$ vs. $O(N^3m^3)$.
    * Note that this ignores the cost of calculating the derivatives - this takes up the bulk of the time in iLQR if finite-differences is used. Ideally, we would have analytic derivatives of dynamics model.
* Intermediate results of shooting methods (of which iLQR is one) are always strictly feasible and "dynamics constraints" are unnecessary.
* Direct methods tend to find optima easier than indirect/shooting methods. Shooting methods are much more sensitive to local minima.
* Shooting methods make it easy to use warm-starts.

### Misc

To build

* In Build directory, 'cmake ..'
* In same directory, 'make'
