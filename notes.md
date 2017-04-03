### Algorithm

* For ilqr, undesired states are usually encoded as repulsive "potential fields" in the cost function. We can't do this for the control input limits, which is why the boxQP is necessary (the "control-limited" part of the paper title)

_Pros/cons vs other trajectory optimization algorithms_

* iLQR/DDP searches in space of control trajectories; it solves an $m$-dimensional problem $N$ times instead of solving a single problem of size $mN$ (as a direct method would). So complexity is $O(Nm^3)$ vs. $O(N^3m^3)$.
    * Note that this ignores the cost of calculating the derivatives - this takes up the bulk of the time in iLQR if finite-differences is used. Ideally, we would have analytic derivatives of dynamics model.
* Intermediate results of shooting methods (of which iLQR is one) are always strictly feasible and "dynamics constraints" are unnecessary.
* Direct methods tend to find optima easier than indirect/shooting methods. Shooting methods are much more sensitive to local minima.
* Shooting methods make it easy to use warm-starts.

### To-do

* add Model with nonlinear dynamics
    * drifting car
    * two-link arm
* remove error codes - use enums
* python script for plotting output

_Second priority_

* integrate cost for change in control input
* add documentation

### Usage

* `mkdir build; cd build`
* `cmake ..` or `cmake -DCMAKE_BUILD_TYPE=Debug ..` to compile with debug flags
* `make`

### Misc references

* [Eigen FAQ](http://eigen.tuxfamily.org/index.php?title=FAQ)
