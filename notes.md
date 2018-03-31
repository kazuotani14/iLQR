### To-do

* Note: get_cost_2nd_derivative_mt doesn't work anymore because there's a (t<T) check in functions
* python script for plotting output - finalize with commandline arguments
* clean up code and formatting

Find TODO:
```
grep -nr TODO .
```

* Optimize
    * Debuggers: gdb, rr, sublime integration
    * replace clamp with better implementation
    * check data structures
    * add openmp 
        * calculating derivatives - use finite_diff2?
        * parallel line search

---


* remove error codes - use enums
* add general features:
    * cost for change in control input
    * repulsor from undesired states (bounds, obstacles) - log or repulsive field?
* additional documentation for future me

### Algorithm

##### Questions 

* Connection between DDP and backprop
* Where does const term come from?
    * How is infinite-horizon term derived?
* Understand pros/cons vs other methods. 

##### Pros/cons vs other trajectory optimization algorithms

* iLQR/DDP searches in space of control trajectories; it solves an $m$-dimensional problem $N$ times instead of solving a single problem of size $mN$ (as a direct method would). So complexity is $O(Nm^3)$ vs. $O(N^3m^3)$.
    * Note that this ignores the cost of calculating the derivatives - this takes up the bulk of the time in iLQR if finite-differences is used. Ideally, we would have analytic derivatives of dynamics model.
    * This also ignores sparsity in direct methods
* Intermediate results of shooting methods (of which iLQR is one) are always strictly feasible and "dynamics constraints" are unnecessary.
* Direct methods tend to find optima easier than indirect/shooting methods. Shooting methods are much more sensitive to local minima.
* Shooting methods make it easy to use warm-starts.
* Another advantage of iLQR is that it finds a set of local linear controllers around the desired trajectory, which can be used to keep the trajectory from drifting

_Misc_

* For iLQR, undesired states are usually encoded as repulsive "potential fields" in the cost function. We can't do this for the control input limits, which is why the boxQP is necessary (the "control-limited" part of the paper title)

### References

* [Good explanation of iLQR](https://studywolf.wordpress.com/2016/02/03/the-iterative-linear-quadratic-regulator-method/)
* http://rll.berkeley.edu/deeprlcourse/f17docs/lecture_8_model_based_planning.pdf
    * Connection between iLQR/DDP and Newton's method
* http://underactuated.csail.mit.edu/underactuated.html?chapter=lqr
* http://underactuated.csail.mit.edu/underactuated.html?chapter=dp
    * Connection between iLQR and value iteration, HJB equation
* https://en.wikipedia.org/wiki/Levenberg%E2%80%93Marquardt_algorithm
* https://people.eecs.berkeley.edu/~pabbeel/cs287-fa15/slides/lecture5-LQR.pdf

##### Misc 

* [Eigen FAQ](http://eigen.tuxfamily.org/index.php?title=FAQ)
