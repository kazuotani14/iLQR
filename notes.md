### To-do

* Test acrobot with constraints enforced in forward pass. Why doesn't it work? 
    * Bug: in forward_pass if we use us (clamped) instead of u_curr (unclamped), we can't solve the problem. Should we even have to clamp a lot? Isn't boxQP supposed to take care of this? See lines for clamping in `forward_pass()`
        * K matrix in direction of clamped inputs should be zero. TODO check this
        * limits are crossed without applying control gains (just feedforward), but only by a bit - TODO check K matrix

* Fix calculate_cxu - why do we calculate the final one? 

* Optimize
    * run in MPC mode on perturbed model - how fast is it with warm-start?
    * speed up derivatives - openmp currently turned off because it was leading to some non-deterministic (?) behavior. See `#pragma` lines in `derivatives.cpp`
    * try callgrind (on linux; having some issues on laptop) 
    * replace virtual function call to dynamics with template/policy
    * Speed up backward pass - what else can be done? 
    * parallelize line search with openmp? how often does it end on the first few iterations?
    * Turn off assertions in a clean way

* Implement other methods, compare

* Save notes of derivations somewhere
* Read extension papers - parallel (ethz thesis)
* Solidify understanding - read thru notes and lectures again 

* clean up code and formatting
* python script for plotting output - finalize with commandline arguments
* additional documentation for future me

---

* Try a more complex model - mujoco? 
* debuggers: gdb, rr, sublime integration
* remove error codes - use enums

### Algorithm

##### Questions 

* Connection between DDP and backprop
* Understand pros/cons vs other methods. 

##### Pros/cons vs other trajectory optimization algorithms

* iLQR/DDP searches in space of control trajectories; it solves an $m$-dimensional problem $N$ times instead of solving a single problem of size $mN$ (as a direct method would). So complexity is $O(Nm^3)$ vs. $O(N^3m^3)$.
    * Note that this ignores the cost of calculating the derivatives - this takes up the bulk of the time in iLQR if finite-differences is used. Ideally, we would have analytic derivatives of dynamics model.
    * This also ignores methods for explotining sparsity in direct methods
* Intermediate results of shooting methods (of which iLQR is one) are always strictly feasible and "dynamics constraints" are unnecessary.
* Direct methods tend to find optima easier than indirect/shooting methods. Shooting methods are much more sensitive to local minima.
* Shooting methods make it easy to use warm-starts.
* Another advantage of iLQR is that it finds a set of local linear controllers around the desired trajectory, which can be used to keep the trajectory from drifting

From Tassa 2014: "The direct approach discards the temporal structure and is forced to search in a constrained space which is slower, however it is far easier to find better optima through continuation. The indirect approach is faster and better suited for warm-starting, but is far more sensitive to local minima."

_Misc_

* For iLQR, undesired states are usually encoded as repulsive "potential fields" in the cost function. We can't do this for the control input limits, which is why the boxQP is necessary (the "control-limited" part of the paper title)

### Misc notes

* `cmake -DCMAKE_BUILD_TYPE=Debug ..` to compile with debug flags
* Find TODO:
```
grep -nr TODO src
grep -n TODO include/*
```

### References

* [Good explanation of iLQR](https://studywolf.wordpress.com/2016/02/03/the-iterative-linear-quadratic-regulator-method/)
* http://rll.berkeley.edu/deeprlcourse/f17docs/lecture_8_model_based_planning.pdf
    * Connection between iLQR/DDP and Newton's method
* http://underactuated.csail.mit.edu/underactuated.html?chapter=lqr
* http://underactuated.csail.mit.edu/underactuated.html?chapter=dp
    * Connection between iLQR and value iteration, HJB equation
* https://en.wikipedia.org/wiki/Levenberg%E2%80%93Marquardt_algorithm
* https://people.eecs.berkeley.edu/~pabbeel/cs287-fa15/slides/lecture5-LQR.pdf

##### Other implementations

* [Yuval Tassa's Matlab code](https://www.mathworks.com/matlabcentral/fileexchange/52069-ilqg-ddp-trajectory-optimization)
* [DDP-Generator](https://github.com/jgeisler0303/DDP-Generator) - fastest implementation I've found. Generates problem-specific code with analytic gradients for dynamics and cost functions.
* [Ben Stephen's implementation](http://www.cs.cmu.edu/~bstephe1/)
* [Tglad's ILQR implementation](https://github.com/TGlad/ILQR)
* [Arun Venkatraman's implementation](https://github.com/LAIRLAB/qr_trees/blob/master/src/ilqr/iLQR.cc)

##### Misc 

* [Eigen FAQ](http://eigen.tuxfamily.org/index.php?title=FAQ)
