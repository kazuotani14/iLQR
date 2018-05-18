## iLQR

Implementation of iLQR (Iterative Linear Quadratic Regulator) algorithm for trajectory optimization, based on [Yuval Tassa's Matlab implementation](https://www.mathworks.com/matlabcentral/fileexchange/52069-ilqg-ddp-trajectory-optimization).

See `notes.md` for comments on algorithm, implementation, to-do's, etc. 

### Current status

Mostly working, except the "control-limited" part... See notes for details. If you're looking for a C++ implementation with control limits, the best one I've found is [Ben Stephen's implementation](http://www.cs.cmu.edu/~bstephe1/) (I haven't tried it though). 

### Usage

* Install eigen or download into `include` (https://eigen.tuxfamily.org/index.php?title=Main_Page#Download)
* `mkdir build; cd build`
* `cmake ..` 
* `make`
* Define a dynamics and cost model based on Model (see double_integrator example), or run with `./run_iLQR acrobot`

### Papers

* Tassa, Yuval, Nicolas Mansard, and Emo Todorov. "Control-limited differential dynamic programming." Robotics and Automation (ICRA), 2014 IEEE International Conference on. IEEE, 2014.
* Li, Weiwei, and Emanuel Todorov. "Iterative linear quadratic regulator design for nonlinear biological movement systems." ICINCO (1). 2004.
