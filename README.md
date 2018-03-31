## iLQR

Implementation of iLQR (Iterative Linear Quadratic Regulator) algorithm for trajectory optimization.

See `notes.md` for comments on algorithm, implementation, etc. 

### Usage

* `mkdir build; cd build`
* `cmake ..` or `cmake -DCMAKE_BUILD_TYPE=Debug ..` to compile with debug flags
* `make`
* Define a dynamics and cost model based on Model (see double_integrator example)

### Other implementations

* [Yuval Tassa's Matlab code](https://www.mathworks.com/matlabcentral/fileexchange/52069-ilqg-ddp-trajectory-optimization)
* [DDP-Generator](https://github.com/jgeisler0303/DDP-Generator) - best iLQR implementation I've found. Generates problem-specific code with analytic gradients for dynamics and cost functions.
* [Tglad's ILQR implementation](https://github.com/TGlad/ILQR)

### Papers

* Tassa, Yuval, Nicolas Mansard, and Emo Todorov. "Control-limited differential dynamic programming." Robotics and Automation (ICRA), 2014 IEEE International Conference on. IEEE, 2014.
* Li, Weiwei, and Emanuel Todorov. "Iterative linear quadratic regulator design for nonlinear biological movement systems." ICINCO (1). 2004.
