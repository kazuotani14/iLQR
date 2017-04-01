### TO-DO

* Test individual functions
* Separate dynamics - make dynamics base class
	* make point mass dynamics

* Test inputs/output:
	* costs (also look for ways to make this faster)
	* compute_derivatives
    * forward_pass
    * backward_pass

* step 2 in iLQR_core
	* check for termination due to small gradient
		* the way they do it in matlab looks complicated - is there an easier way?

- integrate cost for change in control input

* Replace timing with chrono
* Apparently std::function has some overhead - test alternatives
* move to using Eigen::Ref

- write final control sequence and states to text file
	- python script for plotting
- test whole thing and see if output looks reasonable - plot in matlab

_Make stuff faster_

1. calculate_derivatives : this is taking most of the time
	* store auxiliary variables
	* change incrementing method?
	* take analytic derivatives from DDP-generator

2. backward_pass

3. forward_pass
	* vectorize for parallel line search?

_Later_

* fix all the messy stuff
	* extract methods into smaller methods, especially in fw, bw, boxQP

* make it more readable
	* refactor to completely separate iLQR class from Model class
		* make sure everything in iLQR is dimension-agnostic (there are currently some hacks for car model)
	* make skeleton for new Model class, for reference

* add docs
* helper functions are all inline right now - move these to cpp file

### Misc

Files

* iLQR_* : ideally for the algorithm, agnostic to dynamics and cost
	* currently, there are some hacks that hard-code num_controls = 2
* loco_car_* : our specific dynamics model and costs
* main : where optimization is run. can also call individual tester functions to
	verify functionality

To build

* Make sure you have eigen installed
* In Build directory, 'cmake ..'
* In same directory, 'make'
