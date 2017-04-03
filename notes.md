### TO-DO

* read paper, understand algorithm
* Test inputs/output:
	* compute_derivatives
    * backward_pass

* step 2 in iLQR_core
	* check for termination due to small gradient
		* the way they do it in matlab looks complicated - is there an easier way?

* remove error codes - use something else

* integrate cost for change in control input
* make one file that runs all unit tests

* test whole thing and see if output looks reasonable - plot
	* write final control sequence and states to csv file
	* python script for plotting

_Make stuff faster_

* Check for most expensive step - use chrono

1. calculate_derivatives : this is taking most of the time
	* store auxiliary variables
	* change incrementing method?
	* take analytic derivatives from DDP-generator

2. backward_pass

3. forward_pass
	* vectorize for parallel line search?

_Later_

* add docs, explanation of algorithm?
* move to using Eigen::Ref
* helper functions are all inline right now - move these to separate file?
* Apparently std::function has some overhead - test alternatives

### Misc

To build

* In Build directory, 'cmake ..'
* In same directory, 'make'
