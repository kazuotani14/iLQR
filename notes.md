### TO-DO

* Test:
    * backward_pass - check against matlab code
	    * figure out why cu differs starting from iteration 2
	    * Finish this and we're pretty much done!

* step 2 in iLQR_core
	* check for termination due to small gradient
		* the way they do it in matlab looks complicated - is there an easier way?

* read paper, understand algorithm
* remove error codes - use enums?

* make one file that runs all unit tests
* test whole thing and see if output looks reasonable - plot
	* write final control sequence and states to csv file
	* python script for plotting

_Make stuff faster_

* const everything
* Check for most expensive step - use chrono

1. calculate_derivatives : this is prob taking most of the time
	* multi-thread?
2. backward_pass
3. forward_pass
	* vectorize for parallel line search?

* Allow warmstart

_Later_

* Apparently std::function has some overhead - test alternatives
* move to using Eigen::Ref
* helper functions are all inline right now - move these to separate file?
* integrate cost for change in control input

* add docs, short explanation of algorithm?

### Misc

To build

* In Build directory, 'cmake ..'
* In same directory, 'make'
