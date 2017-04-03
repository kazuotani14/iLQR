### TO-DO

* read paper, understand algorithm
* remove error codes - use enums?

* make one file that runs all unit tests
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

* namespace
* change layout of standard includes
* helper functions are all inline right now - move these to separate file?
* move to using Eigen::Ref
* Apparently std::function has some overhead - test alternatives
* integrate cost for change in control input

* Test for back-pass, based on linear system?
* add docs, short explanation of algorithm?

### Misc

To build

* In Build directory, 'cmake ..'
* In same directory, 'make'
