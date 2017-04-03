### TO-DO

* remove error codes from matlab - use enums
* python script for plotting output
* make calculate_derivatives faster: takes >50% of execution time
	* multi-thread?
* Allow warmstart - change interface

_Second priority_

* integrate cost for change in control input
* add docs, short explanation of algorithm?
* move to using Eigen::Ref
* Apparently std::function has some overhead - test alternatives
* Test for back-pass, based on linear system?

### Misc

To build

* In Build directory, 'cmake ..'
* In same directory, 'make'
