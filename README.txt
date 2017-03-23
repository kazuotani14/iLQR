To build
- Make sure you have eigen installed
- In Build directory, 'cmake ..'
- In same directory, 'make'

Files
- iLQR_* : ideally for the algorithm, agnostic to dynamics and cost
	- currently, there are some hacks that hard-code num_controls = 2
- loco_car_* : our specific dynamics model and costs
- main : where optimization is run. can also call individual tester functions to
	verify functionality
