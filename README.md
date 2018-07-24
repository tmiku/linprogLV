# SME gravitational sector coefficient constraints

This code is meant to refine coefficients for Lorentz violation as described in the Standard Model Extension (SME) using differences between arrival times of signals from gravitational waves and their EM counterparts. The codebase's goals, in order, are:

* Calculate the coefficients given a single event, its sky location, and its maximimum and minimum differences in velocity.
* Update a list of these events and the current best set of coefficients.
* Use linear programming to optimize the coefficients further using a set of nine or more events.
* Adapt the program to save and use luminosity distance and travel time to calculate velocity differences automatically.
* Search for differences between time difference at source and at detection using distance from event.
