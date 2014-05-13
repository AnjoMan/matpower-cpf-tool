# Matpower CPF Customization #

This toolbox contains a custom version of the CPF toolbox for [MatPower](http://www.pserc.cornell.edu/matpower/), a Matlab - based power system simulation toolbox. I am adding features to this toolbox in order to evaluate solvability of systems.

Changes so far:
* Implemented participation factors.
* Extended to scaling of all loads instead of just one load.
* Implemented Lagrange Polynomial interpolation for prediction step to improve accuracy and speed of predictions.
* Implemented variable step sizes to improve speed for relatively flat sections.