timeWalk_photodiode_perPMT.cpp
------------------------------
Time-walk corrections by fitting the standard y = A/Sqrt(x) + B to a plot of:
tdc for a given PMT - reference vs. ADC
where reference is the laser photodiode.
The photodiode was plugged in the space that was occupied by V16-A
(which corresponds to sector: 3, layer: 6, component: 6)

These codes were written specifically for run 261.
Run 260 would in principle also work, but ranges will need adjusting.

test_parameters.cpp
------------------------------
Code to check the quality of the timewalk parameters extracted with the
previous code.
