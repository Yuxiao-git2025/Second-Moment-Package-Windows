1. Measurements are intermediate variables, created and written to this variable 
group each time. An example for analysis is provided in the <Input> folder.
Don't worry, once all station waveforms have been fitted, the Measurements variable 
will be complete. You may choose to terminate the waveform fitting at any time, and 
our program will synchronously update the completed process

2. toppinputs and toppoutputs are programs within the hypodd package used to calculate 
takeoff angles. We have slightly modified them to make them more convenient for Windows
systems, and you may need to prepare gfortran and gcc in advance (created via MinGW). 
Of course, if your system is Linux, you can directly use make to recompile and 
obtain the topp program.

3.File Info.txt contains the sampling results for jackknife and bootstrap procedure

Last updated, XA,Yu 2026/10/03