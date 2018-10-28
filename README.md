# Time offsets in MIMO OFDM Systems 
This repository includes the codes used for the FYP titled: Time offsets in MIMO OFDM Systems

The final versions of the simulations are named below and can be run using MATLAB 2017b.

results_SISO_v1.m

results_MISO_v1.m

results_MIMO_v1.m

results_MIMO_OFDMA_v1.m

Inside the folder **all code** are all the versions and edits used to create the scripts above. To test the different systems open MATLAB and run the scripts. Different time offsets can be changed inside the simulation to get different results and to view the sensitivity of the system to these time offsets

It was found that the cyclic prefix in these systems were able to absorb the time offset up to the length of the cyclic prefix. A time offset larger than the set cyclic prefix would reduce the performance of the system and introduce errors.

