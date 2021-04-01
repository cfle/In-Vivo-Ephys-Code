This is the personal MATLAB code repository for C.Walsh, which includes the key
functions used to generate the data presented in "Beta-bursting in the 
Retrosplenial Cortex is a neurophysiological correlate of environmental
novelty which is disrupted in a mouse model of Alzheimer’s disease pathology."
by C.Walsh et al. (WORKING TITLE)
All of the code here was written by C.Walsh from 2017 to 2021.

Some of this code requires the Chronux Toolbox (http://chronux.org/ 
"Observed Brain Dynamics", Partha Mitra and Hemant Bokil, Oxford University
Press, New York, 2008.)

Feel free to use any and all of this code for your own uses, but
When using any of this code in a paper, please, cite this GitHub repository
https://github.com/cfle/In-Vivo-Ephys-Code and the paper "Beta-bursting in the 
Retrosplenial Cortex is a neurophysiological correlate of environmental
novelty which is disrupted in a mouse model of Alzheimer’s disease pathology."
by C.Walsh et al.(WORKING TITLE)

Due to the specific file formatting, data storage, experimental procedure etc,
this code is unlikely to work as is on other data, and will likely need to be
tailored to your data of choice.

The purpose of each script is as such:

-Spectra - Multi-taper spectral analysis for a single channel (Entire 15 minute recording)
-InitialSpectra - Multi-taper spectral analysis for a single channel (First 1 minute of recording)
-BurstDetection - Beta burst detection script (algorithm outlined in Walsh et al 2021)
-PhaseAmplitudeCoupling - Phase-amplitude coupling script based on (Canolty et al 2006) (First 1 minute of recording)
-FinalPhaseAmplitudeCoupling - As before but for the final minute of each session
-BurstTriggeredOscillation - Analysis of oscillatory activity time locked to beta bursts
-SimpleMUADuringBursts - Simple MUA detection script to investigate multi-unit activity during beta bursts

Please report any issues or bugs to CWalsh at callum.f.walsh@gmail.com
