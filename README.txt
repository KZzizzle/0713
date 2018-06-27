This directory contains the code used for data analysis for our paper, "Shared proportional control of a dexterous myoelectric prosthesis." Below is a short description of the purpose of each file. In addition, each file is commented for more readability.


PrelimScriptVMGLabjack is the main script for analyzing MLP predictions.
The script can be used for the data collected either from Experiment 2 or  3. To analyze data, you must provide the cell array of "sessions" which corresponds to a directory of the experimental session. 

svm_prelimscript is called by PrelimScriptVMGLabjack and is used for some state-decoding analysis.

ScriptVMGLabjack_glovesessions is equivalent to PrelimScriptVMGLabjack except that it is only used for data with the kinematic glove. 

parsebag is used to parse bag data obtained from Gazebo simulation in order to analyze trials of grasping and releasing in Experiment 3. 

saveBagInterp is used to align the .bag and the data from the VC++ output for analysis. Also puts both data streams into the same time scale. 

analyzebag_strict does most of the analysis for Experiment 3 including counting contacts, calculating hold times. 

ComplianceEMG is called by analyzebag_strict in order to analyze the EMG activity during grasp and release trials. 