# PC-bSSFP-2Compartment
#Purpose: Simulation code for the description of PC-bSSFP signal profiles of two-compartment systems for different sign conventions.  

#Background: This MATLAB code can be used for the simulation of PC-bSSFP profiles in two-compartment (singlet) systems. 
Originally the code was developed for objective falsification of opposite and aligned phase sign models in two-compartment systems.
A preliminary manuscript can be found on ArXiV https://arxiv.org/abs/2302.12548 .
This manuscript currently gets restructured and if it gets accepted for publication, I will update the link to the 
new published version. 

#Codes: \
1) Literature based simulations:\
i) "S_Simulation_twoCompartmens_PCbSSFP.m": \
Contains the main code. By chaning the parameter "Is_Opposite" the effect of different phase sign parameterizations can be investigated in two-compartment systems. These different parameterizations are all based on different published literature models. The investigation of these models with different predictions is the core and main message of this code. 
Basically this code only implements physical values like T1,T2,... and applies superposition principle "Stot=S1+S2" 
of complex values for the respective single compartment signals. The single compartment signals are simulated in the remaining codes. \
ii) "S_bSSFP_AlignedSign.m": \
Contains the aligned phase-sign parameterization based on published literaturer models.\
iii) "S_bSSFP_OppositeSign.m":\  
Contains the anti-aligned phase-sign parameterization based on published literaturer models.\
2) More general descriptions: \
i) "S_Simulation_twoCompartment_PCbSSFP_parameterized.m": \
Is similar to 1) but here all the phase signs of a single compartment PC-bSSFP profile can be parameterized in order to investigate all different combinations. The sign parameterization can be tuned via the parameters sigmaTR, sigmaTE, sigmaphi and GlobalPhase. They are properly simulated in  "S_bSSFP_SignParameterization.m"
ii) "S_bSSFP_SignParameterization.m":\
Simulates all possible phase sign parameterizations for single compartment PC-bSSFP profiles which are in "S_Simulation_twoCompartment_PCbSSFP_parameterized.m" superimposed. 

Parameters and literature references are included, and explanations are provided within the code
If you have questions or comments on the code, on the theory, phantom or experiment, feel free to contact me: 

plaehn.nils@web.de

