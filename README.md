Phase-cycled bSSFP in two-compartment systems
- Contains Matlab code to simulate phase-cycled bSSFP profiles of two-compartment systems using different sign conventions.  

#Background: The code was developed to emphasize the difference in opposite and aligned phase sign models in two-compartment systems.
The preprint of the manuscript can be found on ArXiV https://arxiv.org/abs/2302.12548 .

#Codes: \
1) Simulation code based on references in the literature :\
i) "S_Simulation_twoCompartmens_PCbSSFP.m": \
Contains the main code. By changing the parameter "Is_Opposite" the effect of different phase sign parameterizations can be investigated in two-compartment systems. These different parameterizations are based on different published works. The comparison of these models with different descriptions is performed using this code. 
The code implements physical values like T1,T2,... and applies superposition principle "Stot=S1+S2" 
of complex values for the respective single compartment signals, which are simulated in the remaining provided code. \
ii) "S_bSSFP_AlignedSign.m": \
Contains the aligned phase-sign parameterization based on published work.\
iii) "S_bSSFP_OppositeSign.m":\  
Contains the anti-aligned phase-sign parameterization based on published work .\
2) More general descriptions: \
i) "S_Simulation_twoCompartment_PCbSSFP_parameterized.m": \
Is similar to 1) but here the phase signs of a single compartment PC-bSSFP profile can be parameterized in order to investigate all different combinations. The sign parameterization can be tuned via the parameters sigmaTR, sigmaTE, sigmaphi and GlobalPhase. They are properly simulated in  "S_bSSFP_SignParameterization.m" \
ii) "S_bSSFP_SignParameterization.m":\
Simulates all possible phase sign parameterizations for single compartment PC-bSSFP profiles which are in "S_Simulation_twoCompartment_PCbSSFP_parameterized.m" superimposed. 

Parameters and  references to literature are included and explanations are provided within the code.
If you have questions or comments on the code, on the theory, phantom or experiment, please contact plaehn.nils@web.de

