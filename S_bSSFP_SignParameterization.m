function profile = S_bSSFP_SignParameterization(M0,T1,T2,alpha,phi,TR,TE,deltaCS,dB0,B0,sigmaTR,sigmaTE,sigmaphi,GlobalPhase)

% Description: generation of PC-bSSFP profiles for parameterized signs

% Author of function: 
% Nils MJ Plähn, Bern, Switzerland
% E-mail: nils.plaehn@students.unibe.ch
% Department of Diagnostic, Interventional and Pediatric Radiology (DIPR), Inselspital, Bern University Hospital, University of Bern, Switzerland
% Translation Imaging Center (TIC), Swiss Institute for Translational and Entrepreneurial Medicine, Bern, Switzerland

% 1) Code based on sign parameterization for all different phase sign
% parameterizations of single-compartment PC-bSSFP profiles

% 2) Used parameters: 
% M0:      polarized magnetization of the substance, i.e. PD ( 1^H proton density)
% T1:      longitudinal relaxation time
% T2:      transversal  relaxation time
% alpha:   excitation angle of RF pulse
% phi:     linear phase increment of RF excitation pulse
% TR:      repetition time of each PC-bSSFP module
% TE:      echo time of each PC-bSSFP module
% theta:   accumulated phase 
% gamma:   gyromagnetic ratio for 1H protons
% B0:      main magnetic field strength
% dB0:     B0 inhomogeneity
% deltaCS: chemical shift

% 3) following parameterization orients base onto the defintions of the Note whose link will be shared
%    as soon as it gets accepted 

        gamma = 2*pi*42.577*10^6; 

        E1    = exp(-TR./T1);
        E2    = exp(-TR./T2);
        % reference for all sign parameterizations 
        theta = -gamma*(dB0-deltaCS*B0)*TR;  
        
        a = M0.*(1-E1).*sin(alpha);
        b = 1-E1*E2^2+(E2^2-E1)*cos(alpha);
        c = 2*(E1-1)*E2*cos(alpha/2)^2;
        % Note: the first exponent theta has opposite sign to the second exponent theta 
        % i.e. ~(1-E2*exp(-1i*(theta-phi)))*exp(1i*theta*TE/TR) is
        % MINUS/PLUS in contrast to the aligned function
        profile = -1i*a/(b+c*cos(sigmaTR*theta+sigmaphi*phi))*(1-E2*exp(1i*(sigmaTR*theta+sigmaphi*phi)))*exp(-TE/T2)*exp(1i*sigmaTE*theta*TE/TR)*exp(1i*GlobalPhase); 
end
