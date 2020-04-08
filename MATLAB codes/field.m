%% Control field generator

%% Introduction
function efield=field(amp,omega,phi,t)
% DESCRIPTION:
    % Calculates a temporal (control) field given the spectral parameters
    
% INPUT:
    % amp: spectral amplitude modes,
    % omega: spectral angular frequency,
    % phi: spectral phase,
    % t: time vector.
    
% OUTPUT:
    % efield: temporal efield vector.
    
% WRITTEN BY: 
    % Andy Koswara,
    % Advisor: Prof. Raj Chakrabarti,
    % School of Chemical Engineering,
    % Purdue University.
    
% VERSION HISTORY
    % ver. 1 (06/2013) - main code developed
    
% REFERENCES:
    % A. Koswara and R. Chakrabarti, Phys. Rev. A. (submitted)
    % A. Koswara and R. Chakrabarti, IEEE trans. evol. comp. (to be
        % submitted)

%% Linear combination of sine waves    
len=length(t);
efield=zeros(1,len);
nomodes=length(amp);
        
for m=1:nomodes
    efield=efield+amp(m)*sin(omega(m)*t+phi(m));
end

                
                