%% Unitary propagator

function [U]=unitpropagator(H0,mu,efield,t,Ui,casetype)

% DESCRIPTION: 
    % unitpropagator calculates the unitary matrix either in the 
    % Schrodinger's or Heisenberg's (interaction) picture via time order
    % expansion of the Dyson's Series

% INPUT:
    % H0: system's free hamiltonian 
    % mu: dipole operator
    % efield: control field
    % t: time vector
    % casetype: types of propagation
    
% OUTPUT:
    % U: unitary matrix
    
% WRITTEN BY: 
    % Andy Koswara
    % Advisor: Prof. Raj Chakrabarti
    % School of Chemical Engineering
    % Purdue University
    
% VERSION HISTORY:
    % ver. 1 - main code developed for the heisenberg's case
    % ver. 2 - schrodinger's case added
        
%% Initialization and Run
U=Ui;    
dt=t(2)-t(1);
im=sqrt(-1);

switch casetype
    
    case 'schrodinger'
        
        for i=1:size(efield,2)
            H = H0 - mu*efield(1,i);
            U = expm(-im*H*dt)*U;
        end
        
    case 'heisenberg'
        
        for i=1:size(efield,2)
            vi = expm(im*H0*t(i))*mu*efield(1,i)*expm(-im*H0*t(i));
            U = expm(-im*vi*dt)*U;
        end
end