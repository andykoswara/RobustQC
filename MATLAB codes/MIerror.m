%% MI error calculation

function [errors_fftUji]=MIerror(Uint,fftUji,finstate,inistate)    

% DESCRIPTION:
    % calculates the error due to MI analysis for a specific final and initial
        % state

% INPUT:
    %Uint: unitary propagator in the interaction picture,
    %fftUji: decoded unitary propagator,
    %finstate: final state (of interest),
    %inistate: initial state (same).

% OUTPUT:
    %errors_fftUji: difference between U_{ji}(T) calculated by time-order
        % expansion and post MI analysis.

% WRITTEN BY:
    %Andy Koswara
    %Advisor: Prof. Raj Chakrabarti
    %School of Chemical Engineering
    %Purdue University     
    
% VERSION HISTORY
    % ver. 1 - main code developed
    % ver. 2 (04/09/14) - finstate and inistate variable added
    
% NOTE
    % May be made more sophisticated by including accurate determination of
        % Nmax

%% Error calculation
N=size(fftUji,3);
errors_fftUji=nan(1,N);
for i=1:N
    errors_fftUji(i)=norm(sum(fftUji(finstate,inistate,1:i),3)-...
        Uint(finstate,inistate));
end
