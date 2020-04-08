%% Amplitude and dipole MI decoding

function [alpha]=decodeMI(mmax,noparams)

%DESCRIPTIONS:
    % Decoding pathways based on non-overlapping modulating frequency
        % scheme. For more info, please look into references.

%INPUT:
    % mmax: pre-determined maximum order in amplitude MI per mode. 
    % noparams: no of encoded parameters in the control field or system's 
        % dipole.
    
%OUTPUT:
    % alpha: matrix containing amplitude pathways sorted in ascending
        % order from 1 to mmax.
        
%WRITTEN BY:
    % Andy Koswara,
    % Advisor: Prof. Raj Chakrabarti,
    % School of Chemical Engineering,
    % Purdue University.
    
%VERSION HISTORY:
    % ver. 1 (2013) - main code written.
    % ver. 2 (2013) - bug fixes.
    % ver. 3 (04/09/2014) - generalized to include dipole decoding.
    
%REFERENCE
    % A. Koswara and R. Chakrabarti, Robustness of controlled quantum 
        % dynamics, Phys. Rev. A (2014).

%NOTE:
    % alpha looks like the set of all possible combinations or 
        % permutations from 1:mmax but not quite.

%% Initialize Decoding parameter 
% This is the $\gamma_k$ in the reference.  

shiftmmax=mmax+1;
N=2^(ceil(log(shiftmmax^noparams)/log(2)));
gammak=nan(1,noparams);
i=1;
for k=noparams:-1:1
    gammak(i)=shiftmmax^(k-1); %encoding frequencies in descending order
    i=i+1;
end

%% Pathway assignment
% Factorization of encoded frequency $\gamma$ based on $\gamma_k$ in order 
% to determine pathways ($\vec\alpha's$)

alpha=zeros(shiftmmax^noparams,noparams);
for n=1:N
    gammatemp=n-1; %starting at the 0th order
    alphatemp=zeros(1,noparams); % $\alpha_k$'s associated with a
        % particular $\gamma$
    if floor(gammatemp/gammak(1))<=mmax
        for k=1:noparams
            temp=mod(gammatemp,gammak(k));
            if temp==0
                alphatemp(k)=gammatemp/gammak(k);
                break; % all subsequent $\alpha_k$'s are zeros
            else
                alphatemp(k)=floor(gammatemp/gammak(k));
                gammatemp=gammatemp-alphatemp(k)*gammak(k); % moving onto
                    % the next set of $\alpha_k$
            end
        end
        
        alpha(n,:)=alphatemp;
    else
        break; %ignore higher order contribution
    end
end

%% Pathway sorting based on orders 
% sorting $\alpha_k$'s according to pathway orders $m$
alphaorder=sum(alpha,2);
alpha=horzcat(alpha,alphaorder);