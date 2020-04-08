%% Expectation of parameter of different orders 

function [paramexp,paramoff,offset]=expectparam(param,sigma,mmax,offs)

% DESCRIPTION:
    % calculates the expectation value of a field's or system's parameters
        % up to certain powers given a Gaussian noise distribution sigma.

% INPUT:
    %param: a vector of parameters (\theta) contributing to the transition 
        % U_ji.
    %sigma: a vector of noisy distribution parameter.
    %mmax: pre-determined significant order in the Dyson Series.
    %offs: with/without offset option. 

% OUTPUT:   
    % paramexp: a matrix of E[(offset_{k}*\theta_{k})^\alpha_i] 
        % for 1<=i<=mmax, and 1<=k<=no of modes.
    % paramoff: parameter with offset (offset_{k}*A_{k}). 
    % offset: a vector of offsets for which E[offset_{k}*A_k^\alpha_i]
        % does not go above or below eps.
    
% WRITTEN BY:
    % Andy Koswara,
    % Advisor: Prof. Raj Chakrabarti,
    % School of Chemical Engineering,
    % Purdue University.
    
% VERSION HISTORY:
    % ver. 1 (2013) - main code developed
    % ver. 2 (2013) - input amp, originally a scalar, is modified into a vector. Output
        % ExpectA is now correspondingly a matrix of
        % [E[A_1^1],...,E[A_1^mmax]];...;[E[A_k^1],...,E[A_k^mmax]]. 
        % In addition, the scalar alpha is replaced with mmax.
    % ver. 3 (2013) - offset scalar added to prevent E[A]'s from going above or below
        % eps.
    % ver. 4 (04/2014) - code modified for generalized expected parameter
        % calculations including for dipole parameters.
        
%REFERENCE
    % A. Koswara and R. Chakrabarti, Robustness of controlled quantum 
        % dynamics, Phys. Rev. A (2014).
        
%NOTE:
    % Code may be improved by calculating offsets for different parameter
        % values.

%% Initizalization        
noparams=size(param,2);
paramexp=nan(noparams,mmax);% E[A^\alpha_i] for all i's are saved for efficient
    % computation
noise=sigma.*param;    

%% Calculating offset
if offs
%     calculating offsets
    errpslb=1E-1;
    errpsub=1E10;
    parammin=min(param);
    parammax=max(param);
    if parammin < 1 && parammax < 1
        check=parammin^mmax;
        i=1;
        offset=ones(noparams,mmax);
        while check < errpslb
            i=i+1;
            check=(i*parammin)^mmax;
            if check > errpslb
                check2=(i*parammax)^mmax;
                if check2 < errpsub
                    offset=i;
                    paramoff=offset*param;
                else
                    error('offset value may need to be different for each \theta_k^\alpha_k');
                end
            else
                continue;
            end
        end
    else
        warning('the case for |\theta_k| >= 1 is not yet accounted for');
        offset=1;
        paramoff=param;
    end
else
    offset=1;
    paramoff=param;
end

%% Calculating E[offset*A_k^\alpha_k]
sumterm=zeros(noparams,mmax);
for k=1:noparams
    paramexp(k,1)=paramoff(k);
    paramexp(k,2)=(paramoff(k))^2+ (offset*noise(k))^2;
    if mmax==1 || mmax==2
        return;
    else
        for m=3:mmax % calculating E[offset*A_k^\alpha_k] recursively 
                % starting from alpha_k=3
            for i=1:m % calculating sum term
                ind=i-1;
                if ind==0
                    sumterm(k,m)=sumterm(k,m)+(-1)^m;
                elseif nchoosek(m,ind) > 1E15
                    error('binomial coefficient is above eps');                 
                elseif paramoff(k)^m < 1E-15
                    warning('A^\alpha is below eps');
                else
                    sumterm(k,m)=sumterm(k,m)+nchoosek(m,ind)*...
                        (paramexp(k,ind)/paramoff(k)^(ind))*(-1)^(m-ind);
                    
                    if sumterm(k,m) > 1E15
                        warning('sumterm is above eps');
                    end
                end
            end
            sumterm(k,m)=sumterm(k,m)*paramoff(k)^m;
            if mod(m,2)==0 % alpha is even
                paramexp(k,m)=(m-1)*(offset*noise(k))^m-sumterm(k,m); 
                    %the sum term lags by 1 iteration
            else % alpha is odd
                paramexp(k,m)=-sumterm(k,m);
            end
        end
    end
end