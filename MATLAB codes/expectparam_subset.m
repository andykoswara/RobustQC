%% Expectation of "a subset" of parameter of different orders 
    
function [paramexp,paramoff,offset]=expectparam_subset(param,subset,...
    sigma,mmaxamp,offs)

% DESCRIPTIONS:
    % calculates the expectation value of "a subset" of field's or system's 
        % parameters up to certain powers given a Gaussian noise distribution 
        % sigma.

% INPUT:
    %param: a vector of parameters (\theta) contributing to the transition 
        % U_ji.
    %subset: 0 and 1 vector indicating which parameters are to be 
        % calculated
    %sigma: a vector of noisy distribution parameter.
    %mmax: pre-determined significant order in the Dyson Series.
    %offs: with/without offset option.   

% OUTPUT:   
    % ampexp: a matrix of E[A_{k}^\alpha_i] for 1<=i<=mmax, and 1<=k<=no of modes
    % paramoff: parameter with offset (offset_{k}*A_{k}). 
    % offset: a vector of offsets for which E[offset_{k}*A_k^\alpha_i]
        % does not go above or below eps.
    
% WRITTEN BY:
    % Andy Koswara
    % Advisor: Prof. Raj Chakrabarti
    % School of Chemical Engineering
    % Purdue University
    
% VERSION HISTORY:
    % ver. 1. - main code developed
    % ver. 2. - input amp, originally a scalar, is modified into a vector. Output
        % ExpectA is now correspondingly a matrix of
        % [E[A_1^1],...,E[A_1^mmax]];...;[E[A_k^1],...,E[A_k^mmax]]. 
        % In addition, the scalar alpha is replaced with mmax.
    % ver. 3. - offset scalar added to prevent E[A]'s from going above or below
        % eps
	% ver. 4. - modified to include a subset of field modes
    % ver. 5 (09/2014) - code modified for generalized expected parameter
        % calculations including for dipole parameters.

%REFERENCE
    % A. Koswara and R. Chakrabarti, (to be published) (2014).

%% Initizalization            
noparams=size(param,2);
nosubparams=nnz(subset);
paramtemp=nan(1,nosubparams);
sigmatemp=nan(1,nosubparams);
ind=1;
for i=1:noparams
	if subset(i)
		paramtemp(ind)=param(i);
		sigmatemp(ind)=sigma(i);
		ind=ind+1;
	end
end

%taking into account only subset of parameters
param=paramtemp;
sigma=sigmatemp;
noparams=nosubparams;
paramexp=nan(noparams,mmaxamp);% E[A^\alpha_i] for all i's are saved 
    % for efficient computation
noise=sigma.*param;    

%% Calculating offset
if offs
%     calculating offsets
    errpslb=1E-1;
    errpsub=1E10;
    ampmin=min(param);
    ampmax=max(param);
    if ampmin < 1 && ampmax < 1
        check=ampmin^mmaxamp;
        i=1;
        offset=ones(noparams,mmaxamp);
        while check < errpslb
            i=i+1;
            check=(i*ampmin)^mmaxamp;
            if check > errpslb
                check2=(i*ampmax)^mmaxamp;
                if check2 < errpsub
                    offset=i;
                    paramoff=offset*param;
                else
                    error('offset value may need to be different for each A_k^\alpha_k');
                end
            else
                continue;
            end
        end
    else
        warning('the case for |A_k| >= 1 is not yet accounted for');
        offset=1;
        paramoff=param;
    end
else
    offset=1;
    paramoff=param;
end

%% Calculating E[offset*A_k^\alpha_k]
sumterm=zeros(noparams,mmaxamp);
for k=1:noparams
    paramexp(k,1)=paramoff(k);
    paramexp(k,2)=(paramoff(k))^2+ (offset*noise(k))^2;
    if mmaxamp==1 || mmaxamp==2
        return;
    else
        for m=3:mmaxamp % calculating E[offset*A_k^\alpha_k] recursively 
                % starting from alpha_k=3
            for i=1:m % calculating sum term
                ind=i-1;
                if ind==0
%                     sumindterm(k,m,i)=(-1)^m;
                    sumterm(k,m)=sumterm(k,m)+(-1)^m;
                elseif nchoosek(m,ind) > 1E15
                    error('binomial coefficient is above eps');                 
                elseif paramoff(k)^m < 1E-15
                    warning('A^\alpha is below eps');
                else
                    sumterm(k,m)=sumterm(k,m)+nchoosek(m,ind)*...
                        (paramexp(k,ind)/paramoff(k)^(ind))*(-1)^(m-ind);
                    
                    if sumterm(k,m) > 1E15
                        error('sumterm is above eps');
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