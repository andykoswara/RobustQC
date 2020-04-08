%% Calculation of upper bounds of P_{ji}^m 

function [pjiexpupbnd,ampprodsum]=pjierrupbnd_rev_rc(H0,mu,T,mmax,...
    alphaord,Aexp)

% DESCRIPTION:
    % Calculation of the upper bound of E[P_ji^m]

% INPUT:
	% H0: system?s Hamiltonian
	% mu: system?s dipole moment
	% mmax: maximum order of perturbation 
	% T: field duration
    % alphaord: amplitude pathways for all order
    % Aexp: expectation values of amplitude modes
    
% OUTPUT:
	% pjiexpupbnd: upper bound on E[P_{ji}^m] for m \in [1,mmax]
    
% WRITTEN BY:
	% Andy Koswara
	% Advisor: Prof. Raj Chakrabarti
	% School of Chemical Engineering
	% Purdue University

% VERSION:
    % ver. 1 - code written   
    
% REFERENCES:
	% ongoing work

%% Calculation of P_ji's upper bounds
d=max(max(mu)); % maximum element of mu
N=size(H0,1); % size of H0;
K=size(Aexp,1); % number of modes

% calculating the first of many upperbound terms
pjiexpupbnd=zeros(1,mmax); % the first term in the calculation of upper 
    % bound
ampprodsum=zeros(mmax,mmax); % sum of product of expected amplitudes

% Calculating the product of E[A^\alpha's]
nopath=size(alphaord,1);
for m=1:mmax% iterate over m
    amptemp=1;
    for i=1:nopath
        if [alphaord(i,:,m)==zeros(1,K+1), i>1] % if m=1, 0-th order 
            % pathways count
            break; % no more meaningful pathways
        else
            for mp=1:m % iterate over m prime
                if m==mp
                    lmax=i;% to prevent from double-counting
                else
                    lmax=nopath;
                end 
                for l=1:lmax
                    if [alphaord(l,:,mp)==zeros(1,K+1), l>1]
                        break; % same as above
                    else
                        for k=1:K
                            if alphaord(l,k,mp)==0 && alphaord(i,k,m)==0
                                continue; % same as above
                            else
                                amptemp=amptemp*...
                                    Aexp(k,alphaord(i,k,m)+...
                                    alphaord(l,k,mp));
                            end
                        end
                    end
                    if mp==m && l==i % second term
                        ampprodsum(m,mp)=ampprodsum(m,mp)+2*amptemp;
                    else % first term
                        ampprodsum(m,mp)=ampprodsum(m,mp)+4*amptemp;
                    end
                    pjiexpupbnd(m)=pjiexpupbnd(m)+ampprodsum(m,mp)*...
                        (N*d*T)^(m+mp-2)/((N^2)* factorial(m-1)*...
                        factorial(mp-1));
%                      pjiexpupbnd(m)=pjiexpupbnd(m)+ampprodsum(m,mp)*...
%                         (N*d*T)^(m+mp)/((N^2)* factorial(m)*...
%                         factorial(mp));
                end
            end
        end
    end
end