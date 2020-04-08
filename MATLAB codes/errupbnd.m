%% Calculation of upper bounds of C_{ji}^m and U_{ji}^m 

function [cjiupbnd, ujiexpupbnd, ampprodsum]=errupbnd(H0,mu,T,mmax,...
    alphaord,ampexp)

% DESCRIPTION:
    % Calculation of c_ji's and E[U_ji^m]'s upper bounds

% INPUT:
	% H0: system’s Hamiltonian
	% mu: system’s dipole moment
	% mmax: assumed maximum order of perturbation 
	% T: field duration
    % alphaord: amplitude pathways for all order
    % ampexp: expectation values of amplitude modes
    
% OUTPUT:
	% cjiupbnd: upper bound error on normalized transition amplitude
	% ujiexpupbnd: upper bound error on E[U_{ji}^m]
    % ampprodsum: sum of product amplitude terms
    
% WRITTEN BY:
	% Andy Koswara
	% Advisor: Prof. Raj Chakrabarti
	% School of Chemical Engineering
	% Purdue University

% VERSION:
    % ver. 1 - code written   
    
% REFERENCES:
	% ongoing work

%% Run    
d=max(max(mu)); % maximum element of mu
N=size(H0,1); % size of H0;
cjiupbnd=zeros(1,mmax); % upper bound(c_ji) 
K=size(ampexp,1); % number of modes
nopath=size(alphaord,1); % max no of pathways
ampprodsum=zeros(1,mmax); % sum of product of expected amplitudes
ujiexpupbnd=zeros(1,mmax);% upper bound (E[U_ji^m])

% Calculation of product of expected amplitudes (ampprodsum)
for m=1:mmax
    ord=m-1;
%     ord=m; % wrong
    cjiupbnd(m)=(N*d*T)^ord/(N*factorial(ord));
    ampprod=ones(1,nopath); % product of expected amplitude
    for i=1:nopath % the first is the 0-th pathway
        if m==1 % 0-th order pathway
            ampprod(2:end)=0;
            break; % no more meaningful pathways
        elseif alphaord(i,:,m)==zeros(1,K+1)
            ampprod(i:end)=0; % no more meaningful pathways
            break;
        else
            for k=1:K
                if alphaord(i,k,m)==0
                    continue; % expectation value of A is 1
                else
                    ampprod(i)=ampprod(i)*ampexp(k,alphaord(i,k,m));
                end
            end
        end
    end
    ampprodsum(m)=sum(ampprod);
    ujiexpupbnd(m)=cjiupbnd(m)*ampprodsum(m);
end