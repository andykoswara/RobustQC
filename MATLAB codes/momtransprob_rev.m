%% Calculation of first moment of quantum transition probability 

function [pjiexp,pathwaynormexp,interfexpm]=momtransprob_rev(Cji,mmax,...
    alphaord,alphaind,ampoffexp,inistate,finstate)

% DESCRIPTION:
    % calculates the expected transition probability given quantum 
        % pathways and first moment of parameters.

% INPUT:
	% Cji: normalized amplitude pathways.
	% mmax: assumed maximum pathway order. 
    % alphaord: list of amplitude pathways of all order up to mmax.
    % alphaind: pathway index associated with Cji.
    % ampoffexp: expectation values of amplitude modes.
    % inistate: initial state.
    % finstate: final state.
    
% OUTPUT:
	% pjiexp: expected transition probability E[P_{ji}^m].
    % pathwaynorm: pathway norm or the first term in expected Pji
        % expression.
    % interfexpm: pathway interference of the second term in the 
        % expression.    
    
% WRITTEN BY:
	% Andy Koswara,
	% Advisor: Prof. Raj Chakrabarti,
	% School of Chemical Engineering,
	% Purdue University.

% VERSION:
    % ver. 1 - code written   
    % ver. 2 - terms involving pathway norms and interference are 
        % separated
    
%REFERENCE
    % A. Koswara and R. Chakrabarti, Robustness of controlled quantum 
        % dynamics, Phys. Rev. A (submitted).

%% Calculation of E[P_ji^m]
pathwaynormexp=zeros(1,mmax); 
interfexpm=zeros(1,mmax);
K=size(ampoffexp,1); % number of modes

nopath=size(alphaord,1);
for m=1:mmax% iterate over m
    for i=1:nopath
        if [alphaord(i,:,m),i==1]==[zeros(1,K+1),0]
            break; % no more meaningful pathways
        else
            for mp=1:m % iterate over m prime
                if m==mp
                    lmax=i;% to prevent from double-counting
                else
                    lmax=nopath; % otherwise account for all combination
                    % of pathways
                end
                for l=1:lmax
                    if [alphaord(l,:,mp),l==1]==[zeros(1,K+1),0]
                        break; % same as above
                    else
                        amptemp=1; %product of expected amplitude is
                        % initialized at every pathway combination
                        for k=1:K
                            if alphaord(l,k,mp)==0 && alphaord(i,k,m)==0
                                continue; % E[A]
                            else
                                amptemp=amptemp*ampoffexp(k,alphaord(i,k,m)+...
                                    alphaord(l,k,mp)); % prod_k E[A_k^\alpha_k + \alpha_k']
                            end
                        end
                        cjirealtemp=...
                            real(Cji(finstate,inistate,alphaind(i,m)))*...
                            real(Cji(finstate,inistate,alphaind(l,mp)));
                        cjiimagtemp=...
                            imag(Cji(finstate,inistate,alphaind(i,m)))*...
                            imag(Cji(finstate,inistate,alphaind(l,mp)));
                        if mp==m && l==i % calculation of the 1st
                            % term
                            pathwaynormexp(m)=pathwaynormexp(m)+amptemp*...
                                (cjirealtemp+cjiimagtemp);
                        else % 2nd term
                            interfexpm(m)=interfexpm(m)+2*amptemp*...
                                (cjirealtemp+cjiimagtemp);
                        end
                    end
                end
            end
        end
    end
end

pjiexp=pathwaynormexp+interfexpm;
