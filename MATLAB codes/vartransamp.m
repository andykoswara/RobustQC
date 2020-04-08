%% Calculation of variance of transition amplitude

function [varujim]=vartransamp(Cji,mmax,alphaord,alphaind,ampoffexp,...
    inistate,finstate)

% DESCRIPTION:
    % Calculation of var(U_ji)

% INPUT:
	% Cji: normalized amplitude pathways
	% mmax: assumed maximum pathway order 
    % alphaord: list of amplitude pathways of all order up to mmax
    % alphaind: pathway index associated with Cji
    % ampoffexp: expectation values of amplitude modes
    
% OUTPUT:
	% varuji: var(U_{ji}^m)
    
% WRITTEN BY:
	% Andy Koswara
	% Advisor: Prof. Raj Chakrabarti
	% School of Chemical Engineering
	% Purdue University

% VERSION:
    % ver. 1 - code written   
    
% REFERENCES:
	% paper to be published

%% Run
% calculating the first of many upperbound terms
varujim=zeros(1,mmax); 
K=size(ampoffexp,1); % number of modes

% Calculating all E[P_ji^m] terms at once
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
                        %product of expected amplitude is
                            %initialized at every pathway combination
                        amptemp=1; 
                        amptemp2=1;
                        % may want to create another function to simplify 
                            % for this
                        for k=1:K
                            if alphaord(l,k,mp)+alphaord(i,k,m)==0
                                continue; % E[A]
                            else
                                amptemp=amptemp*ampoffexp(k,alphaord(i,k,m)+...
                                    alphaord(l,k,mp));% prod_k E[A_k^\alpha_k + \alpha_k']
                                if alphaord(l,k,mp)==0
                                    amptemp2=amptemp2*...
                                        ampoffexp(k,alphaord(i,k,m));
                                elseif alphaord(i,k,m)==0
                                    amptemp2=amptemp2*...
                                        ampoffexp(k,alphaord(l,k,mp));
                                else
                                    amptemp2=amptemp2*...
                                        ampoffexp(k,alphaord(i,k,m))*...
                                        ampoffexp(k,alphaord(l,k,mp));
                                end
                            end
                        end
                        realcjitemp=real(Cji(finstate,inistate,alphaind(i,m)))*...
                            real(Cji(finstate,inistate,alphaind(l,mp)));
                        imagcjitemp=imag(Cji(finstate,inistate,alphaind(i,m)))*...
                            imag(Cji(finstate,inistate,alphaind(l,mp)));
                        if mp==m && l==i % calculation of the 2nd
                            % term
                            varujim(m)=varujim(m)+(amptemp-amptemp2)*...
                                (realcjitemp+sqrt(-1)*imagcjitemp);
                        else % first term
                            varujim(m)=varujim(m)+2*(amptemp-amptemp2)*...
                                (realcjitemp+sqrt(-1)*imagcjitemp);
                        end
                    end
                end
            end
        end
    end
end