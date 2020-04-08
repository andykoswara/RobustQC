%% Moment of quantum pathway interference

function [interfexpmmp,interfexpall]=mominterf(Cjioff,alphaord,alphaind,...
    paramoffexp,mmax,finstate,inistate)

% DESCRIPTION:
    % calculate expected interferences between quantum pathways

% INPUT:
    % Cjioff: normalized unitary matrix broken down by pathways
    % alphaord: list of amplitude pathways ordered according to their 
        % orders
    % alphaind: index of ordered pathways
    % paramoffexp: expected parameters (amplitude or mu) (with offset built
        % -in)
    % mmax: maximum pathway order taken into acccount
    % finstate: final state
    % inistate: initial state

% OUTPUT:
    % interfexpmmp: expected interference between pathways sorted in terms 
        % of orders m and m'
    % interfexpall: expected interference between individual pathways 
    
% WRITTEN BY:
    % Andy Koswara
    % Advisor: Prof. Raj Chakrabarti
    % School of Chemical Engineering
    % Purdue University

% VERSION:
    % ver. 1 - code written
    
% REFERENCES:
    % A. Koswara and R. Chakrabarti, Phys. Rev. A. (submitted)
    % A. Koswara and R. Chakrabarti, IEEE trans. evol. comp. (to be
        % submitted)

% interference calculation
interfexpmmp=zeros(mmax,mmax); %interference btw pathways
nopath=size(alphaord,1);
interfexpall=zeros(nopath,nopath);
K=size(paramoffexp,1);
for m=2:mmax% iterate over m
    for i=1:nopath
        if alphaord(i,:,m)==zeros(1,K+1) % here the 0th order pathway is
            % ignored - should have no values unless for subset
            % encoding (even then it is ignored given definition of quantum
            % pathways)
            break; % no more meaningful pathways
        else
            for mp=1:m % iterate over m prime
                if m==mp
                    tempind=i-1;
                else
                    tempind=nopath;
                end
                for l=1:tempind
                    if alphaord(l,:,mp)==zeros(1,K+1)
                        break; % same as above
                    else
                        amptemp=1; %product of expected amplitude is
                        % initialized at every pathway combination
                        for k=1:K
                            if alphaord(l,k,mp)==0 && alphaord(i,k,m)==0
                                continue; % E[A]
                            else
                                amptemp=amptemp*paramoffexp(k,...
                                    alphaord(i,k,m)+alphaord(l,k,mp)); 
                                % prod_k E[A_k^\alpha_k + \alpha_k']
                            end
                        end                    
                        cjirealtemp=...
                            real(Cjioff(finstate,inistate,alphaind(i,m)))*...
                            real(Cjioff(finstate,inistate,alphaind(l,mp)));
                        cjiimagtemp=...
                            imag(Cjioff(finstate,inistate,alphaind(i,m)))*...
                            imag(Cjioff(finstate,inistate,alphaind(l,mp)));
                        interfexpmmp(m,mp)=interfexpmmp(m,mp)+2*amptemp*...
                            (cjirealtemp+cjiimagtemp);
                        interfexpall(alphaind(i,m),alphaind(l,mp))=2*...
                            amptemp*(cjirealtemp+cjiimagtemp);
                    end
                end
            end
        end
    end
end