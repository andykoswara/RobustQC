%% Quantum pathway (nominal) interference calculations

function [pathwaynorm,interf,interfangle,interfm]=interfere(fftU,...
    alphaord,alphaind,mmax,finstate,inistate)

% DESCRIPTION:
    % Calculates the interference between different order quantum pathways
    
% INPUT:
    % fftU: unitary matrix broken down in terms of orders
    % alphaord: list of amplitude pathways of all order up to mmax
    % alphaind: pathway index associated with Cji
    % mmax: assumed maximum pathway order. 
    % inistate: initial state.
    % finstate: final state.

% OUTPUT:
    % pathwaynorm: norm of pathways (of the same order m)
    % interf: measure of interference
    % interfangle: angle of interference
    % interfm: measure of interferences broken down into different order m 

% WRITTEN BY:
    % Andy Koswara
    % Advisor: Prof. Raj Chakrabarti
    % School of Chemical Engineering
    % Purdue University
    
% VERSION HISTORY:
    % ver. 1 - main code developed
    % ver. 2 - modified to mimic the moment of interference calculation 
        % counterpart
    % ver. 3 (05/2014) - a measure of angle is included and interferences
        % between pathway of the same order (m) is included.

%REFERENCE
    % A. Koswara and R. Chakrabarti, Robustness of controlled quantum 
        % dynamics, Phys. Rev. A (2014).
    
nopath=max(max(alphaind));
K=size(alphaord,2)-1;
interf=zeros(mmax,mmax);
interfangle=zeros(nopath,nopath);
pathwaynorm=zeros(1,mmax);
interfm=zeros(1,mmax);
for m=2:mmax
    for i=1:nopath
        if alphaord(i,:,m)==zeros(1,K+1) % here the 0th order pathway is
            % ignored - should have no values unless for subset
            % encoding, and even then it should not be accounted for per
            % definition of quantum pathways)
            break; % no more meaningful pathways
        else
            for mp=1:m
                if mp==m
                    tempind=i;
                else
                    tempind=nopath;
                end
                for j=1:tempind
                    if alphaord(j,:,mp)==zeros(1,K+1)
                        break; % same as above
                    elseif mp==m && i==j
                        pathwaynorm(m)=pathwaynorm(m)+...
                            real(fftU(finstate,inistate,alphaind(i,m))*...
                            conj(fftU(finstate,inistate,alphaind(j,mp))));
                    else
                        interf(m,mp)=interf(m,mp)+2*...
                            real(fftU(finstate,inistate,alphaind(i,m))*...
                            conj(fftU(finstate,inistate,alphaind(j,mp))));
                        interfangle(alphaind(i,m),alphaind(j,mp))=acosd(...
                            real(fftU(finstate,inistate,alphaind(i,m))*...
                            conj(fftU(finstate,inistate,alphaind(j,mp))))/...
                            (norm(fftU(finstate,inistate,alphaind(i,m)))*...
                            norm(fftU(finstate,inistate,alphaind(j,mp)))));
                        interfm(m)=interfm(m)+2*...
                            real(fftU(finstate,inistate,alphaind(i,m))*...
                            conj(fftU(finstate,inistate,alphaind(j,mp))));
                    end
                end
            end
        end
    end
end