%% Sorting quantum pathway in terms of orders

function [alphaord,alphaind]=sortalpha(alpha)

% DESCRIPTION:
    % The funciton takes a 2D-matrix of alpha containing a list of all 
        % possible amplitude pathways of different orders and convert it 
        % into a 3-D matrix where each order pathway correspond to a 
        % 2D-block

% INPUT:
	% alpha: list of pathways
    
% OUTPUT:
	% alphaord: ordered list of pathways
    % alphaind: corresponding index of alphaord
    
% WRITTEN BY:
	% Andy Koswara
	% Advisor: Prof. Raj Chakrabarti
	% School of Chemical Engineering
	% Purdue University

% VERSION:
    % ver. 1 - code written   
    
% REFERENCES:
	% A. Koswara and R. Chakrabarti, Robustness of controlled quantum 
        % dynamics, Phys. Rev. A (submitted).

%% Group alpha's/amplitude pathways and cji according to their orders
N=size(alpha,1); %size of amplitude pathways of all orders
noparams=size(alpha,2)-1; % the last row stores the info on pathway
    % order
[alphasort,indexsort]=sortrows(alpha,noparams+1); % sort amplitude pathways 
    % based on its order
ordind=0; %index for keeping tab of pathway order
chunkind=1; %starting index for indicating the set of pathways in the 
    % alpha matrix belonging to the same order
mmax=max(alpha(:,noparams+1)); % maximum pathway order
alphaind=zeros(mmax+1); % each order pathway index  
alphaord=zeros(N,noparams+1,mmax+1); % matrix containing pathways sorted 
    % according to their order in each column
for i=1:N
    if alphasort(i,noparams+1)==ordind
        continue; % continue assigning pathways of the same order to the
            % same column
    else
        chunksize=size(indexsort(chunkind:i-1)); % size of pathways of
            % the same order.
        alphaind(1:chunksize,ordind+1)=indexsort(chunkind:i-1);% the 
            % corresponding indices of pathways belonging to the same 
            % order
        alphaord(1:chunksize,:,ordind+1)=alphasort(chunkind:i-1,:); 
            % assignmnet of pathways of the same order to the same column
        ordind=ordind+1;% moving onto the next order
        chunkind=i;%starting index of the next order
    end
end

% assigning last point (since the pathway assignment is always lagging by
% one)
alphaind(1,end)=indexsort(end);
alphaord(1,:,end)=alphasort(end,:);