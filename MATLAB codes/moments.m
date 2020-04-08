%% Calculation of expected transition amplitude

function [Ujigamexp,Ujigam,Cjioff]=moments(fftUamp,alpha,ampoffexp,...
    ampoff,type)

%DESCRIPTION:
    % moments calculate the first or second moment of a quantum transition 
        % amplitude Ujigam given the expectation values of amplitude modes
        % up to certain orders. 

%INPUT:
    %fftUamp: MI amplitude from i to j of order 1 to mmax
        %its index is associated with a unique sum of encoding frequency 
        %gamma i.e. Ujigam(1) corresponds to frequency gamma(1)
    %alpha: encoding frequencies sorted based on orders of encoded pathway
    %ampoffexp: E[A^\alpha_i] for 1<=i<=mmax     
    %ampoff: contributing amplitude modes with offset
    %type: "first" or "second" moment

%OUTPUT:   
    % Ujigamexp: first or second moment of Ujigam (Cjioff in the paper)
    % Ujigam: (for checking with time propagated values)
    % Cjioff: normalized Ujigam(T,\gamma) wrt to amplitude pathways
    % sumterm1prod1: 1st sum term in the equation for E[Ujigam] (for checking
        % purposes)
    % sumterm2prod2: 2nd sum term in the equation for E[Ujigam] (same as
        % above)
    
%WRITTEN BY:
    % Andy Koswara
    % Advisor: Prof. Raj Chakrabarti
    % School of Chemical Engineering
    % Purdue University

%VERSION HISTORY:
    % ver. 1 - written to calculate E[U_{ji}] 
    % ver. 2 - var[Re,Im c_{ji}] included        
    % ver. 3 - calculation of first and second moments are separated.
        % calculation of E[A_k^\alpha_k] are now done outside of moments
        % code 
        
%Note: Need mechanism to show that none of the sum or product terms are 
    % below eps

%% input parameters
nomodes=size(ampoff,2);    
dim=size(fftUamp,1);
N=size(alpha,1);
mmaxamp=alpha(end,nomodes+1)+1;
Ujigamexp=zeros(dim,dim,mmaxamp);
Ujigam=zeros(dim,dim,mmaxamp);       

% U_ji(T,m)
for n=1:N
    Ujigam(:,:,alpha(n,nomodes+1)+1)=Ujigam(:,:,alpha(n,nomodes+1)+1)+...
        fftUamp(:,:,n);
end

% C_ji(T,\gamma)
	%amplitude product terms  
amprod=ones(1,N);
for n=1:N
    for k=1:nomodes
        if alpha(n,k)~=0
            amprod(n)=amprod(n)*ampoff(k)^alpha(n,k); 
                % A_1^\alpha_1*...*A_k^\alpha_k
        else
            continue;
        end
    end
end
        
Cjioff=zeros(dim,dim,N);
for n=1:N
    Cjioff(:,:,n)=fftUamp(:,:,n)/amprod(n);
end

switch type

%% first moment
    case 'first'% 1st moment
        
        sumterm1prod1=ones(1,N);
        for n=1:N
            for k=1:nomodes
                if alpha(n,k)~=0
                    sumterm1prod1(n)=sumterm1prod1(n)*...
                    	ampoffexp(k,alpha(n,k)); % E[A_1^\alpha_1]*...*E[A_k^\alpha_k]
                else
                    continue;
                end
            end
        end
        
        for n=1:N
            Ujigamexp(:,:,alpha(n,nomodes+1)+1)=...
            	Ujigamexp(:,:,alpha(n,nomodes+1)+1)+...
            		Cjioff(:,:,n)*sumterm1prod1(n); % E[U_ji(T,m)]
        end
        
%% second moment
    case 'second'% 2nd moment (calculated in separate file for now)
   
end