%% Quantum control mechanism identification (MI) - amplitude pathways

function [fftUamp,gammamod]=MIamp(H0,mu,amp,omega,phi,t,mmax,dispbool)

% DESCRIPTION:
    % Determine combinations of contributing field's spectral amplitude in 
        % the Dyson's Series (amplitude pathways) using Fourier encoding 
        % method followed by FFT.

% INPUT:
    % H_0: System's time-independent Hamiltonian,
    % mu: System's dipole moment,
    % amp: Control field's amplitude modes,
    % omega: Control field's frequency modes,
    % phi: Control field's phase modes,
    % t: Control field's time vector,
    % mmax: (assumed) maximum perturbation order,
    % dispbool: display option.

% OUTPUT:
    % fftUamp: decoded amplitude pathways sorted increasing order.
    % gammamod: amplitude modulation frequency

% WRITTEN BY: 
    %Andy Koswara,
    %Advisor: Prof. Raj Chakrabarti,
    %School of Chemical Engineering,
    %Purdue University.     
    
% VERSION HISTORY
    % ver. 1 (05/2013) - main code developed.
    % ver. 2 (05/2013) - spectral field's profile and maximum order are 
        % defined as inputs.
    % ver. 3 (06/2013) - error calculation is moved onto a different 
        % function.   
    % ver. 4 (04/09/2013) - definition of mmax is modified (shifted with
        % +1) for consistency purposes only

%REFERENCE
    % A. Koswara and R. Chakrabarti, Robustness of controlled quantum 
        % dynamics, Phys. Rev. A (2014) (submitted).

%% MI parameter settings
tic
nomodes=size(amp,2);
n=ceil(log((mmax+1)^nomodes)/log(2)); % for e.g. assuming a maximum pathway 
    % order of 9 associated with each of the 2 possible modes, 10^2=100 
    % covers \gamma of all amplitude pathways without overlapping 
    % (i.e. \alpha_1*\gamma_1 is in the set [1,2,...,9] and \alpha_2*\gamma_2
    % in [10,20,...,90]. Note that there is 10 extra n's not actually needed. So
    % really should have been 10^2-10.
N=2^n; % encoding frequency range
gammacons=2*pi/N; 
len=size(t,2);
dt=t(2)-t(1);
dim=size(H0,1);
im=sqrt(-1);
Uamps=zeros(dim,dim,N);

%% Amplitude modulating frequency
gammamod=zeros(1,nomodes);
ct=1;
for m=1:nomodes
    gammamod(m)=ct*gammacons;
    ct=ct*(mmax+1);% $\gamma_k$ needs to be "far" enough from each other
end

%% Amplitude encoding
for s=1:N% iterate over dummy encoding variable
    if dispbool
        disp(['working on s=' num2str(s)]);
    end
    
    Uamps(:,:,s)=eye(dim);
    for i=1:len % time propagation
        
        % amplitude encoding
        viamps=zeros(dim);
        for m=1:nomodes 
            viamps=viamps+expm(im*H0*t(i))*mu*amp(m)*...
                exp(im*gammamod(m)*(s-1))*...
                sin(omega(m)*t(i)+phi(m))*expm(-im*H0*t(i));
        end
        
        Uamps(:,:,s)=expm(-im*viamps*dt)*Uamps(:,:,s); % U(T,s)
    end
end

%% Pathway extraction via FFT
fftUamp=zeros(dim,dim,N);
for i=1:dim
   for j=1:dim
       fftUamp(i,j,:)=fft(Uamps(i,j,:))/N; 
   end
end
disp(['time consumed in amplitude MI: ' num2str(toc)]);