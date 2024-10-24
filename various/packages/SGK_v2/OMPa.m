function [S]=OMPa(Phi,V,m,alpha,errorGoal)
%==========================================================================
% Author: Sujit Kumar Sahoo, School of Electrical and Electronic
% Engineering, Nanyang Technological University, Singapore.
% -----------------------------------------------------------------------
% Sparse coding of a group of signals based on a given
% dictionary and specified number of atoms to use.
% input arguments: Phi - the sensing matrix( or dictionary).
%                  V - the signals to represent.
%                  m - the sparsity of the signals (No. of no zero 
%					coefficients)
% 				   alpha - the extended run factor, it should be 
%					0<alpha<1. Larger 'alpha' may increase the run 
%					time. 'alpha = 0.25' is sufficient to achieve 
%					the optimum performance. 
%                  errorGoal - the maximal allowed representation 
%					error for each signal, which should be greater 
%					than the noise variance. It can be set to  
%					10^-16 (~0) by default for noise-free case. 
% output arguments: S - sparse coefficient matrix.

% Reference: 
% [1] Sahoo, S.K.; Makur, A., "Signal recovery from random measurements via 
% extended orthogonal matching pursuit", Transactions on Signal Processing, 
% IEEE , vol.63, no.10, pp.2572-2581, May 2015.
%
% [2] Y. C. Pati, R. Rezaiifar, and P. S. Krishnaprasad,
% "Orthogonal matching pursuit: recursive function approximation with
% applications to wavelet decomposition", in Conference Record of the
% Asilomar Conference on Signals, Systems & Computers, vol. 1, 1993,
% pp. 40-44.
%===========================================================================
[N,d] = size(Phi);

if errorGoal == 0
    errorGoal = 1e-16;
end

maxNumCoef = floor((1+alpha).*m);
maxNumCoef(maxNumCoef>N)=N;

indx = zeros(1,maxNumCoef);
a = zeros(1,maxNumCoef);
Ai = zeros(maxNumCoef,maxNumCoef);

P = size(V,2);
S = zeros(d,P);

for j=1:1:P,
    x=V(:,j);
    if mean(abs(x).^2)< errorGoal
        continue;
    end
    %-- intializations--
    k=1;
    proj= Phi'*x;
    [~,pos]= max(abs(proj));
    pos=pos(1);
    indx(k)=pos;
    Ai(1,1) = 1;
    gama = sum(abs(Phi(:,indx(k))).^2);
    a(k) = proj(indx(k))/gama;
    %--residue calculation---
    residual=x-Phi(:,indx(1:k))*a(1:k).';
    currResNorm2 = mean(abs(residual).^2);
    
    while currResNorm2>errorGoal && k < maxNumCoef,
        k = k+1;
        %-- Atom Selection---
        proj= Phi'*residual;
        [~,pos]= max(abs(proj));
        pos=pos(1);
        indx(k)=pos;
        
        %--Intera atomic dependancy calculation--
        v = Phi(:,indx(1:k-1))'*Phi(:,indx(k));
        b = Ai(1:k-1, 1:k-1)*v;
        
        beta = 1/(1-v'*b);
        Ai(1:k,1:k) = [Ai(1:k-1, 1:k-1)+beta*(b*b'), -beta*b
            -beta*b', beta];
        %--coeficeient calculation--
        gama = sum(abs(Phi(:,indx(k))- Phi(:,indx(1:k-1))*b).^2);
        a(k) = proj(indx(k))/gama;
        a(1:k-1) = a(1:k-1) - a(k)*b.';
        
        %--residue calculation---
        residual=x-Phi(:,indx(1:k))*a(1:k).';
        currResNorm2 = mean(abs(residual).^2);
    end
    if k>0
        S(indx(1:k),j)=a(1:k);
    end
    
    a(1:k)=0;
    indx(1:k)=0;
end;
return;