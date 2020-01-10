%-SYNTHETIC EXPERIMENT-%
% -----------------------------------------------------------------------
% Author: Sujit Kumar Sahoo, School of Electrical and Electronic
% Engineering, Nanyang Technological University, Singapore.
% -----------------------------------------------------------------------
% Results for the paper "Dictionary Training for Sparse Representation as
% Generalization of K-means Clustering", written by S.K. Sahoo and A. Makur
% and to be appeared in the IEEE Signal Processing Letter, 2013.
% -----------------------------------------------------------------------
% Courtesy : M. Aharon and M. Elad, Department of Computer Science,
% Technion, Haifa, Israel
% -----------------------------------------------------------------------
% in this file a synthetic test of the dictionary training algorithms are
% performed. First, a random dictionary with normalized columns is being
% generated, and then a set of data signals, each as a linear combination
% of 'm' dictionary element is created, with SNR of 10, 20, 30 and inf dB.
% this set is given as input to the dictionary training algorithms.

% a different mode for activating the dictionary training algorithms is to
% run Sparse coding stage until a fixed error is reached in the , instead
% until a fixed number of coefficients is found (it was used by us for the
% denoising experiments). in order to switch between those two modes just
% change the param.errorFlag (0 - for fixed number of coefficients, 1 -
% until a certain error is reached).
%-------------------------------------------------------------------------
NoTest = 1; %- N. of Test Cases

SNRdB = 20; %-Signal to Noise Raion (SNR) in dB

%---creating the data to train on ----%
N = 1500; % number of signals to generate
n = 20;   % dimension of each data
m = 3;      % no of sparse coefficients

%---initial dictionary: Dictionary elements---%
param.K = 50; % number of dictionary elements
Dictionary = randn(n,param.K);
Dictionary = Dictionary*diag(1./sqrt(sum(Dictionary.*Dictionary)));
param.TrueDictionary = Dictionary;
param.errorFlag = 0; % decompose signals to a fix number of coefficients,
% do not continue until a certain error is reached.
param.errorGoal = 0;
param.preserveDCAtom = 0;
param.InitializationMethod =  'DataElements';
param.displayProgress = 1;

param.m = m;   % number of elements in each linear combination.
param.numIteration = floor(9*m^2); % number of iteration for training.

% Initializing Computation Time
TimeKSVD=0;
TimeSGK=0;

% Initializing Atoms identitfied in each Iteration for each test
% cases for diffrent SNR
RecoveryKSVD = zeros(NoTest,param.numIteration+1);
RecoverySGK = zeros(NoTest,param.numIteration+1);


for test = 1:NoTest
    % Training Data Generated with aditive noise
    [~, ~, data] = gererateNoiseAddedSyntheticData(N, ...
        param.m, param.TrueDictionary, SNRdB);
    
    % Dictionary training algorithms starts
    disp(' ')
    disp('Starting K-SVD');
%     [~,output]  = KSVDnoiseS(data,param);
    tic;
    [~,output]  = KSVD(data,param);
    toc;
    TimeKSVD = TimeKSVD + output.time;
    disp(['K-SVD retrived ',num2str(output.ratio(end)),...
        ' atoms from the original dictionary']);
    RecoveryKSVD(test,:) = output.ratio;
  
    disp(' ')
    disp('Starting SGK');
%     [~,output]  = SGKnoiseS(data,param);
    tic;
    [~,output]  = SGK(data,param);
    toc;
    TimeSGK = TimeSGK + output.time;
    disp(['SGK retrived ',num2str(output.ratio(end)),...
        ' atoms from the original dictionary']);
    RecoverySGK(test,:) = output.ratio;
    
end

%% Displaying the Tabulated Result
% Calculating average Computation Time for each Training algorithm
TimeKSVD = TimeKSVD/NoTest;
TimeSGK = TimeSGK/NoTest;

% Displaying the obtained Results
disp(['Per iteratioon average Compuation Time in Milli Second, for m = ', num2str(m)])
disp(['KSVD = ',num2str(TimeKSVD),' SGK = ',num2str(TimeSGK)])
disp(' ')
disp(['AVERAGE NO. OF ATOMS RETRIEVED BY DICTIONARY TRAINING for m = ',...
    num2str(m)])

disp(['For SNR  = ', num2str(SNRdB)])
avgKSVD = mean(RecoveryKSVD(:,end))*param.K/100;
avgSGK = mean(RecoverySGK(:,end))*param.K/100;
disp(['KSVD = ',num2str(avgKSVD),' SGK = ',num2str(avgSGK)])
disp(' ')


%% Convergence Plot
figure,hold; 
disp(['Ploting average % of atoms retrieved in each iteration, for m = ',...
    num2str(m)])
disp(' ')

x = zeros(1,size(RecoveryKSVD,2));

x(:)=mean(RecoveryKSVD);
plot(x,'-g')
x(:)=mean(RecoverySGK);
plot(x,'-.b')

hold;
title(['SNR = ',num2str(SNRdB),'dB']);
xlabel('Iteration No.'),ylabel('Average % of atoms recovered');
legend('K-SVD','SGK');