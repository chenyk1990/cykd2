function [coefs,dataC,data] = gererateNoiseAddedSyntheticData(N, m, Dictionary, SNRdB)
%-Generates Synthetic Training Data for give SNR-%
%------------------------------------------------------------------
% Courtesy : M. Aharon and M. Elad, Department of Computer Science, 
% Technion, Haifa, Israel
%------------------------------------------------------------------
% randn('state',sum(100*clock));
% rand('state',sum(100*clock))

[dataC,coefs] = CreateDataFromDictionarySimple(Dictionary, N, m);

if (SNRdB==0||SNRdB>=80) 
    data=dataC;
    return;
else
    noise = randn(size(dataC));
    actualNoise = calcNoiseFromSNR(SNRdB,dataC, noise);
    data =  dataC + actualNoise;
end

function [D,xOrig] = CreateDataFromDictionarySimple(dictionary, numElements, numCoef)

xOrig = zeros(size(dictionary,2),numElements);

maxRangeOfCoef = 1;
coefs = randn(numCoef,numElements)*maxRangeOfCoef;
xOrig(1:numCoef,:) = coefs;

for i=1:size(xOrig,2)
    xOrig(:,i) = xOrig(randperm(size(xOrig,1)),i);
end
% for i=1:size(xOrig,2)
%     xOrig(2:end,i) = xOrig(1+randperm(size(xOrig,1)-1),i);
% end

D = dictionary*xOrig;

function  actualNoise = calcNoiseFromSNR(TargerSNR, signal, randomNoise)
signal_2 = sum(signal.^2);
ActualNoise_2 = signal_2./(10^(TargerSNR/10));
noise_2 = sum(randomNoise.^2);
ratio = ActualNoise_2./noise_2;
actualNoise = randomNoise.*repmat(sqrt(ratio),size(randomNoise,1),1);

% function SNR = calcSNR(origSignal, noisySignal)
% errorSignal = origSignal-noisySignal;
% signal_2 = sum(origSignal.^2);
% noise_2 = sum(errorSignal.^2);
% 
% SNRValues = 10*log10(signal_2./noise_2);
% SNR = mean(SNRValues);
