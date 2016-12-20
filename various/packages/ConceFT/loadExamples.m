TN1 = round(N*.3);
TN2 = round(N*.4);
%% Generate the signal
% switch ExampleID
%     case 1
%         %% the amplitude modulation of the simulated signal
%         am1 = 2 + smooth(cumsum(randn(N,1)) ./ Hz, 200, 'loess') ;
%         am2 = 2 + smooth(cumsum(randn(N,1)) ./ Hz, 200, 'loess') ;
%         am1(1:TN1) = 0 ;
%         am2(end-TN2:end) = 0 ;
%         %% the instantaneous frequency of the simulated signal
%         if1 = 17 + 6 * smooth(cumsum(randn(N,1)) ./ Hz, 400, 'loess') ;
%         if2 = pi*2 + 2 * smooth(cumsum(randn(N,1)) ./ Hz, 300, 'loess') ;
%     case 2
%% the amplitude modulation of the simulated signal
am1 = smooth(cumsum(randn(N,1)) ./ Hz, 200, 'loess') ;
am1 = 2 + am1 ./ max(abs(am1)) ;
am2 = smooth(cumsum(randn(N,1)) ./ Hz, 200, 'loess') ;
am2 = 2 + am2 ./ max(abs(am2)) ;
am1(1:TN1) = 0 ;
am2(end-TN2:end) = 0 ;
%% the instantaneous frequency of the simulated signal
if1 = smooth(cumsum(randn(N,1)) ./ Hz, 400, 'loess') ;
if1 = tau + 6 * if1 ./ max(abs(if1)) ;
if2 = smooth(cumsum(randn(N,1)) ./ Hz, 300, 'loess') ;
if2 = pi*2 + 2 * if2 ./ max(abs(if2)) ;
% end
%%
phi1 = cumsum(if1) / Hz ;
phi2 = cumsum(if2) / Hz ;
if1(1:TN1) = nan ;
if2(end-TN2:end) = nan ;
%% the simulated signal.
clean = am1 .* cos(2*pi*phi1) + am2 .* cos(2*pi*phi2) ;

if ~exist('sigma','var')
    sigma = sqrt( var(clean)*10.^( -snrdb /10 ) );
end

%% Add noise
switch NoiseID
    
    case 1
        %% add noise (Gaussian white noise)
        disp('Gaussian')
        noise = randn(N, 1) ;
        xm = clean  + sigma * noise ;
        
    case 2
        %% ARMA
        disp('ARMA')
        dis = random('t',4,N, 1) ;
        e = armaxfilter_simulate(dis, .5, 1, .5, 1, -.5) ;
        noise = e ./ std(e) ;
        xm = clean + sigma * noise ;
        
    case 3
        %% poisson
        disp('Poisson')
        noise = random('poiss',1,N,1);
        xm = clean  + sigma * noise ;
end

actualSnr = 20 * log10(std(clean)./std(sigma * noise)) ;
fprintf(['actualSnr = ',num2str(actualSnr),'\n']) ;
fprintf(['actualSigma = ',num2str(sigma),'\n']) ;