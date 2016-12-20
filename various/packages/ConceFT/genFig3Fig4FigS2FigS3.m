clear ; close all ;
scrsz = get(0,'ScreenSize');
ftsz = 22;
Hz = 160;
T = 60;
time = [1/Hz:1/Hz:T]' ;
N = length(time) ;
ID = 2; % use 1 for Figure 3 and 4; 2 for EMS
seeds = [111 4]; % use the first seed for Figure 3 and 4; use the second for those of EMS
taus = [10 8];
%==================
snrdb = 0;
clear sigma
initstate(seeds(ID));
tau = taus(ID);

TN1 = round(N*.3);
TN2 = round(N*.4);
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

noise1 = randn(N, 1) ;
xm1 = clean  + sigma * noise1 ;

dis = random('t',4,N, 1) ;
e = armaxfilter_simulate(dis, .5, 1, .5, 1, -.5) ;
noise2 = e ./ std(e) ;
xm2 = clean  + sigma * noise2 ;

noise3 = random('poiss',1,N,1);
xm3 = clean  + sigma * noise3 ;


%=============================
	%% Figure 3 or S2
figure('Position',[1 scrsz(4)/2 scrsz(3) scrsz(4)/3])
plot(time, clean+18, 'color', 'k', 'linewidth', 1) ; hold on ;
plot(time, if1, 'color', [.7 .7 .7], 'linewidth', 3) ;
plot(time, if2, 'color', [.7 .7 .7], 'linewidth', 3) ;
axis([15 40 2 inf]) ; set(gca, 'fontsize', 22) ;
xlabel('Time (Sec)') ; ylabel('Frequency (Hz)') ;


%=============================
	%% figure 4 or S3
figure('Position',[1 scrsz(4)/2 scrsz(3) scrsz(4)/2])
a1 = 1.1*range(xm1)/2 ; a2 = 1.1*range(xm2)/2 ; a3 = 1.1*range(xm3)/2 ; a4 = 1.1*range(clean)/2 ;
a = a1+a2 ; b = a1+2*a2+a3 ; c = a1+2*a2+2*a3+a4 ;
plot(time, clean + c, 'color', 'k', 'linewidth', 1.2) ; hold on ; 
plot(time, xm1 + b, 'color', [.5 .5 .5], 'linewidth', 1.2) ; hold on ; 
plot(time, xm2 + a, 'color', [.5 .5 .5], 'linewidth', 1.2) ;
plot(time, xm3, 'color', [.5 .5 .5], 'linewidth', 1.2) ;
set(gca, 'ytick', [0 a b c]) ; set(gca, 'fontsize', 22) ;
set(gca, 'yticklabel', ['  Poisson'; 'ARMA(1,1)'; ' Gaussian'; '    Clean']) ;
axis([15 40 -inf inf]) ;
