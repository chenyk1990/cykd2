clear ; close all ;

addpath('/Users/hautiengwu/Dropbox/ConCeft/SubmissionCode/Conceft') ;
scrsz = get(groot,'ScreenSize');

	%% generate the simulated data
	%% she sampling time (100Hz sampling rate)
	%% high sampling rate to avoid sampling issue
Hz = 100 ;
time = [1/Hz:1/Hz:16]' ;
N = length(time) ;

        %% Setup parameters
        %% the number of random Multitaper (rMT)
MT = 100 ;
DDD = 1200*1000*5 ;

initstate(1) ;
dim = 2 ;


	%% the amplitude modulation of the simulated signal
am1 = smooth(cumsum(randn(N,1)) ./ Hz, 200, 'loess') ;
am1 = 2 + am1 ./ max(abs(am1)) ;
am2 = smooth(cumsum(randn(N,1)) ./ Hz, 200, 'loess') ;
am2 = 2 + am2 ./ max(abs(am2)) ;
am1(1:500) = 0 ;
am2(end-600:end) = 0 ;


    %% the instantaneous frequency of the simulated signal
if1 = smooth(cumsum(randn(N,1)) ./ Hz, 400, 'loess') ;
if1 = 10 + 6 * if1 ./ max(abs(if1)) ;
if2 = smooth(cumsum(randn(N,1)) ./ Hz, 300, 'loess') ;
if2 = pi + 3 * if2 ./ max(abs(if2)) ;
phi1 = cumsum(if1) / Hz ; 
phi2 = cumsum(if2) / Hz ; 

	%% the simulated signal.
s1 = am1 .* cos(2*pi*phi1) ; 
s2 = am2 .* cos(2*pi*phi2) ; 
clean = s1 + s2 ;

if1(1:500) = nan ;
if2(end-600:end) = nan ;
am1(1:500) = nan ;
am2(end-600:end) = nan ;


	%% add noise (Gaussian white noise)
sigma = 1 ;%sqrt( var(clean)*10.^( -snrdb /10 ) );
noise = random('T',4,N,1) ;
noise = sigma * noise ; 
var(noise)
snrdb = 20 * log10(std(clean)./std(noise)) ;
fprintf(['snrdb = ',num2str(snrdb),'\n']) ;

	%% simulated observed time series
xm = clean + noise ;

Smooth = 0 ;
Hemi = 0 ;

%%%% Multitapering %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        	%% generate the window for short time Fourier transform (STFT)

LEN = 377 ;
    [h, Dh, ~] = hermf(LEN, 2, 10) ;


    clear Q ; clear R ; clear v ;



	[tfr, tfrtic, tfrsqClean, tfrsqtic] = sqSTFTbase(clean, 0, 1/5, 2e-4, 1, h(1,:)', Dh(1,:)', Smooth, 0);
	[~, ~, ~, tfrsqMTclean, tfrsqtic] = ConceFT_STFT(clean, 0, 1/5, 2e-4, 1, LEN, dim, 10, MT, Smooth, Hemi) ;


    [tfrsqMTold, tfrsqMToldAll, tfrsqtic] = sqSTFTmultitaper(xm, 0, 1/5, 2e-4, LEN, 6, 10, Smooth) ;
    tfrsqnoisy = tfrsqMToldAll(:, :, 1) ;

	[~, ~, tfrsq, tfrsqMT, tfrsqtic] = ConceFT_STFT(xm, 0, 1/5, 2e-4, 1, LEN, dim, 10, MT, Smooth, Hemi) ;



	tfrsqClean = tfrsqClean(:, 201:1400) ;
	tfrsqnoisy = tfrsqnoisy(:, 201:1400) ;
	tfrsqMTold = tfrsqMTold(:, 201:1400) ;
	tfrsqMTclean = tfrsqMTclean(:, 201:1400) ;
	tfrsqMT = tfrsqMT(:, 201:1400) ;
	time = time(201:1400) - time(200); xm = xm(201:1400) ; clean = clean(201:1400) ;
	am1 = am1(201:1400) ; am2 = am2(201:1400) ; 
	if1 = if1(201:1400) ; if2 = if2(201:1400) ;
	s1 = s1(201:1400) ; s2 = s2(201:1400) ;
	noise = noise(201:1400) ;



%===========================================
	%% plot the result for comparison. This is Figure 1
hq = figure('Position',[1 scrsz(4)/2 scrsz(3)*2 scrsz(4)]) ;

subplot(3, 2, 1) ;
plot(time, s1, 'color', [.7 .7 .7]) ; axis([0 12 -3 3]) ; set(gca, 'fontsize', 20) ; 
hold on; plot(time, am1, 'k', 'linewidth', 2) ;
set(gca, 'xtick', []) ; ylabel('s_1(t)') ;

subplot(3, 2, 2) ;
plot(time, if1, 'k', 'linewidth', 2) ; set(gca, 'fontsize', 20) ; hold on ;
plot(time, if2, 'k--', 'linewidth', 2) ; axis([0 12 0 18]) ; 
text(4, if1(400)+3, '$$\varphi''_1(t)$$', 'Interpreter','latex', 'fontsize', 24) ;
text(6, if2(600)+2.5, '$$\varphi''_2(t)$$', 'Interpreter','latex', 'fontsize', 24) ; set(gca, 'xtick', []) ;

subplot(3, 2, 3) ;
plot(time, s2, 'color', [.7 .7 .7]) ; axis([0 12 -3 3]) ; set(gca, 'fontsize', 20) ; 
hold on; plot(time, am2, 'k', 'linewidth', 2) ; set(gca, 'xtick', []) ; ylabel('s_2(t)') ;

subplot(3, 2, 4) ;
plot(time, clean, 'k', 'linewidth', 2) ; set(gca, 'fontsize', 20) ;
set(gca, 'xtick', []) ; ylabel('s_1(t)+s_2(t)') ; axis tight ; axis([0 inf -8 10])

subplot(3, 2, 5) ;
plot(time, noise, 'color', [.7 .7 .7]) ; axis tight ; set(gca, 'fontsize', 20) ;
xlabel('Time (sec)') ; ylabel('\xi(t)') ; axis([0 inf -8 10])

subplot(3, 2, 6) ; hold off ;
plot(time, xm, 'color', [.3 .3 .3]) ; axis tight ; set(gca, 'fontsize', 20) ;
xlabel('Time (sec)') ; ylabel('s_1(t)+s_2(t)+\xi(t)') ; axis([0 inf -8 10])
set(hq,'PaperPositionMode','auto');
saveas(hq,['Fig1.eps'],'epsc')




alpha = (tfrsqtic(2)-tfrsqtic(1))*Hz ;
itvPS = zeros(size(tfrsqMTold)) ;

for kk = 1: size(tfrsqMTold, 2)
    if ~isnan(if1(kk)./alpha) ; itvPS(round(if1(kk)./alpha), kk) = am1(kk) ; end
    if ~isnan(if2(kk)./alpha) ; itvPS(round(if2(kk)./alpha), kk) = am2(kk) ; end
end



	%% energy normalization. For the real data, we don't have itvPS,
	%% and the normalization could be done by thresholding or others
tfrsqClean9 = abs((tfrsqClean./( 2 * (tfrtic(2)-tfrtic(1))))).^2 ;
tfrsqnoisy9 = abs((tfrsqnoisy./( 2 * (tfrtic(2)-tfrtic(1))))).^2 ;
tfrsqMTold9 = abs((tfrsqMTold./( 2 * (tfrtic(2)-tfrtic(1))))).^2 ;
tfrsqMT9 = abs((tfrsqMT./( 2 * (tfrtic(2)-tfrtic(1))))).^2 ;
tfrsqMTclean9 = abs((tfrsqMTclean./( 2 * (tfrtic(2)-tfrtic(1))))).^2 ;

%tfrsqClean9 = abs(tfrsqClean).^2 ./ sum(abs(tfrsqClean(:)).^2) ;
%tfrsqnoisy9 = abs(tfrsqnoisy).^2 ./ sum(abs(tfrsqnoisy(:)).^2) ;
%tfrsqMTold9 = abs(tfrsqMTold).^2 ./ sum(abs(tfrsqMTold(:)).^2) ;
%tfrsqMT9 = abs(tfrsqMT).^2 ./ sum(abs(tfrsqMT(:)).^2) ;
%tfrsqMTclean9 = abs(tfrsqMTclean).^2 ./ sum(abs(tfrsqMTclean(:)).^2) ;




%==========================================
	%% determine the common threshold for a fair comparison
LLcwt = zeros(5, 1) ; HHcwt = zeros(5, 1) ;

LLcwt(1) = quantile(log(1+tfrsqClean9(:)),0.002); 
HHcwt(1) = quantile(log(1+tfrsqClean9(:)),0.998);
LLcwt(2) = quantile(log(1+tfrsqnoisy9(:)),0.002); 
HHcwt(2) = quantile(log(1+tfrsqnoisy9(:)),0.998);
LLcwt(3) = quantile(log(1+tfrsqMTold9(:)),0.002); 
HHcwt(3) = quantile(log(1+tfrsqMTold9(:)),0.998);
LLcwt(4) = quantile(log(1+tfrsqMT9(:)),0.002); 
HHcwt(4) = quantile(log(1+tfrsqMT9(:)),0.998);
LLcwt(5) = quantile(log(1+tfrsqMTclean9(:)),0.002); 
HHcwt(5) = quantile(log(1+tfrsqMTclean9(:)),0.998);

minV = 0 %max(LLcwt) ; 
maxV = min(HHcwt) ;



%===========================================
	%% this is Figure 2
hf = figure('Position',[1 scrsz(4)/2 scrsz(3) scrsz(4)]) ;

subplot(9, 2, [1 3 5]) ;
imagesc(time, tfrsqtic*Hz, log(1+abs(tfrsqClean9)), [LLcwt(1) HHcwt(1)]) ; colormap(1-gray) ;
axis xy ; set(gca,'fontsize', 20) ; axis([0 inf 0 20]) ;
set(gca, 'xtick', []); ylabel('Freq (Hz)') ;

subplot(9, 2, [7 9 11]) ;
imagesc(time, tfrsqtic*Hz, log(1+abs(tfrsqnoisy9)), [LLcwt(2) HHcwt(2)]) ; colormap(1-gray) ;
axis xy ; set(gca,'fontsize', 20) ; axis([0 inf 0 20]) ;
set(gca, 'xtick', []); ylabel('Freq (Hz)') ;

subplot(9, 2, [13 15 17]) ;
imagesc(time, tfrsqtic*Hz, log(1+abs(tfrsqMTold9)), [LLcwt(3) HHcwt(3)]) ; colormap(1-gray) ;
axis xy ; set(gca,'fontsize', 20) ; axis([0 inf 0 20]) ;
xlabel('Time (sec)') ; ylabel('Freq (Hz)') ;

subplot(9, 2, [8 10 12]) ;
imagesc(time, tfrsqtic*Hz, log(1+abs(tfrsqMT9)), [LLcwt(4) HHcwt(4)]) ; colormap(1-gray) ;
axis xy ; set(gca,'fontsize', 20) ; axis([0 inf 0 20]) ;
set(gca, 'xtick', []) ; ylabel('Freq (Hz)') ;

subplot(9, 2, [2 4 6]) ;
imagesc(time, tfrsqtic*Hz, log(1+abs(tfrsqMTclean9)), [LLcwt(5) HHcwt(5)]) ; colormap(1-gray) ;
axis xy ; set(gca,'fontsize', 20) ; axis([0 inf 0 20]) ;
set(gca, 'xtick', []) ; ylabel('Freq (Hz)') ;

subplot(9, 2, [14 16 18]) ;
imagesc(time, tfrsqtic*Hz, log(1+abs(tfrsqMT9)), [LLcwt(4) HHcwt(4)]) ; colormap(1-gray) ;
axis xy ; set(gca,'fontsize', 20) ; axis([0 inf 0 20]) ;
hold on; plot(time, if1, 'r', 'linewidth', 3) ;
plot(time, if2, 'b', 'linewidth', 3) ;
xlabel('Time (sec)') ; ylabel('Freq (Hz)') ;

set(hf,'PaperPositionMode','auto');
saveas(hf,['Fig2.eps'],'epsc')

