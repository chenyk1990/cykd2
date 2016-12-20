clear ; close all ;

addpath('/Users/hautiengwu/Dropbox/ConCeft/SubmissionCode/Conceft') ;
scrsz = get(groot,'ScreenSize');

	%% generate the simulated data
	%% she sampling time (100Hz sampling rate)
	%% high sampling rate to avoid sampling issue
Hz = 160 ;
time = [1/Hz:1/Hz:70]' ;
N = length(time) ;

snrdb = 0 ;
DDD = 1800*9600*5 ;

	%% the amplitude modulation of the simulated signal
am1 = 1 + 0.3*cos((pi*(time-10)./20).^2) ;
am2 = 0.4 + 0.9*sin(pi*time/60).^2 ;
am3 = 1.2*ones(size(time)) ;
am1(1:10*Hz) = 0 ; am1(48*Hz+1:end) = 0 ; am3(1:15*Hz) = 0 ;



    %% the instantaneous frequency of the simulated signal

	%% the simulated signal.
s1 = am1 .* cos(2*pi*(pi/3+5*time+time.^2/50)) ; 
s2 = am2 .* cos(2*pi*(12*time+sin(pi*time/6))) ; 
s3 = am3 .* cos(2*pi*(17*time + (time-35).^3/800)) ;
clean = s1 + s2 + s3 ;

if1 = 5 + time/25 ;
if2 = 12 + cos(pi*time/6) * pi / 6 ;
if3 = 17 + 3*(time-35).^2 / 800 ;
if1(1:10*Hz) = nan ; if1(48*Hz+1:end) = nan ; if3(1:15*Hz) = nan ;

sigma = sqrt( var(clean)*10.^( -snrdb /10 ) );

	%% add noise (Gaussian white noise)
noise1 = randn(N, 1) ;
xm1 = clean  + sigma * noise1 ;

dis = random('t',4,N, 1) ;
e = armaxfilter_simulate(dis, .5, 1, .5, 1, -.5) ;
noise2 = e ./ std(e) ;
xm2 = clean  + sigma * noise2 ;

noise3 = random('poiss',1,N,1);
xm3 = clean  + sigma * noise3 ;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

MT = 20;
opts.k = 0;
opts.dim = 2;
opts.beta = 30;
opts.gam = 9;
opts.rrnd = 0;
opts.motherwavelet = 'morse-c' ;
Smooth = 0 ;
Hemi = 1 ;

[tfrsq1, tfrsqMT1, tfrsqtic] = ConceFT_CWT(time, xm1, 2, 20, 0.01, MT, opts, Smooth, Hemi) ;
[tfrsq2, tfrsqMT2, tfrsqtic] = ConceFT_CWT(time, xm2, 2, 20, 0.01, MT, opts, Smooth, Hemi) ;
[tfrsq3, tfrsqMT3, tfrsqtic] = ConceFT_CWT(time, xm3, 2, 20, 0.01, MT, opts, Smooth, Hemi) ;

tfrsqMT1 = tfrsqMT1(:, Hz*5+1: length(xm1)-Hz*5) ;
tfrsqMT2 = tfrsqMT2(:, Hz*5+1: length(xm2)-Hz*5) ;
tfrsqMT3 = tfrsqMT3(:, Hz*5+1: length(xm3)-Hz*5) ;




%=======================================
	%% determined overall threshold
alpha = tfrsqtic(2) - tfrsqtic(1) ;
itvPS = zeros(size(tfrsqMT1)) ;
trueIF = zeros(size(tfrsqMT1)) ;

if1 = if1 - 2 ; if2 = if2 - 2 ; if3 = if3 - 2 ;
for ii = 1: size(tfrsqMT1, 2)

    if ~isnan(if1(Hz*5+ii)) ;
        itvPS(round(if1(Hz*5+ii)./alpha)-12:round(if1(Hz*5+ii)./alpha)+12, ii) = am1(Hz*5+ii) ;
        trueIF(round(if1(Hz*5+ii)./alpha), ii) = am1(Hz*5+ii) ;
    end

    if ~isnan(if2(Hz*5+ii)) ;
        itvPS(round(if2(Hz*5+ii)./alpha)-12:round(if2(Hz*5+ii)./alpha)+12, ii) = am2(Hz*5+ii) ;
        trueIF(round(if2(Hz*5+ii)./alpha), ii) = am2(Hz*5+ii) ;
    end

    if ~isnan(if3(Hz*5+ii)) ;
        itvPS(round(if3(Hz*5+ii)./alpha)-12:round(if3(Hz*5+ii)./alpha)+12, ii) = am3(Hz*5+ii) ;
        trueIF(round(if3(Hz*5+ii)./alpha), ii) = am3(Hz*5+ii) ;
    end

end


LLcwt = zeros(3, 1) ; HHcwt = zeros(3, 1) ;

YY = abs(tfrsqMT1).^2 ;
YY1 = DDD*YY / sum(YY(:)) ; %* sum(trueIF(:)) ;
LLcwt(1) = quantile(log(1+YY1(:)),0.002); HHcwt(1) = quantile(log(1+YY1(:)),0.998);

YY = abs(tfrsqMT2).^2 ;
YY2 = DDD*YY / sum(YY(:)) ; %* sum(trueIF(:)) ;
LLcwt(2) = quantile(log(1+YY2(:)),0.002); HHcwt(2) = quantile(log(1+YY2(:)),0.998);

YY = abs(tfrsqMT3).^2 ;
YY3 = DDD*YY / sum(YY(:)) ; %* sum(trueIF(:)) ;
LLcwt(3) = quantile(log(1+YY3(:)),0.002); HHcwt(3) = quantile(log(1+YY3(:)),0.998);

minV = 0 ; %max(LLcwt) ; 
maxV = 5.718 ; %min(HHcwt) ;






%===========================================
	%% this is Figure 9
hf = figure('Position',[1 scrsz(4)/2 scrsz(3) scrsz(4)/3]) ;


subplot(1, 4, 1) ; hold off
imagesc(time, tfrsqtic, log(1+abs(itvPS).^2)) ; colormap(1-gray) ;
axis xy ; set(gca,'fontsize', 20) ;
xlabel('Time (sec)') ; ylabel('Frequency (Hz)') ;

subplot(1, 4, 2) ;
imagesc(time, tfrsqtic, log(1+YY1), [minV maxV]) ; colormap(1-gray) ;
axis xy ; set(gca,'fontsize', 20) ;
xlabel('Time (sec)') ; ylabel('Frequency (Hz)') ;

subplot(1, 4, 3) ;
imagesc(time, tfrsqtic, log(1+YY2), [minV maxV]) ; colormap(1-gray) ;
axis xy ; set(gca,'fontsize', 20) ;
xlabel('Time (sec)') ; ylabel('Frequency (Hz)') ;

subplot(1, 4, 4) ;
imagesc(time, tfrsqtic, log(1+YY3), [minV maxV]) ; colormap(1-gray) ;
axis xy ; set(gca,'fontsize', 20) ;
xlabel('Time (sec)') ; ylabel('Frequency (Hz)') ;

set(hf,'PaperPositionMode','auto');
saveas(hf,['Fig9.eps'],'eps')

