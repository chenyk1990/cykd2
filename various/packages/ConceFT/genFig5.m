clear ;
% close all ;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Plotting clean and noisy signals
%==================
scrsz = get(0,'ScreenSize');
Hz = 64;
T = 80;
time = [1/Hz:1/Hz:T]' ;
N = length(time) ;
freqLow = 0;
freqHigh = 20;
alpha = 0.01;
ftsz = 22;
%==================

ExampleID = 2;
seeds = [1 111 23400];
dnSampleID = 1:1:(Hz*T);


NoiseID = 1;
initstate(seeds(ExampleID));
tau = 10 ;
clear sigma
snrdb = 0;
loadExamples
time = time(dnSampleID);
if1 = if1(dnSampleID);
if2 = if2(dnSampleID);
am1 = am1(dnSampleID);
am2 = am2(dnSampleID);

scaling = 1 / alpha;
trueIF = sparse(round(if1((TN1+1):N)*scaling),(TN1+1):N,am1((TN1+1):N).^2,(freqHigh-freqLow)*scaling,N)...
    + sparse(round(if2( (1:(N-TN2-1)))*scaling ),1:(N-TN2-1),am2(1:(N-TN2-1)).^2,(freqHigh-freqLow)*scaling,N);


%========================
	%% plot Figure 5
hf = figure('Position',[1 scrsz(4)/2 scrsz(3) scrsz(4)*2/3]) ;

subplot(4, 4, [1 5]) ;
plot(time,if1,'k','LineWidth',2),    hold on
plot(time,if2,'k','LineWidth',2)
axis([0 T freqLow freqHigh])
set(gca,'fontsize',ftsz)
%xlabel('time') ; ylabel('frequency') ; %title('ground truth')
%saveas(hf,'otDemo_true','epsc')

%%
maxv1 = 3.4;
tmp = linspace(0,maxv1,300);
offset1 = [tmp fliplr(tmp)];
initID = 2500;
idx = initID:(initID+length(offset1)-1);
if1_test = if1; if2_test = if2;
if1_test(TN1+idx) = if1(TN1+idx) + offset1';
if1_test = smooth(if1_test, 200, 'loess');
if1_test(1:TN1) = nan;

%hf = figure;
subplot(4, 4, [2 6]) ;
plot(time,if1_test,'color', [.7 .7 .7], 'LineWidth',3)
hold on
plot(time,if2_test,'color', [.7 .7 .7], 'LineWidth',3)
plot(time,if1,'k','LineWidth', 3),
plot(time,if2,'k','LineWidth', 3)
axis([0 T freqLow freqHigh])
set(gca,'fontsize',ftsz) ; %xlabel('time') ; ylabel('frequency')


testIF1 = sparse(round(if1_test((TN1+1):N)*scaling),(TN1+1):N,am1((TN1+1):N).^2,(freqHigh-freqLow)*scaling,N)...
    + sparse(round(if2_test( (1:(N-TN2-1)))*scaling ),1:(N-TN2-1),am2(1:(N-TN2-1)).^2,(freqHigh-freqLow)*scaling,N);

ot = slicedOT(trueIF, testIF1) * 100;
% ot = slicedOT(trueIF, trueIF) * 100;
text(40, 18, sprintf('OT=%.2f',ot), 'fontsize', 22)
%saveas(hf,'otDemo1','epsc')
%%
offset2 = 2;
if1_test = if1 + randn(size(if1))*offset2; if2_test = if2 + randn(size(if2))*offset2;
if1_test = smooth(if1_test, 200, 'loess');
if1_test(1:TN1) = nan;
if2_test = smooth(if2_test, 200, 'loess');
if2_test(end-TN2:end) = nan ;

subplot(4, 4, [3 7]) ;
plot(time,if1_test,'color', [.7 .7 .7], 'LineWidth', 3),
hold on
plot(time,if2_test,'color', [.7 .7 .7], 'LineWidth', 3),
plot(time,if1,'k','LineWidth', 3),
plot(time,if2,'k','LineWidth', 3)
axis([0 T freqLow freqHigh])
set(gca,'fontsize',ftsz) ; %xlabel('time') ; ylabel('frequency')

testIF1 = sparse(round(if1_test((TN1+1):N)*scaling),(TN1+1):N,am1((TN1+1):N).^2,(freqHigh-freqLow)*scaling,N)...
    + sparse(round(if2_test( (1:(N-TN2-1)))*scaling ),1:(N-TN2-1),am2(1:(N-TN2-1)).^2,(freqHigh-freqLow)*scaling,N);

ot = slicedOT(trueIF, testIF1) * 100;
%title(sprintf('demo2,OT=%.2f',ot))
text(40, 18, sprintf('OT=%.2f',ot), 'fontsize', 22)
%saveas(hf,'otDemo2','epsc')

%%
offset2 = .2;
if1_test = if1 + 2*offset2; if2_test = if2 + 2*offset2;
if2_test(round(5120/3):end) = if2_test(round(5120/3):end) - 8*offset2 ;
if1_test = smooth(if1_test, 200, 'loess');
if1_test(1:TN1) = nan;
if2_test = smooth(if2_test, 200, 'loess');
if2_test(end-TN2:end) = nan ;

%hf = figure;
subplot(4, 4, [4 8]) ;
plot(time,if1_test,'color', [.7 .7 .7], 'LineWidth', 3),
hold on
plot(time,if2_test,'color', [.7 .7 .7], 'LineWidth', 3),
plot(time,if1,'k','LineWidth', 3),
plot(time,if2,'k','LineWidth', 3)
axis([0 T freqLow freqHigh])
set(gca,'fontsize',ftsz) ; %xlabel('time') ; ylabel('frequency')

testIF1 = sparse(round(if1_test((TN1+1):N)*scaling),(TN1+1):N,am1((TN1+1):N).^2,(freqHigh-freqLow)*scaling,N)...
    + sparse(round(if2_test( (1:(N-TN2-1)))*scaling ),1:(N-TN2-1),am2(1:(N-TN2-1)).^2,(freqHigh-freqLow)*scaling,N);

ot = slicedOT(trueIF, testIF1) * 100;
%title(sprintf('demo3,OT,=%.2f',ot))
text(40, 18, sprintf('OT=%.2f',ot), 'fontsize', 22)

%saveas(hf,'otDemo3','epsc')



% ======================
	%% amplitude difference

tfrsqtic = linspace(0, Hz/2, 1000) ;
itvPS = zeros(length(tfrsqtic), length(am1)) ;
alpha = tfrsqtic(2) ;
for ii = 1: length(am1)
    if ~isnan(if1(ii)./alpha) ; itvPS(round(if1(ii)./alpha), ii) = am1(ii) ; end
    if ~isnan(if1(ii)./alpha) ; itvPS(round(if1(ii)./alpha)+1, ii) = am1(ii) ; end
    if ~isnan(if1(ii)./alpha) ; itvPS(round(if1(ii)./alpha)+2, ii) = am1(ii) ; end
    if ~isnan(if1(ii)./alpha) ; itvPS(round(if1(ii)./alpha)+3, ii) = am1(ii) ; end
    if ~isnan(if1(ii)./alpha) ; itvPS(round(if1(ii)./alpha)-1, ii) = am1(ii) ; end
    if ~isnan(if1(ii)./alpha) ; itvPS(round(if1(ii)./alpha)-2, ii) = am1(ii) ; end
    if ~isnan(if1(ii)./alpha) ; itvPS(round(if1(ii)./alpha)-3, ii) = am1(ii) ; end
    if ~isnan(if2(ii)./alpha) ; itvPS(round(if2(ii)./alpha), ii) = am2(ii) ; end
    if ~isnan(if2(ii)./alpha) ; itvPS(round(if2(ii)./alpha)+1, ii) = am2(ii) ; end
    if ~isnan(if2(ii)./alpha) ; itvPS(round(if2(ii)./alpha)+2, ii) = am2(ii) ; end
    if ~isnan(if2(ii)./alpha) ; itvPS(round(if2(ii)./alpha)+3, ii) = am2(ii) ; end
    if ~isnan(if2(ii)./alpha) ; itvPS(round(if2(ii)./alpha)-1, ii) = am2(ii) ; end
    if ~isnan(if2(ii)./alpha) ; itvPS(round(if2(ii)./alpha)-2, ii) = am2(ii) ; end
    if ~isnan(if2(ii)./alpha) ; itvPS(round(if2(ii)./alpha)-3, ii) = am2(ii) ; end
end



%%
mksz = 2;
%hf = figure;
subplot(4, 4, [10 14]) ;
imageSQ(time, tfrsqtic, itvPS.^2) ; colormap(1-gray) ;
%scatter(time,if1,mksz,am1),hold on
%scatter(time,if2,mksz,am2)
%axis([0 T freqLow freqHigh])
%caxis([1 3])
colorbar('North') 
axis([0 T freqLow freqHigh])
set(gca,'fontsize',ftsz) ; %xlabel('time') ; ylabel('frequency')
%title('ground truth with amplitude')
%box on
%saveas(hf,'otDemo4_true','epsc')

%%
initstate(1)
offset4 = 5.2;
am1_test = am1 + 3/4 ;
am2_test = am2 - 1/4 ;% offset4*randn(size(am2));
%am1_test = smooth(am1_test, 800, 'loess');
%am2_test = smooth(am2_test, 800, 'loess');
am1_test(1:TN1) = 0 ;
am2_test(end-TN2:end) = 0 ;
testIF1 = sparse(round(if1((TN1+1):N)*scaling),(TN1+1):N,am1_test((TN1+1):N).^2,(freqHigh-freqLow)*scaling,N)...
    + sparse(round(if2( (1:(N-TN2-1)))*scaling ),1:(N-TN2-1),am2_test(1:(N-TN2-1)).^2,(freqHigh-freqLow)*scaling,N);
ot = slicedOT(trueIF, testIF1) * 100;


itvPS = zeros(length(tfrsqtic), length(am1)) ;
for ii = 1: length(am1)
    if ~isnan(if1(ii)./alpha) ; itvPS(round(if1(ii)./alpha), ii) = am1_test(ii) ; end
    if ~isnan(if1(ii)./alpha) ; itvPS(round(if1(ii)./alpha)+1, ii) = am1_test(ii) ; end
    if ~isnan(if1(ii)./alpha) ; itvPS(round(if1(ii)./alpha)+2, ii) = am1_test(ii) ; end
    if ~isnan(if1(ii)./alpha) ; itvPS(round(if1(ii)./alpha)+3, ii) = am1_test(ii) ; end
    if ~isnan(if1(ii)./alpha) ; itvPS(round(if1(ii)./alpha)-1, ii) = am1_test(ii) ; end
    if ~isnan(if1(ii)./alpha) ; itvPS(round(if1(ii)./alpha)-2, ii) = am1_test(ii) ; end
    if ~isnan(if1(ii)./alpha) ; itvPS(round(if1(ii)./alpha)-3, ii) = am1_test(ii) ; end
    if ~isnan(if2(ii)./alpha) ; itvPS(round(if2(ii)./alpha), ii) = am2_test(ii) ; end
    if ~isnan(if2(ii)./alpha) ; itvPS(round(if2(ii)./alpha)+1, ii) = am2_test(ii) ; end
    if ~isnan(if2(ii)./alpha) ; itvPS(round(if2(ii)./alpha)+2, ii) = am2_test(ii) ; end
    if ~isnan(if2(ii)./alpha) ; itvPS(round(if2(ii)./alpha)+3, ii) = am2_test(ii) ; end
    if ~isnan(if2(ii)./alpha) ; itvPS(round(if2(ii)./alpha)-1, ii) = am2_test(ii) ; end
    if ~isnan(if2(ii)./alpha) ; itvPS(round(if2(ii)./alpha)-2, ii) = am2_test(ii) ; end
    if ~isnan(if2(ii)./alpha) ; itvPS(round(if2(ii)./alpha)-3, ii) = am2_test(ii) ; end
end

%%
%hf = figure;
subplot(4, 4, [11 15]) ;
imageSQ(time, tfrsqtic, itvPS.^2) ; colormap(1-gray) ;
%scatter(time,if1,mksz,am1),hold on
%scatter(time,if2,mksz,am2)
%axis([0 T freqLow freqHigh])
%caxis([1 3])
colorbar('North')
%scatter(time,if1,mksz,am1_test),hold on
%scatter(time,if2,mksz,am2_test)
axis([0 T freqLow freqHigh])
%caxis([1 3])
%     colormap(gray)
%colorbar
set(gca,'fontsize',ftsz) ; %xlabel('time') ; ylabel('frequency') ;
text(40, 14, sprintf('OT=%.2f',ot), 'fontsize', 22)
%box on
%saveas(hf,'otDemo4','epsc')

set(hf,'PaperPositionMode','auto');
saveas(hf,['Fig5.eps'],'epsc')

