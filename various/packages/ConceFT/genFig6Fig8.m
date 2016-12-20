% clear ;
% close all ;
PC = 2;
switch PC
    case 0
        addpath(genpath('/gtmp/YWang/Ingrid/synchroSqueezing'))
    case 1
        addpath(genpath('F:\Dropbox\Research\Ingrid\Component Detection\codes'))
    case 2
        addpath(genpath('/dscrgrps/harerlab/yw112/Conceft/'))
end
%% fix the seed
scrsz = get(0,'ScreenSize');
ftsz = 22;
Hz = 64;
T = 80;
time = [1/Hz:1/Hz:T]' ;
N = length(time) ;
freqLow = 0;
freqHigh = 20;
alpha = 0.01;
seeds = [1 111 23400];
%==================
opts.motherwavelet = 'morse-c' ;
opts.k = 0;
opts.beta = 30;
opts.gam = 3;
opts.dim = 8;
opts.rrnd = 0;

MT = 1000;

ExampleID = jj;
NoiseID = ii;
initstate(seeds(ExampleID));
snrdb = 0;
clear sigma
loadExamples

rng('default')
rng('shuffle')
disp('random seed')
%%
scaling = 1 / alpha;
trueIF = sparse(round(if1((TN1+1):N)*scaling),(TN1+1):N,abs(am1((TN1+1):N)).^2,(freqHigh-freqLow)*scaling,N)...
    + sparse(round(if2( (1:(N-TN2-1)))*scaling ),1:(N-TN2-1),abs(am2(1:(N-TN2-1))).^2,(freqHigh-freqLow)*scaling,N);
trueIF = full(trueIF);
%%
OT_N = zeros(MT,1);
tfrsq = 0;

for kk = 1: MT
    if mod(kk,100)==0
        fprintf('kk =%d\n',kk)
    end
    [~, tfrsqX, ~, tfrsqtic] = sqCWT(time, xm, 1, 2, 32, 1e-8, alpha,freqLow,freqHigh, opts);
    tfrsq = tfrsq + tfrsqX ;
    tmp = tfrsq ./ kk;
    OT_N(kk) = slicedOT(trueIF, abs(tmp').^2)*100;
end

save(['/netscratch/yw112/N_CWT_rrnd0_rs/N_exID',num2str(ExampleID),'_nsID',num2str(NoiseID),...
    '_arrayID',getenv('SLURM_ARRAY_TASK_ID')], 'OT_N')

%%
cd('/netscratch/yw112/N_CWT_rrnd0_rs/new')
for jj = 1:2
    for ii = 1:3
        
        ExampleID = jj;
        NoiseID = ii;

        myDir = dir(['N_exID',num2str(ExampleID),'_nsID',num2str(NoiseID),'*.mat']);
        ot = zeros(1000,length(myDir));
        for tt = 1:length(myDir)
            load(myDir(tt).name)
            ot(:,tt) = OT_N;
        end
        ot_avg = mean(ot,2);
        ot_std = std(ot,0,2);
        
        save(['/dscrgrps/harerlab/yw112/Conceft/Conceft/workingfiles/N/N1000_exID',...
            num2str(ExampleID),'_nsID',num2str(NoiseID)], 'ot_avg','ot_std','ot')
    end
end

%%
scrsz = get(0,'ScreenSize');
ftsz = 22;

for jj=1:2
    for ii = 1:3
        
        ExampleID = jj;
        NoiseID = ii;
        
        load(['N1000_exID',num2str(ExampleID),'_nsID',num2str(NoiseID)])
        xx = 1:length(ot_avg);
        data = ot_avg;
        dataSTD = ot_std;
        
        hf = figure;
%         hf = figure('OuterPosition',[1 1 scrsz(3) scrsz(4)]) ;

        xx = xx';
%         xx = log(xx)';
%         xtck = 100:200:1000;
%         xtck = log(xtck);
        xmax = max(xx);
        fill([xx;flipud(xx)],[data-dataSTD;flipud(data+dataSTD)],[.9 .9 .9],'linestyle','none');
        hold on
        plot(xx,data,'b-','linewidth',2)
        set(gca,'fontsize',ftsz)
        if jj==1
            axis([0 xmax 13 15]),
        else
            axis([0 xmax 11 13])
        end
%         axis([-1 inf 10 18])
        xlabel('the number of random projections')
        ylabel('OT')
        saveas(hf,['OT_N1000_ExID',num2str(ExampleID),'_NsID',num2str(NoiseID)],'epsc')
    end
end


%%
% scrsz = get(0,'ScreenSize');
% ftsz = 22;
% 
% for jj=1:2
%     for ii = 1:3
%         
%         ExampleID = jj;
%         NoiseID = ii;
%         
%         load(['N_exID',num2str(ExampleID),'_nsID',num2str(NoiseID)])
%         data = ot_avg;
%         dataSTD = ot_std;
%         
%         hf = figure('OuterPosition',[1 1 scrsz(3) scrsz(4)]) ;
%         xx = [1 2 10 20 50:50:1000];
% %         xx = xx';
%         xx = log(xx)';
%         xtck = 100:200:1000;
%         xtck = log(xtck);
%         xmax = 1000;
%         xmax = log(xmax);
%         fill([xx;flipud(xx)],[data-dataSTD;flipud(data+dataSTD)],[.9 .9 .9],'linestyle','none');
%         hold on
%         plot(xx,data,'b-','linewidth',2)
%         set(gca,'fontsize',ftsz)
%         axis([0 xmax 10 18]),
% %         axis([-1 inf 10 18])
%         xlabel('the number of random projections')
%         ylabel('OT')
%         %         saveas(hf,['OT_N_ExID',num2str(ExampleID),'_NsID',num2str(NoiseID)],'epsc')
%     end
% end
