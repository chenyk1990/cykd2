clear ; close all ;

addpath('/Users/hautiengwu/Dropbox/ConCeft/SubmissionCode/Conceft') ;


%%
scrsz = get(0,'ScreenSize');
Hz = 160;
T = 60;
time = [1/Hz:1/Hz:T]' ;
N = length(time) ;
freqLow = 2;
freqHigh = 20;
alpha = 0.01;
ftsz = 22;
%==================
seeds = [111 4];
taus = [10 8];
%==================
dim = 4 ; 
MT = 20 ;
Smooth = 0 ;
Hemi = 1 ;
%==================
opts.motherwavelet = 'morse-c' ;
opts.k = 0;
opts.dim = 2;
opts.beta = 30;
opts.gam = 9;
opts.rrnd = 0;


%% obtaining the range from CWT results
for ii = 1:3
    for jj = 1:2
        ExampleID = jj;
        tau = taus(jj);
        initstate(seeds(ExampleID));
        NoiseID = ii;
        clear sigma
        snrdb = 0;
        loadExamples
        
        [~, tmp, tfrsqticCWT] = ConceFT_CWT(time, xm, freqLow, freqHigh, alpha, MT, opts, Smooth, Hemi) ;

		eval(['tmpCWT',num2str(ii),num2str(jj),' = tmp ;']) ;

    end
end


%% obtaining the range from STFT results
for ii = 1:3;
	for jj = 1:2
	    ExampleID = jj;
	    tau = taus(jj);
	    initstate(seeds(ExampleID));
	    NoiseID = ii;
	    clear sigma
	    snrdb = 0;
	    loadExamples
    
	    [tfrsqOrig, tfrsqMToldAll, tfrsqtic] = sqSTFTmultitaper(xm, 0, 0.5/4, 6.25e-5, 601, 6, 6, Smooth) ;
	    tfrsqSST = tfrsqMToldAll(:, :, 1) ;
	    [~, ~, ~, tmp, tfrsqticSTFT] = ConceFT_STFT(xm, 0, 0.5/4, 6.25e-5, 1, 601, dim, 6, MT, Smooth, Hemi) ;

		eval(['tfrsqOrig',num2str(ii),num2str(jj),' = tfrsqOrig ;']) ;
		eval(['tfrsqSST',num2str(ii),num2str(jj),' = tfrsqSST ;']) ;
		eval(['tmpSTFT',num2str(ii),num2str(jj),' = tmp ;']) ;

	end    
end




%=================================================
	%% generating Figure S8

hf = figure('OuterPosition',[1 scrsz(4)/2 scrsz(3) scrsz(4)*2/3]);
for ii = 1:3
    for jj=1:2
        ExampleID = jj;
        tau = taus(jj);
        initstate(seeds(ExampleID));
        NoiseID = ii;
        clear sigma
        snrdb = 0;
        loadExamples
        
        eval(['tmp = tmpCWT',num2str(ii),num2str(jj),' ;']) ;

        scaling = 1 / alpha;
        if1 = if1 - freqLow; if2 = if2 - freqLow;
        trueIF = sparse(round(if1((TN1+1):N)*scaling),(TN1+1):N,am1((TN1+1):N).^2,(freqHigh-freqLow)*scaling,N)...
            + sparse(round(if2( (1:(N-TN2-1)))*scaling ),1:(N-TN2-1),am2(1:(N-TN2-1)).^2,(freqHigh-freqLow)*scaling,N);
        trueIF = full(trueIF);
 
		OT = slicedOT(trueIF, abs(tmp).^2)*100;       
        
        itvPS = zeros(size(tmp)) ;
        
        for kk = 1: size(itvPS, 2)
            
            if ~isnan(if1(kk)./alpha) ;
                itvPS(round(if1(kk)./alpha)-12:round(if1(kk)./alpha)+12, kk) = am1(kk) ;
            end
            
            if ~isnan(if2(kk)./alpha) ;
                itvPS(round(if2(kk)./alpha)-12:round(if2(kk)./alpha)+12, kk) = am2(kk) ;
            end
            
        end
        
        subplot(2, 4, (jj-1)*4+1) ; hold off
        imagesc(time, tfrsqticCWT, log(1+abs(itvPS).^2)) ;
        axis xy ; colormap(1-gray) ; set(gca,'fontsize', 20) ;
        axis([0 inf freqLow freqHigh]) ; xlabel('Time (sec)') ; ylabel('Frequency (Hz)') ;
        
        YYnm = abs(tmp).^2;
        subplot(2, 4, (jj-1)*4+ii+1) ;
        imagesc(time, tfrsqticCWT, log(1+YYnm),[0 quantile(log(1+YYnm(:)),.99)] ); colormap(1-gray) ;
        axis xy
        set(gca,'fontsize', 20) ; axis([0 inf freqLow freqHigh]) ;
        xlabel('Time (sec)') ; ylabel('Frequency (Hz)') ;
        title(sprintf('OT=%.2f',OT(1)))
        
    end
end
set(hf,'PaperPositionMode','auto');
saveas(hf,'FigS8','eps')







%================================================
	%% generating Figure S9 and S10
for jj=1:2

    hf = figure('OuterPosition',[1 scrsz(4)/2 scrsz(3) scrsz(4)]);
    
    for ii = 1:3
        ExampleID = jj;
        tau = taus(jj);
        initstate(seeds(ExampleID));
        
        NoiseID = ii; clear sigma
        snrdb = 0;
        loadExamples
	
    	scaling = 1 / (tfrsqticSTFT(2)-tfrsqticSTFT(1)) / Hz ;
        if1 = if1 - freqLow; if2 = if2 - freqLow;
        trueIF = sparse(round(if1((TN1+1):N)*scaling),(TN1+1):N,am1((TN1+1):N).^2,length(tfrsqtic)-200,N)...
            + sparse(round(if2( (1:(N-TN2-1)))*scaling ),1:(N-TN2-1),am2(1:(N-TN2-1)).^2,length(tfrsqtic)-200,N);
        trueIF = full(trueIF);

        eval(['tfrsqOrig = tfrsqOrig',num2str(ii),num2str(jj),' ;']) ;
        eval(['tfrsqSST = tfrsqSST',num2str(ii),num2str(jj),' ;']) ;
        eval(['tmp = tmpSTFT',num2str(ii),num2str(jj),' ;']) ;


        YY1 = abs(tmp(201:2000,:)).^2;   % ConceFT
        YY2 = abs(tfrsqOrig(201:2000,:)).^2;  %oldMT
        YY3 = abs(tfrsqSST(201:2000,:)).^2;  %SST
        
        itvPS = zeros(size(YY1)) ;
        
        for kk = 1: size(itvPS, 2)           
            if ~isnan(if1(kk)./alpha) ;
                itvPS(round(if1(kk)./alpha)-12:round(if1(kk)./alpha)+12, kk) = am1(kk) ;
            end
            
            if ~isnan(if2(kk)./alpha) ;
                itvPS(round(if2(kk)./alpha)-12:round(if2(kk)./alpha)+12, kk) = am2(kk) ;
            end           
        end
        
        subplot(3, 4, 5) ; hold off
        imagesc(time, tfrsqticSTFT(201:2000)*Hz, log(1+abs(itvPS).^2)) ; colormap(1-gray) ;
        axis xy ; set(gca,'fontsize', 20) ;
        axis([0 inf freqLow freqHigh]) ; xlabel('Time (sec)') ; ylabel('Frequency (Hz)') ;
               
        ot = slicedOT(trueIF,YY3);
        subplot(3, 4, (ii-1)*4+2) ;
        imagesc(time, tfrsqticSTFT(201:2000)*Hz, log(1+YY3), [0 quantile(log(1+YY3(:)),.99)]) ; colormap(1-gray) ;
        axis xy ; set(gca,'fontsize', 20) ; axis([0 inf freqLow freqHigh]) ;
        xlabel('Time (sec)') ; ylabel('Frequency (Hz)') ;
        title(sprintf('OT=%.2f',ot*100))
        
        ot = slicedOT(trueIF,YY2);
        subplot(3, 4, (ii-1)*4+3) ;
        imagesc(time, tfrsqticSTFT(201:2000)*Hz, log(1+YY2),[0 quantile(log(1+YY2(:)),.99)]) ; colormap(1-gray) ;
        axis xy ; set(gca,'fontsize', 20) ; axis([0 inf freqLow freqHigh]) ;
        xlabel('Time (sec)') ; ylabel('Frequency (Hz)') ;
        title(sprintf('OT=%.2f',ot*100))
        
        ot = slicedOT(trueIF,YY1);
        subplot(3, 4, (ii-1)*4+4) ;
        imagesc(time, tfrsqticSTFT(201:2000)*Hz, log(1+YY1),[0 quantile(log(1+YY1(:)),.99)]) ; colormap(1-gray) ;
        axis xy ; set(gca,'fontsize', 20) ; axis([0 inf freqLow freqHigh]) ;
        xlabel('Time (sec)') ; ylabel('Frequency (Hz)') ;
        title(sprintf('OT=%.2f',ot*100))
        
    end
    set(hf,'PaperPositionMode','auto');
    saveas(hf,['FigS',num2str(jj+8)],'eps')
    
end
