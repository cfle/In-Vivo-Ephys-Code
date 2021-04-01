%Speed-binned LFP analysis using the SpeedSnip method. Suitable for one
%channel at a time.

function [EventData] = BurstDetectionWYW(Folder,Channel)

%First selects the folder with the data, and locates all the LFP data
folder = Folder;
folder = dir(folder);
foldername = folder.folder;
subfolders = {folder.name};
subfolders = subfolders(contains(subfolders,'Day'));
subfolders = strcat(foldername,'\',subfolders);
sessionnames = ["Day1a","Day1b","Day2a","Day2b","Day3a","Day3b","Day4a","Day4b","Day5a","Day5b"];
cd(foldername);

%Starts a megaloop to analyse and cycle through each session in an animals
%folder.
for ind0 = 1:length(subfolders)
    
    files = dir(subfolders{ind0});
    files = files(arrayfun(@(x) ~strcmp(x.name(1),'.'),files));
    channels = {files.name};
    channels = channels(contains(channels,Channel));
    channels = strcat(subfolders{ind0},'\',channels);
    channelname = string(channels{1}(end-14:end-11));
    
    if sum(contains(["CH13","CH14","CH19","CH20"],channelname)) == 1
        location = "RSA";
    elseif sum(contains(["CH45","CH46","CH51","CH52"],channelname)) == 1
        location = "MEC";
    end
    
    animalnumber = Folder(end-1:end);
    animal = strcat('Mouse',animalnumber);
    session = files.folder;
    session = session(end-4:end);
    
    [~,~,extension] = fileparts(channels{1});
    
    if strcmp(extension,'.continuous')
        %Performs tracking, calculates speed and total distance travelled in m
        [x,y,~,speed_out,distance] = TrackingNovelFamiliar(subfolders{ind0},1,'two');
        tsf = fillmissing(speed_out.ts,'linear');
        speed = speed_out.speed;
        %Imports the .continuous data into MatLab
        signals = [];
        [signals(:,1),~,info] = load_open_ephys_data_faster(channels{1});
    elseif strcmp(extension,'.mat')
        [x,y,ts,speed_out,distance] = TrackingNovelFamiliarSplitSession(subfolders{ind0});
        tsf = fillmissing(speed_out.ts,'linear');
        speed = speed_out.speed;
        %Imports the .continuous data into MatLab
        load(channels{1});
        signals(:,1) = sig;
        info.header.sampleRate = 30000;
    end
    
    %Cuts the signals to the experiment duration ie 15minutes in this case
    ExpDurationMinutes = 15;
    Fs = info.header.sampleRate;
    ExpDurationFrames = (ExpDurationMinutes*60)*Fs;
    ExpDurationSeconds = ExpDurationMinutes*60;
    signals = signals(1:ExpDurationFrames,:);
    
    %Downsamples the signals and timestamps, to reduce computation time,
    %and detrends to remove any baseline fluctuations
    downsampleX = 10;
    Fs = Fs/downsampleX;
    ExpDurationFrames = ExpDurationFrames/downsampleX;
    signalsds = downsample(signals,downsampleX);
    signalsdsdt = detrend(signalsds);
    signalmagnitude = abs(hilbert(signalsdsdt));
    %     noisethreshold = prctile(signalmagnitude,99.99);
    %     noisyframes = signalmagnitude > noisethreshold;
    %     noisysecs = reshape(noisyframes,[Fs,ExpDurationSeconds]);
    %     noisysecs = sum(noisysecs,1);
    %     noiseindex = noisysecs>0;
    %     signalmatrix = reshape(signalsdsdt,[Fs,ExpDurationSeconds]);
    %     signalmatrix(:,noiseindex) = [];
    %     signalsdsdt = signalmatrix(:);
    
    %Set the parameters for Chronux analysis
    params.Fs = Fs;
    params.pad = 2;
    params.tapers = [1 1];
    params.fpass = [3 120];
    params.err = 0;
    params.trialave = 0;
    params.segave = 0;
    
    %% Creates a filter in the beta frequency band (20-30Hz) and applies
    %it to the signal. Then uses a hilbert transform to find the envelope
    %amplitude of the signal, and the instantaneous phase.
    betafilter = designfilt('bandpassiir','FilterOrder',2,'HalfPowerFrequency1',20,'HalfPowerFrequency2',30,'DesignMethod','butter','SampleRate',Fs);
    betasignal = filtfilt(betafilter,signalsdsdt);
    betasignalmagnitude = abs(hilbert(betasignal));
    zscorebetasignalmagnitude = zscore(betasignalmagnitude);
    betasignalphase = rad2deg(angle(hilbert(betasignal)));
    
    %Sets the beta burst detection parameters. The burst must cross 3 standard
    %deviations from the mean for at least 150ms, to be classified.
    burststdthreshold = 2;
    minbetaburstdurationms = 150;
    minbetaburstsamples = Fs*(minbetaburstdurationms/1000);
    minburstcycles = 3;
    
    %The amplitude detection parameter is applied to the beta signal, in order
    %to find indexes of when the burst reaches the detection threshold. A second threshold of half the detection threshold
    % is also used to find the starts and stops of each putative burst, in
    % order to include the sides of each burst.
    burstdetectionthreshold = zscorebetasignalmagnitude>=burststdthreshold;
    burstsidethreshold = zscorebetasignalmagnitude>=burststdthreshold/2;
    
    %The two thresholds are combined so that a burst would show as 01210,
    %ie below threshold 1, above threshold 1, above thresholds 1+2, and
    %then back down.
    BdtBst = burstdetectionthreshold + burstsidethreshold;
    difBdtBst = diff(BdtBst);
    %Find is used to find where these changes between threshold levels
    %occur in the signal.
    [blocations,~,bchanges] = find(difBdtBst);
    revblocations = flipud(blocations);
    %strfind is used to find the starts of sequences with the change
    %+1+1-1-1, ie 01210, and also the starts of sequences with the change
    %-1-1+1+1, ie the reverse, to find the stops of these 01210 sequences.
    bpatternseek = strfind(bchanges',[1 1 -1 -1]);
    revbpatternseek = strfind(flipud(bchanges)',[-1 -1 1 1]);
    betaburststarts = blocations(bpatternseek)+1;
    betaburststops = flipud(revblocations(revbpatternseek));
    %These stretches are then removed from the data, to turn it into 01110,
    %so that bursts that don't return to 0, such as 0121210, can be
    %detected without redetecting the previously discovered bursts. These
    %bursts that don't return to 0 are referred to as "waves" for now.
    for idy = 1:length(betaburststarts)
        burstdetectionthreshold(betaburststarts(idy):betaburststops(idy)) = 0;
    end
    %Bursts that started before recording are discarded, as are those that stopped after recording.
    if betaburststops(1)<betaburststarts(1)
        betaburststops = betaburststops(2:end);
    end
    if betaburststarts(end)>betaburststops(end)
        betaburststarts = betaburststarts(1:end-1);
    end
    %In order to detect "waves", we can just look at the detection
    %threshold and only pick out the 2's of 0121's, rather than 121's. This
    %method prevents excessive combination of adjacent waves into extra long
    %events.
    difBdt = diff(burstdetectionthreshold);
    [wlocations,~,wchanges] = find(difBdt);
    revwlocations = flipud(wlocations);
    wpatternseek = strfind(wchanges',[1 -1]);
    revwpatternseek = strfind(flipud(wchanges)',[-1 1]);
    wavestarts = wlocations(wpatternseek);
    wavestops = flipud(revwlocations(revwpatternseek)+1);
    
    %"Waves" that started before recording are discarded, as are those that stopped after recording.
    if length(wavestops)>1
        if wavestops(1)<wavestarts(1)
            wavestops = wavestops(2:end);
        end
        if wavestarts(end)>wavestops(end)
            wavestarts = wavestarts(1:end-1);
        end
        
        wavelength = wavestops-wavestarts;
        wavestarts = wavestarts(wavelength>minbetaburstsamples);
        wavestops = wavestops(wavelength>minbetaburstsamples);
    else
    end
    
    %A previosly assigned MINIMUM DURATION parameter is applied to find bursts with
    %sufficient duration to be classified, while all shorter events are removed.
    betaburstlength = betaburststops-betaburststarts;
    betaburststarts = betaburststarts(betaburstlength>minbetaburstsamples);
    betaburststops = betaburststops(betaburstlength>minbetaburstsamples);
    
    %The burst starts and "wave" starts are combined and put in
    %ascending order, as are the stops. These starts and stops are now the indexes of all putative bursts.
    ratiob2w = length(betaburststarts)/length(wavestarts);
    betaburststarts = [betaburststarts;wavestarts];
    betaburststops = [betaburststops;wavestops];
    [betaburststarts,sortI] = sort(betaburststarts);
    betaburststops = betaburststops(sortI);
    betaburstlength = betaburststops-betaburststarts;
    
    %     burstphase = cell(length(betaburststarts),1);
    %     burstcycles = zeros(length(betaburststarts),1);
    %     for indbcycles = 1:length(betaburststarts)
    %         burstphase{indbcycles} = betasignalphase(betaburststarts(indbcycles):betaburststops(indbcycles),:);
    %         minBValue1 = min(burstphase{indbcycles});
    %         paddedBSignal1 = [minBValue1;burstphase{indbcycles};minBValue1];
    %         %As the phase signal is constantly increasing, findpeaks picks up
    %         %where each cycle ends. As there is no built in "findtroughs"
    %         %function, in order to find the start of each cycle, this is all
    %         %repeated for the -ve of the phase signal.
    %         [bcyclestart,bcyclestartloc] = findpeaks(paddedBSignal1);
    %         minBValue2 = min(-burstphase{indbcycles});
    %         paddedBSignal2 = [minBValue2;-burstphase{indbcycles};minBValue2];
    %         [bcyclestop,bcyclestoploc] = findpeaks(paddedBSignal2);
    %         bcyclestop = -bcyclestop;
    %         if bcyclestoploc(1)<bcyclestartloc(1)
    %             bcyclestop = bcyclestop(2:end);
    %         end
    %         if bcyclestartloc(end)>bcyclestoploc(end)
    %             bcyclestart = bcyclestart(1:end-1);
    %         end
    %         %The length of each cycle is calculated, and then the total length
    %         %of all these cycles is calculated, and then divided by 360degrees
    %         %(a full cycle phase), and rounded down to get the number of
    %         %complete burst cycles.
    %         bcyclelength = bcyclestart-bcyclestop;
    %         totalbcycles = sum(bcyclelength);
    %         burstcycles(indbcycles) = floor(totalbcycles/360);
    %     end
    
    %Artefact removal on a burst by burst basis, by detecting large almost
    %instantaneous changes in the unfiltered LFP.
    for ib = 1:length(betaburststarts)
        pbm(ib) = max(betasignalmagnitude(betaburststarts(ib):betaburststops(ib),:));
    end
    discard = isoutlier(pbm);
    params.discardrate = sum(discard)/length(betaburststarts);
    betaburststarts(discard) = [];
    betaburststops(discard) = [];
    betaburstlength(discard) = [];
    
    %Now the burst start and stop indexes of all putative bursts that do
    %not meet these 2 criteria have been removed, all remaining can
    %be counted as bursts, so we can calculate the total number of bursts
    %in this session, as well as the duration of each burst in
    %milliseconds.
    betaburstnumber = length(betaburstlength);
    burstduration = (betaburstlength/Fs)*1000;
    meanburstduration = mean(burstduration);
    stdburstduration = std(burstduration);
    semburstduration = stdburstduration/sqrt(length(burstduration));
    betaburststartsds = ceil(betaburststarts./Fs);
    betaburststopsds = ceil(betaburststops./Fs);
    
    %Finally the burst start and stop indexes are used to isolate the
    %filtered and unfiltered bursts, and the peak magnitude of each burst
    %is ascertained. These are averaged, and SEM calculated to gain an
    %idea of the burst dynamics for the session.
    betabursts = cell(betaburstnumber,1);
    betaburstsunfiltered = cell(betaburstnumber,1);
    peakburstmagnitude = zeros(betaburstnumber,1);
    if betaburstnumber>0
        for ind6 = 1:betaburstnumber
            betabursts{ind6} = betasignal(betaburststarts(ind6):betaburststops(ind6),:);
            betaburstsunfiltered{ind6} = signalsdsdt(betaburststarts(ind6):betaburststops(ind6),:);
            peakburstmagnitude(ind6) = max(betasignalmagnitude(betaburststarts(ind6):betaburststops(ind6),:));
            burstphase{ind6} = betasignalphase(betaburststarts(ind6):betaburststops(ind6),:);
            burstspeed(ind6) = mean(speed(:,betaburststartsds(ind6):betaburststopsds(ind6)));
            %         burstx(ind6) = nanmean(x(betaburststartscds(ind6):betaburststopscds(ind6)));
            %         bursty(ind6) = nanmean(y(betaburststartscds(ind6):betaburststopscds(ind6)));
            
            burstmiddle = length(betabursts{ind6})/2;
            [~,locs] = findpeaks(betabursts{ind6},'MinPeakDistance',100);
            [~,closestpeak] = min(abs(locs-burstmiddle));
            middlepeak(ind6) = locs(closestpeak);
        end
        meanburstmagnitude = mean(peakburstmagnitude);
        stdburstmagnitude = std(peakburstmagnitude);
        semburstmagnitude = stdburstmagnitude/sqrt(length(peakburstmagnitude));
        biggestburstleftlength = middlepeak(find(betaburstlength == max(betaburstlength)));
        biggestburstrightlength = max(betaburstlength)-biggestburstleftlength;
        
        %As in (Shin et al 2017), to attempt to determine the mechanism
        %underlying beta bursts (bursty vs dynamic amplitude modulation),
        %calculate the time lag modulus for the burst and non-burst segments.
        for indb = 1:betaburstnumber
            [~,locs] = findpeaks(betabursts{indb},'MinPeakDistance',100);
            tla(indb) = locs(end)-locs(1);
            %Also calculate average burst frequency
            burstfreq(indb) = Fs/((locs(end)-locs(1))/(length(locs)-1));
        end
        tlmb = (tla - 120.*round(tla./120))/120;
        tlmbh = histcounts(tlmb,-0.5:0.1:0.5,'Normalization','probability');
        for indb = 1:betaburstnumber
            if indb == 1
                nonburstsignal = betasignal(1:betaburststarts(indb),:);
            else
                nonburstsignal = betasignal(betaburststops(indb-1):betaburststarts(indb),:);
            end
            if length(nonburstsignal)>=200
                [~,locs] = findpeaks(nonburstsignal,'MinPeakDistance',100);
                tlb(indb) = locs(end)-locs(1);
                %Also calculate average burst frequency
                nonburstfreq(indb) = Fs/((locs(end)-locs(1))/(length(locs)-1));
            else
                tlb(indb) = NaN;
            end
        end
        tlmn = (tlb - 120.*round(tlb./120))/120;
        tlmnh = histcounts(tlmn,-0.5:0.1:0.5,'Normalization','probability');
        
        meanburstfreq = mean(burstfreq);
        stdburstfreq = std(burstfreq);
        semburstfreq = stdburstfreq/sqrt(length(burstfreq));
        meannonburstfreq = mean(nonburstfreq);
        stdnonburstfreq = std(nonburstfreq);
        semnonburstfreq = stdnonburstfreq/sqrt(length(nonburstfreq));
        
        %As bursts vary in length, in order to compare them time locked to
        %eachother, they are padded with varying numbers of NaNs on either side
        %to center them all.
        burstsegpad = cell(1,betaburstnumber);
        burstrawsegpad = cell(1,betaburstnumber);
        for ind7 = 1:betaburstnumber
            burstprepadsize = floor(biggestburstleftlength - middlepeak(ind7));
            burstpostpadsize = round(biggestburstrightlength - (betaburstlength(ind7)-middlepeak(ind7)));
            if any([burstprepadsize burstpostpadsize] < 0)
                burstsegpad{ind7}(max(betaburstlength)+1,1) = NaN;
                burstrawsegpad{ind7}(max(betaburstlength)+1,1) = NaN;
            else
                burstsegpad{ind7} = padarray(betabursts{ind7},burstprepadsize,NaN,'pre');
                burstsegpad{ind7} = padarray(burstsegpad{ind7},burstpostpadsize,NaN,'post');
                burstrawsegpad{ind7} = padarray(betaburstsunfiltered{ind7},burstprepadsize,NaN,'pre');
                burstrawsegpad{ind7} = padarray(burstrawsegpad{ind7},burstpostpadsize,NaN,'post');
            end
        end
        burstsegpad = cell2mat(burstsegpad);
        burstrawsegpad = cell2mat(burstrawsegpad);
        meanburstseg = nanmean(burstsegpad,2)';
        burstvar = nanstd(burstsegpad,1,2)';
        burstvarlow = meanburstseg-burstvar;
        burstvarhigh = meanburstseg+burstvar;
        [~,meanburstpeaks] = findpeaks(meanburstseg,'MinPeakDistance',100);
        meanburstpeaks = meanburstpeaks./(Fs/1000);
        
        %In order to investigate oscillatory dynamics before, during and after
        %bursts, power spectra of each burst segment, as well as equal length
        %segments before and after the burst are generated and averaged across
        %bursts. Any bursts that less than a bursts length from the start or
        %end of the session are discarded from this particular analysis.
        fFinterp = linspace(3,120,234);
        bbstarts = betaburststarts(betaburststarts > betaburstlength);
        bblength = betaburstlength(betaburststarts > betaburstlength);
        bbstops = betaburststops(betaburststarts > betaburstlength);
        %The start of a burst is called its start, the start of the whole
        %segment, including the pre-burst, is called the beginning. The same
        %applies to its stop, and its end, respectively.
        bbbegin = bbstarts-bblength;
        bbstarts(bbbegin<=0) = [];
        bblength(bbbegin<=0) = [];
        bbstops(bbbegin<=0) = [];
        bbend = bbstops+bblength;
        bbstarts(bbend>length(signalsdsdt)) = [];
        bblength(bbend>length(signalsdsdt)) = [];
        bbstops(bbend>length(signalsdsdt)) = [];
        newbetaburstnumber = length(bbstarts);
        
        %The signal segments are indexed and then undergo spectral analysis.
        %Each burst spectrum has 50Hz and 100Hz removed and linearly
        %interpolated, then it is logged.
        signalbeforebursts = cell(newbetaburstnumber,1);
        signalduringbursts = cell(newbetaburstnumber,1);
        signalafterbursts = cell(newbetaburstnumber,1);
        spectrumbeforebursts = cell(newbetaburstnumber,1);
        spectrumduringbursts = cell(newbetaburstnumber,1);
        spectrumafterbursts = cell(newbetaburstnumber,1);
        for ind8 = 1:newbetaburstnumber
            signalbeforebursts{ind8} = signalsdsdt(bbstarts(ind8)-bblength(ind8):bbstarts(ind8),:);
            signalduringbursts{ind8} = signalsdsdt(bbstarts(ind8):bbstops(ind8),:);
            signalafterbursts{ind8} = signalsdsdt(bbstops(ind8):bbstops(ind8)+bblength(ind8),:);
            fullsignal{ind8} = signalsdsdt(bbstarts(ind8)-bblength(ind8):bbstops(ind8)+bblength(ind8));
            filtsignal{ind8} = betasignal(bbstarts(ind8)-bblength(ind8):bbstops(ind8)+bblength(ind8));
            [spectrumbeforebursts{ind8},fF{ind8}] = mtspectrumc(signalbeforebursts{ind8},params);
            [spectrumduringbursts{ind8}] = mtspectrumc(signalduringbursts{ind8},params);
            [spectrumafterbursts{ind8}] = mtspectrumc(signalafterbursts{ind8},params);
            spectrumbeforeburstsi(:,ind8) = interp1(fF{ind8},spectrumbeforebursts{ind8},fFinterp);
            spectrumbeforeburstsi(:,ind8) = 10*log10(spectrumbeforeburstsi(:,ind8));
            %         spectrumbeforeburstsi(fFinterp>40 & fFinterp<60,ind8) = NaN;
            %         spectrumbeforeburstsi(fFinterp>90 & fFinterp<110,ind8) = NaN;
            %         spectrumbeforeburstsi(:,ind8) = fillmissing(spectrumbeforeburstsi(:,ind8),'linear');
            spectrumduringburstsi(:,ind8) = interp1(fF{ind8},spectrumduringbursts{ind8},fFinterp);
            spectrumduringburstsi(:,ind8) = 10*log10(spectrumduringburstsi(:,ind8));
            %         spectrumduringburstsi(fFinterp>40 & fFinterp<60,ind8) = NaN;
            %         spectrumduringburstsi(fFinterp>90 & fFinterp<110,ind8) = NaN;
            %         spectrumduringburstsi(:,ind8) = fillmissing(spectrumduringburstsi(:,ind8),'linear');
            spectrumafterburstsi(:,ind8) = interp1(fF{ind8},spectrumafterbursts{ind8},fFinterp);
            spectrumafterburstsi(:,ind8) = 10*log10(spectrumafterburstsi(:,ind8));
            %         spectrumafterburstsi(fFinterp>40 & fFinterp<60,ind8) = NaN;
            %         spectrumafterburstsi(fFinterp>90 & fFinterp<110,ind8) = NaN;
            %         spectrumafterburstsi(:,ind8) = fillmissing(spectrumafterburstsi(:,ind8),'linear');
        end
        %These burst spectra are then averaged across all bursts.
        meansbb = mean(spectrumbeforeburstsi,2);
        meansdb = mean(spectrumduringburstsi,2);
        meansab = mean(spectrumafterburstsi,2);
        
    else
    end
    
    %% Stores all data of interest as a non-scalar structure, making it easy
    %to look at and compare data from different sessions.
    %First the parameters and settings.
    EventData(ind0).Session = session;
    EventData(ind0).Parameters = params;
    EventData(ind0).Parameters.MinimumBurstDurationms = minbetaburstdurationms;
    EventData(ind0).Parameters.MinimumBurstCycles = minburstcycles;
    EventData(ind0).Parameters.BurstThreshold = burststdthreshold;
    
    %Then Burst Data
    EventData(ind0).NumberOfBetaBursts = betaburstnumber;
    if betaburstnumber>0
        EventData(ind0).BurstStarts = betaburststarts;
        EventData(ind0).BurstStops = betaburststops;
        EventData(ind0).BurstMagnitude = peakburstmagnitude;
        EventData(ind0).BurstLength = betaburstlength;
        EventData(ind0).BurstDuration = burstduration;
        EventData(ind0).BurstFreq = burstfreq;
        EventData(ind0).BurstTimeLagModulusHist = tlmbh;
        EventData(ind0).NonBurstTimeLagModulusHist = tlmnh;
        EventData(ind0).AllBurstSegments = burstsegpad;
        EventData(ind0).AllBurstRawSegments = burstrawsegpad;
        EventData(ind0).BurstRS = burstspeed;
        EventData(ind0).FullBurstSignal = fullsignal;
        EventData(ind0).FiltBurstSignal = filtsignal;
        EventData(ind0).SignalBeforeBursts = signalbeforebursts;
        EventData(ind0).SignalDuringBursts = signalduringbursts;
        EventData(ind0).SignalAfterBursts = signalafterbursts;
        EventData(ind0).SpectrumBeforeBursts = meansbb;
        EventData(ind0).SpectrumDuringBursts = meansdb;
        EventData(ind0).SpectrumAfterBursts = meansab;
        EventData(ind0).Signal = signalsdsdt;
    else
    end
    
    clearvars -except Folder Channel folder foldername subfolders sessionnames ind0 animal channelname EventData
    
end



RipplesAndBursts = strcat(foldername,'\','RipplesAndBursts');
mkdir(RipplesAndBursts);
cd(RipplesAndBursts);

NumberOfBetaBursts = [EventData(:).NumberOfBetaBursts];
figure2name = strcat(animal,'-',channelname,'-','NumberOfBetaBursts');
figure('Name',figure2name);
subplot(2,2,[1,3])
bar(NumberOfBetaBursts,'b');
set(gca,'LineWidth',1.5);
set(gca,'FontSize',12);
xlabel('Session','fontsize',12);
ylabel('Number Of Beta Bursts','fontsize',14);
set(gca,'box','off');
BurstDuration = vertcat(EventData(:).BurstDuration);
subplot(2,2,2);
histogram(BurstDuration,'FaceColor','b');
set(gca,'LineWidth',1.5);
set(gca,'FontSize',12);
set(gca,'TickLength',[0.02,0.02]);
xlabel('Burst Length (ms)','fontsize',12);
set(gca,'box','off');
BurstMagnitude = vertcat(EventData(:).BurstMagnitude);
subplot(2,2,4);
histogram(BurstMagnitude,'FaceColor','b');
set(gca,'LineWidth',1.5);
set(gca,'FontSize',12);
set(gca,'TickLength',[0.02,0.02]);
xlabel('Burst Magnitude (µV)','fontsize',12);
set(gca,'box','off');
set(gcf, 'Position',  [250, 100, 750, 500])
hgsave(figure2name)

for indz2 = 1:10
    Fs = EventData(indz2).Parameters.Fs;
    BetaBurstStarts{indz2} = (EventData(indz2).BurstStarts/(Fs*60));
    BurstDist{indz2} = cumsum(histcounts(BetaBurstStarts{indz2},0:(1/60):15)');
    EventData(indz2).BurstDistribution = BurstDist{indz2};
end

NovelBurstDist = BurstDist{9};
FamiliarBurstDist = BurstDist{8};
NormNovelBurstDist = NovelBurstDist/max(NovelBurstDist)*100;
NormFamiliarBurstDist = FamiliarBurstDist/max(FamiliarBurstDist)*100;

figure6name = strcat(animal,'-',channelname,'-','BurstProfileNovelty');
figure('Name',figure6name);
subaxis3 = subplot(2,1,1);
plot(NovelBurstDist,'LineWidth',2,'Color','r')
hold on
plot(FamiliarBurstDist,'LineWidth',2,'Color','b')
xlim([0 900])
set(gca,'xTick',0:60:900)
set(gca,'Xticklabel',[])
ylab1 = ylabel('Total Number of Bursts','fontsize',14);
set(ylab1, 'Units', 'Normalized', 'Position', [-0.08, 0.5, 0]);
set(gca,'FontSize',12);
legend({'Novel','Familiar'},'Location','northwest');
legend boxoff
set(gca,'box','off');
subaxis4 = subplot(2,1,2);
plot(NormNovelBurstDist,'LineWidth',2,'Color','r')
hold on
plot(NormFamiliarBurstDist,'LineWidth',2,'Color','b')
xlim([0 900])
set(gca,'xTick',0:60:900)
set(gca,'xTickLabel',0:1:15)
alpha(gca,0.75)
xlabel('Time (min)','fontsize',14);
ylab2 = ylabel('Percentage of Bursts (%)','fontsize',14);
set(ylab2, 'Units', 'Normalized', 'Position', [-0.08, 0.5, 0]);
set(gca,'FontSize',12);
set(gca,'box','off');
set(gcf, 'Position',  [300, 100, 500, 630])
set(subaxis3,'Position',[0.13 0.55 0.775 0.4]);
set(subaxis4,'Position',[0.13 0.11 0.775 0.4]);
hgsave(figure6name)

%%
SpectrumBeforeBursts = [EventData(:).SpectrumBeforeBursts]';
fFinterp = linspace(3,120,234);
MMSBB = mean(SpectrumBeforeBursts);
stdSBB = std(SpectrumBeforeBursts);
semSBB = stdSBB/sqrt(size(SpectrumBeforeBursts,1));
SBBLow = MMSBB-semSBB;
SBBHigh = MMSBB+semSBB;
SpectrumDuringBursts = [EventData(:).SpectrumDuringBursts]';
MMSDB = mean(SpectrumDuringBursts);
stdSDB = std(SpectrumDuringBursts);
semSDB = stdSDB/sqrt(size(SpectrumDuringBursts,1));
SDBLow = MMSDB-semSDB;
SDBHigh = MMSDB+semSDB;
SpectrumAfterBursts = [EventData(:).SpectrumAfterBursts]';
MMSAB = mean(SpectrumAfterBursts);
stdSAB = std(SpectrumAfterBursts);
semSAB = stdSAB/sqrt(size(SpectrumAfterBursts,1));
SABLow = MMSAB-semSAB;
SABHigh = MMSAB+semSAB;

figure7name = strcat(animal,'-',channelname,'-','BeforeDuringAfterBurstSpectrums');
figure('Name',figure7name);
SBBbar = patch([fFinterp fFinterp(end:-1:1) fFinterp(1)],[SBBLow SBBHigh(end:-1:1) SBBLow(1)], 'r');
hold on;
SBBline = line(fFinterp,MMSBB,'LineWidth',1.5);
SDBbar = patch([fFinterp fFinterp(end:-1:1) fFinterp(1)],[SDBLow SDBHigh(end:-1:1) SDBLow(1)], 'r');
SDBline = line(fFinterp,MMSDB,'LineWidth',1.5);
SABbar = patch([fFinterp fFinterp(end:-1:1) fFinterp(1)],[SABLow SABHigh(end:-1:1) SABLow(1)], 'r');
SABline = line(fFinterp,MMSAB,'LineWidth',1.5);
alpha(0.3);
set(SBBbar, 'facecolor', [0 1 0.8667], 'edgecolor', 'none');
set(SDBbar, 'facecolor', [0 (12/255) 1], 'edgecolor', 'none');
set(SABbar, 'facecolor', [1 0 (233/255)], 'edgecolor', 'none');
set(SBBline, 'color', [0 0.8 0.8]);
set(SDBline, 'color', [0 0 0.8]);
set(SABline, 'color', [0.8 0 0.8]);
set(gca,'LineWidth',1.5);
set(gca,'FontSize',12)
xlabel('Frequency (Hz)','fontsize',14);
ylabel('Power (DB)','fontsize',14);
legend([SBBline SDBline SABline],{'Before Burst','During Burst','After Burst'});
legend boxoff
set(gca,'box','off');
set(gcf, 'Position',  [400, 100, 600, 500])
% xlim([3 40])
hgsave(figure7name)

%Save the finished data structure as a .mat file for later use
EventFolder = strcat(foldername(1:end-2),'Analysed\','EventData');
mkdir(EventFolder);
cd(EventFolder);
save(strcat('EventData','-',(animal),'-',(channelname),'.mat'),'EventData');

end

