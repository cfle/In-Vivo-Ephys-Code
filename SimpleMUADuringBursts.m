%Speed-binned LFP analysis using the SpeedSnip method. Suitable for one
%channel at a time.

function [BurstMUAData] = SimpleMUADuringBursts(Folder)

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
    
    cd(Folder)
    files = dir(subfolders{ind0});
    files = files(arrayfun(@(x) ~strcmp(x.name(1),'.'),files));
    channels = {files.name};
    channels = channels(contains(channels,'CH'));
    channels = strcat(subfolders{ind0},'\',channels);

    animalnumber = Folder(end-1:end);
    animal = strcat('Mouse',animalnumber);
    session = files.folder;
    session = session(end-4:end);
    
    [~,~,extension] = fileparts(channels{1});
    signals = [];
    
    for ich = 1:length(channels)
        if strcmp(extension,'.continuous')
            %Imports the .continuous data into MatLab
            [signals(:,ich),~,info] = load_open_ephys_data_faster(channels{ich});
            channelnames{ich} = channels{ich}(end-14:end-11);
        elseif strcmp(extension,'.mat')
            %Imports the .continuous data into MatLab
            load(channels{ich});
            signals(:,ich) = sig;
            info.header.sampleRate = 30000;
            channelnames{ich} = channels{ich}(end-14:end-11);
        end
    end
    
    %Cuts the signals to the experiment duration ie 15minutes in this case
    ExpDurationMinutes = 15;
    Fs = info.header.sampleRate;
    ExpDurationFrames = (ExpDurationMinutes*60)*Fs;
    ExpDurationSeconds = ExpDurationMinutes*60;
    signals = signals(1:ExpDurationFrames,:);
    
    %Downsamples the signals and timestamps, to reduce computation time,
    %and detrends to remove any baseline fluctuations
    signalsdt = detrend(signals);
    
    %Set the parameters for Chronux analysis
    params.fpass = [500 5000];
    
    for ich = 1:length(channels)
        
        %Common Average Referencing
        refchannels = 1:length(channelnames);
        refchannels(ich) = [];
        refsignal = mean(signalsdt(:,refchannels),2);
        signaldt = signalsdt(:,ich);
        signal = signaldt-refsignal;
        
        % Creates a filter in the MUA frequency band (600-6000Hz) and applies
        %it to the signal. Then uses a hilbert transform to find the envelope
        %amplitude of the signal, and the instantaneous phase.
        muafilter = designfilt('bandpassiir','FilterOrder',4,'HalfPowerFrequency1',params.fpass(1),'HalfPowerFrequency2',params.fpass(2),'DesignMethod','butter','SampleRate',Fs);
        muasignal = filtfilt(muafilter,signal);
        
        threshold = 4*median(abs(muasignal)/0.6745);
        
        %Sets the beta burst detection parameters. The burst must cross 3 standard
        %deviations from the mean for at least 150ms, to be classified.
        [~,peaklocations] = findpeaks(abs(muasignal),'MinPeakHeight',threshold,'MinPeakDistance',Fs/2000);
        
        peakstarts = peaklocations - Fs/1000;
        peakstops = peaklocations + Fs/1000;
        
        peaksbeforesession = peakstarts <= Fs;
        peakstarts(peaksbeforesession) = [];
        peakstops(peaksbeforesession) = [];
        peaksaftersession = peakstops >= (ExpDurationFrames-Fs);
        peakstarts(peaksaftersession) = [];
        peakstops(peaksaftersession) = [];
        
        peaklength = peakstops-peakstarts;
        peaknumber = length(peaklength);
                
%         for indpeaks = 1:peaknumber
%             waveforms(:,indpeaks) = muasignal(peakstarts(indpeaks):peakstops(indpeaks));
%             [~,maxind] = max(abs(waveforms(:,indpeaks)));
%             wfamplitude(indpeaks) = waveforms(maxind,indpeaks);
%         end
        
        for indpeaks = 1:peaknumber
            waveforms(:,indpeaks) = muasignal(peakstarts(indpeaks):peakstops(indpeaks));
            [~,maxind(indpeaks)] = max(abs(waveforms(:,indpeaks)));
            peakstarts(indpeaks) = peakstarts(indpeaks) + (maxind(indpeaks)-31);
            peakstops(indpeaks) = peakstops(indpeaks) + (maxind(indpeaks)-31);
            if waveforms(maxind(indpeaks),indpeaks)<0
                waveforms(:,indpeaks) = muasignal(peakstarts(indpeaks):peakstops(indpeaks));
            else
                waveforms(:,indpeaks) = -muasignal(peakstarts(indpeaks):peakstops(indpeaks));
            end
            wfamplitude(indpeaks) = waveforms(31,indpeaks);
        end
        peaks = [peakstarts peakstops];
        [~,idx,~] = unique(peaks,'rows');
        peakstarts = peakstarts(idx);
        peakstops = peakstops(idx);
        waveforms = waveforms(:,idx);
        wfamplitude = wfamplitude(idx);
        
        % Creates a filter in the beta frequency band (20-30Hz) and applies
        %it to the signal. Then uses a hilbert transform to find the envelope
        %amplitude of the signal, and the instantaneous phase.
        betafilter = designfilt('bandpassiir','FilterOrder',2,'HalfPowerFrequency1',20,'HalfPowerFrequency2',30,'DesignMethod','butter','SampleRate',Fs);
        betasignal = filtfilt(betafilter,signaldt);
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
        
        %Artefact removal on a burst by burst basis, by detecting large almost
        %instantaneous changes in the unfiltered LFP.
        for ib = 1:length(betaburststarts)
            pbm(ib) = max(betasignalmagnitude(betaburststarts(ib):betaburststops(ib),:));
        end
        discard = isoutlier(pbm);
        betaburststarts(discard) = [];
        betaburststops(discard) = [];
        betaburstlength(discard) = [];
        
        %% Now the burst start and stop indexes of all putative bursts that do
        %not meet these 2 criteria have been removed, all remaining can
        %be counted as bursts, so we can calculate the total number of bursts
        %in this session, as well as the duration of each burst in
        %milliseconds.
        SegLengthF = Fs;
        betaburststarts(betaburststarts < SegLengthF) = [];
        betaburstlength(betaburststarts < SegLengthF) = [];
        betaburststarts(betaburststarts + Fs > ExpDurationFrames) = [];
        betaburstlength(betaburststarts + Fs > ExpDurationFrames) = [];
        
        % Detect beta burst segments (burststart-Fs:burststart+Fs) that
        % overlap with eachother and remove them
        bbstarts = betaburststarts-Fs;
        bbstops = betaburststops+Fs;
        for i = 1:length(bbstarts)
            burstseg{i} = bbstarts(i):bbstops(i);
        end
        for i = 1:length(bbstarts)
            for i2 = 1:length(bbstarts)
                overlap(i,i2) = ~isempty(intersect(burstseg{i},burstseg{i2}));
            end
        end
        OL = sum(overlap,1)>1;
        
%         OverlappingBursts = diff(betaburststarts)<=SegLengthF;
%         OL = ([0;OverlappingBursts]+[OverlappingBursts;0]>0);
        betaburststarts(logical(OL)) = [];
        betaburststops(logical(OL)) = [];
        betaburstlength(logical(OL)) = [];
        numberofbetabursts = length(betaburststarts);
        
        if numberofbetabursts > 0
             
            rasteredges = [-SegLengthF SegLengthF];
            spikeraster = zeros(ExpDurationFrames,1);
            spikeraster(peakstarts) = 1;
            
            for indbr = 1:numberofbetabursts
                periburstspikeR(indbr,:) = spikeraster(betaburststarts(indbr)-Fs:betaburststarts(indbr)+Fs);
                burstspikes(indbr,:) = peakstarts>=betaburststarts(indbr)&peakstarts<=betaburststops(indbr);
            end
            periburstspikeRall = sum(periburstspikeR,1);
            burstspikes = logical(sum(burstspikes,1));
            nonburstspikes = ~burstspikes;
            meanburstspike = mean(waveforms(:,burstspikes),2);
            meannonburstspike = mean(waveforms(:,nonburstspikes),2);
            
            histogramedges = 1:1500:60001;
            for ibin = 1:length(histogramedges)-1
                periburstspikeH(ibin) = sum(periburstspikeRall(histogramedges(ibin):histogramedges(ibin+1)));
            end
            periburstspikeH = periburstspikeH./(sum(periburstspikeH));
            
        else
        end
        
        %% Stores all data of interest as a non-scalar structure, making it easy
        %to look at and compare data from different sessions.
        %First the parameters and settings.
        BurstMUAData(ich).Channel = channelnames{ich};
        %Then Burst Data
        BurstMUAData(ich).NumberOfBetaBursts = numberofbetabursts;
        BurstMUAData(ich).TimeStamps = peakstarts;
        BurstMUAData(ich).Waveforms = waveforms;
        if numberofbetabursts>0
            BurstMUAData(ich).BurstLengths = betaburstlength;
            BurstMUAData(ich).RasterEdges = rasteredges;
            BurstMUAData(ich).PeriBurstSpikeRaster = periburstspikeR;
            BurstMUAData(ich).HistogramEdges = histogramedges;
            BurstMUAData(ich).PeriBurstSpikeHistogram = periburstspikeH;
            BurstMUAData(ich).MeanBurstSpike = meanburstspike;
            BurstMUAData(ich).MeanNonBurstSpike = meannonburstspike;
        else
        end
        
    clearvars -except Folder folder foldername subfolders sessionnames ind0 ich animal channelnames signalsdt Fs params ExpDurationFrames BurstMUAData    
        
    end

    %Save the finished data structure as a .mat file for later use
    MUAFolder = strcat(foldername(1:end-2),'Analysed\','MUAData');
    warning('off', 'MATLAB:MKDIR:DirectoryExists');
    mkdir(MUAFolder);
    cd(MUAFolder);
    save(strcat('BurstMUAData','-',(animal),'-',sessionnames{ind0},'.mat'),'BurstMUAData');
    
end

end

