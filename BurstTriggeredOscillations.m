%Speed-binned LFP analysis using the SpeedSnip method. Suitable for one
%channel at a time.

function [BTOData] = BurstTriggeredOscillations(Folder,DetectorChannel,ComparatorChannel)

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
    channel1 = channels(contains(channels,DetectorChannel));
    channel1 = strcat(subfolders{ind0},'\',channel1);
    channel1name = string(channel1{1}(end-14:end-11));
    channel2 = channels(contains(channels,ComparatorChannel));
    channel2 = strcat(subfolders{ind0},'\',channel2);
    channel2name = string(channel2{1}(end-14:end-11));
    channels = [channel1 channel2];
    
    animalnumber = Folder(end-1:end);
    animal = strcat('Mouse',animalnumber);
    session = files.folder;
    session = session(end-4:end);
    
    for indchannel = 1:2
        [~,~,extension] = fileparts(channel1{1});
        if strcmp(extension,'.continuous')
            %Imports the .continuous data into MatLab
            [signal(:,indchannel),~,info] = load_open_ephys_data_faster(channels{indchannel});
        elseif strcmp(extension,'.mat')
            %Imports the .continuous data into MatLab
            load(channels{indchannel});
            signal(:,indchannel) = sig;
            info.header.sampleRate = 30000;
        end
    end
    
    signal1 = signal(:,1);
    signal2 = signal(:,2);
    
    %Cuts the signals to the experiment duration ie 15minutes in this case
    ExpDurationMinutes = 15;
    Fs = info.header.sampleRate;
    ExpDurationFrames = (ExpDurationMinutes*60)*Fs;
    ExpDurationSeconds = ExpDurationMinutes*60;
    signal1 = signal1(1:ExpDurationFrames,:);
    signal2 = signal2(1:ExpDurationFrames,:);
    
    %Downsamples the signals and timestamps, to reduce computation time,
    %and detrends to remove any baseline fluctuations
    downsampleX = 10;
    Fs = Fs/downsampleX;
    ExpDurationFrames = ExpDurationFrames/downsampleX;
    signal1ds = downsample(signal1,downsampleX);
    signal1dsdt = detrend(signal1ds);
    signal2ds = downsample(signal2,downsampleX);
    signal2dsdt = detrend(signal2ds);
    
    path = regexp(Folder,filesep,'split');
    path = path(1:end-1);
    StudyFolder = strjoin(path,'\');
    BurstDataFile = strcat('EventData-',animal,'-',DetectorChannel,'.mat');
    BurstDataPath = strcat(StudyFolder,'\Analysed\EventData\',BurstDataFile);
    load(BurstDataPath);
    BurstStarts = EventData(ind0).BurstStarts;
    BurstStops = EventData(ind0).BurstStops;
    betaburstnumber = length(BurstStarts);
    
    SegLengthS = 1;
    SegLengthmS = SegLengthS*1000;
    SegLengthF = SegLengthS*Fs;
    t = -SegLengthmS:SegLengthmS/SegLengthF:SegLengthmS;
    
    BurstStarts = EventData(ind0).BurstStarts;
    BurstLengths = EventData(ind0).BurstLength;
    NumberOfBetaBursts = EventData(ind0).NumberOfBetaBursts;
    
    BurstSegStarts = BurstStarts-SegLengthF;
    BurstStarts(BurstSegStarts<0) = [];
    BurstLengths(BurstSegStarts<0) = [];
    BurstSegStops = BurstStarts+SegLengthF;
    BurstStarts(BurstSegStops>length(signal1dsdt)) = [];
    BurstLengths(BurstSegStops>length(signal1dsdt)) = [];
    
    OverlappingBursts = diff(BurstStarts)<=SegLengthF;
    OLs = logical([0;OverlappingBursts]);
    OL = logical([OverlappingBursts;0]);
    BurstStarts(OL) = NaN;
    BurstStarts(OLs) = NaN;
    BurstStarts(isnan(BurstStarts)) = [];
    NumberOfBetaBursts = length(BurstStarts);
    
    deltafilter = designfilt('bandpassiir','FilterOrder',2,'HalfPowerFrequency1',1,'HalfPowerFrequency2',5,'DesignMethod','butter','SampleRate',Fs);
    deltasignal1 = filtfilt(deltafilter,signal1dsdt);
    deltasignal2 = filtfilt(deltafilter,signal2dsdt);
    
    thetafilter = designfilt('bandpassiir','FilterOrder',2,'HalfPowerFrequency1',5,'HalfPowerFrequency2',12,'DesignMethod','butter','SampleRate',Fs);
    thetasignal1 = filtfilt(thetafilter,signal1dsdt);
    thetasignal2 = filtfilt(thetafilter,signal2dsdt);
    
    alphafilter = designfilt('bandpassiir','FilterOrder',2,'HalfPowerFrequency1',12,'HalfPowerFrequency2',20,'DesignMethod','butter','SampleRate',Fs);
    alphasignal1 = filtfilt(alphafilter,signal1dsdt);
    alphasignal2 = filtfilt(alphafilter,signal2dsdt);
    
    betafilter = designfilt('bandpassiir','FilterOrder',2,'HalfPowerFrequency1',20,'HalfPowerFrequency2',30,'DesignMethod','butter','SampleRate',Fs);
    betasignal1 = filtfilt(betafilter,signal1dsdt);
    betasignal2 = filtfilt(betafilter,signal2dsdt);
    
    if betaburstnumber>0
        for ib = 1:NumberOfBetaBursts
            %Delta
            BurstDelta(:,ib) = deltasignal1(BurstStarts(ib)-SegLengthF:BurstStarts(ib)+SegLengthF);
            TriggeredDelta(:,ib) = deltasignal2(BurstStarts(ib)-SegLengthF:BurstStarts(ib)+SegLengthF);
            %Theta
            BurstTheta(:,ib) = thetasignal1(BurstStarts(ib)-SegLengthF:BurstStarts(ib)+SegLengthF);
            TriggeredTheta(:,ib) = thetasignal2(BurstStarts(ib)-SegLengthF:BurstStarts(ib)+SegLengthF);
            %Alpha
            BurstAlpha(:,ib) = alphasignal1(BurstStarts(ib)-SegLengthF:BurstStarts(ib)+SegLengthF);
            TriggeredAlpha(:,ib) = alphasignal2(BurstStarts(ib)-SegLengthF:BurstStarts(ib)+SegLengthF);
            %Beta
            BurstBeta(:,ib) = betasignal1(BurstStarts(ib)-SegLengthF:BurstStarts(ib)+SegLengthF);
            TriggeredBeta(:,ib) = betasignal2(BurstStarts(ib)-SegLengthF:BurstStarts(ib)+SegLengthF);
        end
        
        meanBurstDelta = mean(BurstDelta,2);
        semBurstDelta = std(BurstDelta,[],2)/sqrt(length(BurstStarts));
        meanTriggeredDelta = mean(TriggeredDelta,2);
        semTriggeredDelta = std(TriggeredDelta,[],2)/sqrt(length(BurstStarts));
        
        meanBurstTheta = mean(BurstTheta,2);
        semBurstTheta = std(BurstTheta,[],2)/sqrt(length(BurstStarts));
        meanTriggeredTheta = mean(TriggeredTheta,2);
        semTriggeredTheta = std(TriggeredTheta,[],2)/sqrt(length(BurstStarts));
        
        meanBurstAlpha = mean(BurstAlpha,2);
        semBurstAlpha = std(BurstAlpha,[],2)/sqrt(length(BurstStarts));
        meanTriggeredAlpha = mean(TriggeredAlpha,2);
        semTriggeredAlpha = std(TriggeredAlpha,[],2)/sqrt(length(BurstStarts));
        
        meanBurstBeta = mean(BurstBeta,2);
        semBurstBeta = std(BurstBeta,[],2)/sqrt(length(BurstStarts));
        meanTriggeredBeta = mean(TriggeredBeta,2);
        semTriggeredBeta = std(TriggeredBeta,[],2)/sqrt(length(BurstStarts));
    else
    end
    
    %% Stores all data of interest as a non-scalar structure, making it easy
    %to look at and compare data from different sessions.
    %First the parameters and settings.
    BTOData(ind0).Session = session;
    %Then Burst Data
    BTOData(ind0).NumberOfBetaBursts = NumberOfBetaBursts;
    if betaburstnumber>0
        BTOData(ind0).MeanBurstDelta = meanBurstDelta;
        BTOData(ind0).SEMBurstDelta = semBurstDelta;
        BTOData(ind0).MeanTriggeredDelta = meanTriggeredDelta;
        BTOData(ind0).SEMTriggeredDelta = semTriggeredDelta;
        
        BTOData(ind0).MeanBurstTheta = meanBurstTheta;
        BTOData(ind0).SEMBurstTheta = semBurstTheta;
        BTOData(ind0).MeanTriggeredTheta = meanTriggeredTheta;
        BTOData(ind0).SEMTriggeredTheta = semTriggeredTheta;
        
        BTOData(ind0).MeanBurstAlpha = meanBurstAlpha;
        BTOData(ind0).SEMBurstAlpha = semBurstAlpha;
        BTOData(ind0).MeanTriggeredAlpha = meanTriggeredAlpha;
        BTOData(ind0).SEMTriggeredAlpha = semTriggeredAlpha;
        
        BTOData(ind0).MeanBurstBeta = meanBurstBeta;
        BTOData(ind0).SEMBurstBeta = semBurstBeta;
        BTOData(ind0).MeanTriggeredBeta = meanTriggeredBeta;
        BTOData(ind0).SEMTriggeredBeta = semTriggeredBeta;
    else
    end
    
    clearvars -except Folder DetectorChannel ComparatorChannel ind0 BTOData foldername subfolders animal channel1name channel2name
    
end

%Save the finished data structure as a .mat file for later use
BTOFolder = strcat(foldername(1:end-2),'Analysed\','BTOData');
mkdir(BTOFolder);
cd(BTOFolder);
save(strcat('BTOData','-',(animal),'-',(channel1name),(channel2name),'.mat'),'BTOData');

end

