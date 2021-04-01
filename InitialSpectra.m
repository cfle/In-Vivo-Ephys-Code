%Speed-binned LFP analysis using the SpeedSnip method. Suitable for one
%channel at a time.

function [InitialSpectraStructure] = InitialSpectraWYW(Folder,Channel)

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

    animalnumber = foldername(end-1:end);
    animal = strcat('Mouse',animalnumber);
    session = files.folder;
    session = session(end-4:end);
    
    [~,~,extension] = fileparts(channels{1});
    signals = [];
    
    for indch = 1
        if strcmp(extension,'.continuous')
            %Performs tracking, calculates speed and total distance travelled in m
            [x,y,~,speed_out,distance] = TrackingNovelFamiliar(subfolders{ind0},1,'two');
            tsf = fillmissing(speed_out.ts,'linear');
            speed = speed_out.speed;
            %Imports the .continuous data into MatLab
            [signals(:,indch),~,info] = load_open_ephys_data_faster(channels{indch});
        elseif strcmp(extension,'.mat')
            [x,y,ts,speed_out,distance] = TrackingNovelFamiliarSplitSession(subfolders{ind0});
            tsf = fillmissing(speed_out.ts,'linear');
            speed = speed_out.speed;
            %Imports the .continuous data into MatLab
            load(channels{indch});
            signals(:,indch) = sig;
            info.header.sampleRate = 30000;
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
    downsampleX = 30;
    signalsds = downsample(signals,downsampleX);
    signaldsdt = detrend(signalsds);
    
    %Set the parameters for Chronux analysis
    Fs = Fs/downsampleX;
    params.Fs = Fs;
    params.pad = 2;
    params.tapers = [2 3];
    params.err = 0;
    params.trialave = 0;
    params.segave = 0;
    
    %Performs multi-taper spectral analysis, firstly, for the unbinned
    %signal (ie the whole session) then for the binned signal, for each
    %frequency band. As it is not possible to manually set the bandwidth,
    %this allows us to cut off each band at the exact start and end.
    params.fpass = [1 120];
    [SpecgFull,t,fF] = mtspecgramc(signaldsdt,[1 1],params);
    SpecgFull = 10*log10(SpecgFull);

    MeanSpec = mean(SpecgFull(1:60,:),1);
    
    MeanSpec(fF>48 & fF<52) = NaN;
    MeanSpec = fillmissing(MeanSpec,'linear');
%     params.smoothfactorlow = length(fF(fF>=3 & fF<30))/9;
%     MeanCoh(fF>=3 & fF<30) = smoothdata(MeanCoh(fF>=3 & fF<30),'gaussian',params.smoothfactorlow);
%     params.smoothfactorhigh = length(fF(fF>=30 & fF<=120))/6;
%     MeanCoh(fF>=30 & fF<=120) = smoothdata(MeanCoh(fF>=30 & fF<=120),'gaussian',params.smoothfactorhigh);
    
    fD = fF(fF>=1 & fF<5);
    fT = fF(fF>=1 & fF<12);
    fA = fF(fF>=12 & fF<20);
    fB = fF(fF>=20 & fF<30);
    fLG = fF(fF>=30 & fF<65);
    fHG = fF(fF>=65 & fF<=120);
    MeanDelta = MeanSpec(fF>=1 & fF<5);
    MeanTheta = MeanSpec(fF>=5 & fF<12);
    MeanAlpha = MeanSpec(fF>=12 & fF<20);
    MeanBeta = MeanSpec(fF>=20 & fF<30);
    MeanLowGamma = MeanSpec(fF>=30 & fF<65);
    MeanHighGamma = MeanSpec(fF>=65 & fF<=120);
    
    clearvars PowerD PowerT PowerA PowerB PowerG PowerLG PowerHG
    clearvars FreqD FreqT FreqA FreqB FreqG FreqLG FreqHG
    
    %Calculates the raw power of each speed bin as the mean power in each
    %frequency band, and the freqency, as the frequency at which the power
    %is maximal.
    PowerD = mean(MeanDelta);
    PowerT = mean(MeanTheta);
    PowerA = mean(MeanAlpha);
    PowerB = mean(MeanBeta);
    PowerLG = mean(MeanLowGamma);
    PowerHG = mean(MeanHighGamma);
    
    PeakDC = max(MeanDelta,[],2);
    PeaksDC = find(MeanDelta==PeakDC);
    FreqD = median(fD(PeaksDC));
    PeakTC = max(MeanTheta,[],2);
    PeaksTC = find(MeanTheta==PeakTC);
    FreqT = median(fT(PeaksTC));
    PeakAC = max(MeanAlpha,[],2);
    PeaksAC = find(MeanAlpha==PeakAC);
    FreqA = median(fA(PeaksAC));
    PeakBC = max(MeanBeta,[],2);
    PeaksBC = find(MeanBeta==PeakBC);
    FreqB = median(fB(PeaksBC));
    PeakLGC = max(MeanLowGamma,[],2);
    PeaksLGC = find(MeanLowGamma==PeakLGC);
    FreqLG = median(fLG(PeaksLGC));
    PeakHGC = max(MeanHighGamma,[],2);
    PeaksHGC = find(MeanHighGamma==PeakHGC);
    FreqHG = median(fHG(PeaksHGC));
    
    %Stores all data of interest as a non-scalar structure, making it easy
    %to look at and compare data from different sessions.
    InitialSpectraStructure(ind0).Session = session;
    InitialSpectraStructure(ind0).Parameters = params;
    InitialSpectraStructure(ind0).fF = fF;
    InitialSpectraStructure(ind0).MeanSpectrum = MeanSpec;
    InitialSpectraStructure(ind0).SpecgFull = SpecgFull;
    InitialSpectraStructure(ind0).X = x;
    InitialSpectraStructure(ind0).Y = y;
    InitialSpectraStructure(ind0).Distance = distance;
    InitialSpectraStructure(ind0).Speed = speed;
    InitialSpectraStructure(ind0).DeltaPower = PowerD;
    InitialSpectraStructure(ind0).ThetaPower = PowerT;
    InitialSpectraStructure(ind0).AlphaPower = PowerA;
    InitialSpectraStructure(ind0).BetaPower = PowerB;
    InitialSpectraStructure(ind0).LowGammaPower = PowerLG;
    InitialSpectraStructure(ind0).HighGammaPower = PowerHG;
    InitialSpectraStructure(ind0).DeltaFrequency = FreqD;
    InitialSpectraStructure(ind0).ThetaFrequency = FreqT;
    InitialSpectraStructure(ind0).AlphaFrequency = FreqA;
    InitialSpectraStructure(ind0).BetaFrequency = FreqB;
    InitialSpectraStructure(ind0).LowGammaFrequency = FreqLG;
    InitialSpectraStructure(ind0).HighGammaFrequency = FreqHG;
    
end

%Session Analysis
MeanSpec = cell2mat({InitialSpectraStructure.MeanSpectrum}');
NovelSpecs = MeanSpec(1:8:9,:);
FamiliarSpecs = MeanSpec([2 3 4 5 6 7 8 10],:);
MeanNovelSpec = mean(NovelSpecs);
MeanFamiliarSpec = mean(FamiliarSpecs);
NovelStd = std(NovelSpecs,0,1);
FamiliarStd = std(FamiliarSpecs,0,1);
NovelSEM = NovelStd/sqrt(size(NovelSpecs,1));
FamiliarSEM = FamiliarStd/sqrt(size(FamiliarSpecs,1));

% %Creates the mean spectrums for Novel and Familiar with SEM.
Spectrums = strcat(foldername,'\','Spectrums');
mkdir(Spectrums);
cd(Spectrums);
figure2name = strcat(animal,'-',Channel,'-','InitialSpectrumFamiliarity');
figure('Name',figure2name);
NovelSpec = shadedErrorBar(fF,MeanNovelSpec,NovelSEM);
set(NovelSpec.edge,'LineStyle','none')
NovelSpec.mainLine.LineWidth = 1.5;
NovelSpec.mainLine.Color = [1 0 0];
NovelSpec.patch.FaceColor = [0.6 0 0];
hold on;
FamiliarSpec = shadedErrorBar(fF,MeanFamiliarSpec,FamiliarSEM);
set(FamiliarSpec.edge,'LineStyle','none')
FamiliarSpec.mainLine.LineWidth = 1.5;
FamiliarSpec.mainLine.Color = [0 0 1];
FamiliarSpec.patch.FaceColor = [0 0 0.6];
alpha(0.3);
set(gca,'LineWidth',1.5);
set(gca,'FontSize',12)
xlabel('Frequency (Hz)','fontsize',14);
ylabel('Power (DB)','fontsize',14);
legend({'Novel','Familiar'});
legend boxoff
set(gca,'box','off');
hgsave(figure2name)

%Save the finished data structure as a .mat file for later use
SpectraFolder = strcat(foldername(1:end-2),'Analysed\','Spectra');
mkdir(SpectraFolder);
cd(SpectraFolder);
save(strcat('InitialSpectra','-',(animal),'-',Channel,'.mat'),'InitialSpectraStructure');

end

