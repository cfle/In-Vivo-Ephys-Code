function [PACData] = PhaseAmplitudeCoupling(Folder,Channel)

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
    
    %Imports the .continuous data into MatLab
    signal = [];
    [signal(:,1),~,info] = load_open_ephys_data_faster(channels{1});
    
    %Cuts the signals to the experiment duration ie 15minutes in this case
    ExpDurationMinutes = 15;
    Fs = info.header.sampleRate;
    ExpDurationFrames = (ExpDurationMinutes*60)*Fs;
    ExpDurationSeconds = ExpDurationMinutes*60;
    signal = signal(1:ExpDurationFrames,:);
    
    %Downsamples the signals and timestamps, to reduce computation time,
    %and detrends to remove any baseline fluctuations
    downsampleX = 10;
    Fs = Fs/downsampleX;
    signalsds = downsample(signal,downsampleX);
    signalsdsdt = detrend(signalsds);
    signalsdsdt(isnan(signalsdsdt)) = [];
   
    signalsdsdt = signalsdsdt(1:60*Fs);
    
    %Sets the theta and gamma frequency bins to investigate, and
    %preallocates the modulation index matrix.
    phasebinstart = 2;
    phasebinstep = 0.25;
    phasebinstop = 12;
    phasebins = (phasebinstart:phasebinstep:phasebinstop);
    numpbins = (length(phasebins)-1);
    amplitudebinstart = 10;
    amplitudebinstep = 2;
    amplitudebinstop = 100;
    amplitudebins = (amplitudebinstart:amplitudebinstep:amplitudebinstop);
    numabins = (length(amplitudebins)-1);
    
    MI = zeros(numel(amplitudebins)-1,numel(phasebins)-1);
    for indpf = 1:numpbins
        phasefilter{indpf} = designfilt('bandpassiir','FilterOrder',2,'HalfPowerFrequency1',phasebins(indpf),'HalfPowerFrequency2',phasebins(indpf+1),'DesignMethod','butter','SampleRate',Fs);
    end
    for indaf = 1:numabins
        amplitudefilter{indaf} = designfilt('bandpassiir','FilterOrder',2,'HalfPowerFrequency1',amplitudebins(indaf),'HalfPowerFrequency2',amplitudebins(indaf+1),'DesignMethod','butter','SampleRate',Fs);
    end
    
    %For the calculation of a normalised MI value, a number of surrogate
    %MI values are calculated by circularly shifting the signals with
    %random magnitude delays. Here, the number of surrogates calculated, as
    %well as the maximum and minimum delays are given. Fs is used to
    %prevent the data being shifted only one second forward or backward.
    %This method for calculating MI was adapted from (Canolty et al. 2006),
    %but changes were made to the script, including vectorisation.
    NumberOfPoints = length(signalsdsdt);
    NumberOfSurrogates = 10;
    MinSkip = Fs;
    MaxSkip = NumberOfPoints-Fs;
    
    clear surrogate_amplitude
    for indp = 1:numpbins
        phasesignal = filtfilt(phasefilter{indp},signalsdsdt);
        phase = angle(hilbert(phasesignal));
        for inda = 1:numabins
            amplitudesignal = filtfilt(amplitudefilter{inda},signalsdsdt);
            amplitude = abs(hilbert(amplitudesignal));
            mraw = abs(mean(amplitude.*exp(1i*phase)));
            
            % compute surrogate values
            skip = randi([MinSkip MaxSkip],[NumberOfSurrogates,1]);
            surrogateamplitude = zeros(length(amplitude),NumberOfSurrogates);
            for inds = 1:NumberOfSurrogates
                surrogateamplitude(:,inds) = circshift(amplitude,-skip(inds));
            end
            surrogateM = abs(mean(surrogateamplitude(:,1:NumberOfSurrogates).*exp(1i.*phase)));
            
            %Fits gaussian to surrogate data and normalize length using surrogate data (z-score)
            [surrogatemean,surrogatestd] = normfit(surrogateM);
            mnormlength = (mraw-surrogatemean)/surrogatestd;
            mnormphase = angle(mraw);
            mnorm = mnormlength*exp(1i*mnormphase);
            MI(inda,indp) = abs(mnorm);
        end
    end
    
    pbt = (phasebinstart+phasebinstep/2:phasebinstep:phasebinstop-phasebinstep/2);
    abt = (amplitudebinstart+amplitudebinstep/2:amplitudebinstep:amplitudebinstop-amplitudebinstep/2);
    numphaseticks = numpbins;
    numamplitudeticks = numabins;
    interpfactor = 2;
    pbi = phasebinstart:phasebinstep/interpfactor:phasebinstop;
    abi = amplitudebinstart:amplitudebinstep/interpfactor:amplitudebinstop;
    pbti = (phasebinstart+(phasebinstep/2)/interpfactor:phasebinstep/interpfactor:phasebinstop-(phasebinstep/2)/interpfactor);
    abti = (amplitudebinstart+(amplitudebinstep/2)/interpfactor:amplitudebinstep/interpfactor:amplitudebinstop-(amplitudebinstep/2)/interpfactor);
    numphaseticksint = interpfactor*numphaseticks;
    numamplitudeticksint = interpfactor*numamplitudeticks;
    
    phaseticksi = linspace(pbt(1),pbt(end),numphaseticksint);
    amplitudeticksi = linspace(abt(1),abt(end),numamplitudeticksint);
    [PT,AT] = ndgrid(pbt,abt);
    F = griddedInterpolant(PT,AT,MI');
    [PTI,ATI] = ndgrid(pbti,abti);
    MIinterp = F(PTI,ATI)';
    
    %%
    %Sets the theta and gamma frequency bins to investigate, and
    %preallocates the modulation index matrix.
    phasefocus = [4 12];
    amplitudefocus = [30 120];
    phasebars = (0:10:360);
    phasebarstr = string(phasebars);
    phasefocusfilter = designfilt('bandpassiir','FilterOrder',2,'HalfPowerFrequency1',phasefocus(1),'HalfPowerFrequency2',phasefocus(2),'DesignMethod','butter','SampleRate',Fs);
    amplitudefocusfilter = designfilt('bandpassiir','FilterOrder',2,'HalfPowerFrequency1',amplitudefocus(1),'HalfPowerFrequency2',amplitudefocus(2),'DesignMethod','butter','SampleRate',Fs);
    
    phasefocussignal = filtfilt(phasefocusfilter,signalsdsdt);
    focusphase = angle(hilbert(phasefocussignal));
    amplitudefocussignal = filtfilt(amplitudefocusfilter,signalsdsdt);
    focusamplitude = abs(hilbert(amplitudefocussignal));
    phasedeg = rad2deg(focusphase)+180;
    [phasesorted,phasesort] = sort(phasedeg);
    amplitudesorted = focusamplitude(phasesort);
    for izz3 = 1:36
        meanamplitude(izz3) = mean(amplitudesorted(phasesorted >= phasebars(izz3) & phasesorted < phasebars(izz3+1)));
    end
    totalamplitude = sum(meanamplitude);
    meanamplitude = meanamplitude./totalamplitude;
    meanmeanamplitude = mean(meanamplitude);
    meanamplituderep = [meanamplitude meanamplitude];

    %%
    
    TPACFigureFolder = strcat(foldername,'\','ThetaPACHeatMaps');
    mkdir(TPACFigureFolder);
    cd(TPACFigureFolder);
    figurename = strcat(animal,'-',session(end-1:end),'-',channelname,'-','ThetaPAC');
    figure('Name',figurename);
    imagesc(pbti,abti,MIinterp);
    axis xy
    set(gca,'LineWidth',1.5);
    set(gca,'FontSize',12);
    xLabel = xlabel('Phase Frequency (Hz)','fontsize',14);
    yLabel = ylabel('Amplitude Frequency (Hz)','fontsize',14);
    set(gca,'box','off');
    grid off
    CB = colorbar;
    caxis([0 10]);
    CBtitle = get(CB,'Title');
    set(CBtitle,'String','MI');
    CBpos = get(CBtitle,'Position');
    set(CBtitle,'Position',[CBpos(1) CBpos(2)-320]);
    set(CBtitle,'FontSize',14);
    hgsave(figurename)
    
    PACData(ind0).PhaseBins = phasebins;
    PACData(ind0).AmplitudeBins = amplitudebins;
    PACData(ind0).NumberOfSurrogates = NumberOfSurrogates;
    PACData(ind0).MI = MI;
    PACData(ind0).PhaseBinsinterp = pbi;
    PACData(ind0).AmplitudeBinsinterp = abi;
    PACData(ind0).MIinterp = MIinterp;
    
end

%Save the finished data structure as a .mat file for later use

PACFolder = strcat(foldername(1:end-2),'Analysed\','PACData');
mkdir(PACFolder);
cd(PACFolder);
save(strcat('PACData','-',(animal),'-',(channelname),'.mat'),'PACData');

end
