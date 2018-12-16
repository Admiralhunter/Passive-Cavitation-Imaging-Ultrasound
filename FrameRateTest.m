% PCI beamforming code for use with Verasonics. Requires both the receive
% data file and the Verasonics .mat File that defines all of the parameter
% settings.
% 20131209: Modified for Vantage code that can use units of 'mm'; also add
% in a frequency-sum beamforming method based on Abadi 2013 POMA. Also
% converted data to all doubles rather than singles because the maximum
% number handled by single was insufficiently large when dealing with
% a 4x frequency-sum
% 20140106: Based on PCIVS_v04, which included frequency-sum beamforming.
% Only do quadratic frequency-sum.
%20170605: This script is specifically for post processing ultrasound data.
%This can not run in real time, and instead will run through a folder
%containing the necessary data. After processing a specific data file the
%program will save a .mat file with the workspace variables saved as well
%as a .png figure for the subharmonic as well as the Fundamental plot.
% Written by KJH & HAP
%NOTE: Requires function "Analyze_Frame.m" on the computer to be called.
%Ensure to download this crucial file from appropiate folder. Ask Kevin
%where it's located.
%20170717: Added code to convert RData from digitized units (14-Bit) to volts

%xres = xlsread('FrameRateVSNumPixels.xlsx','B2:B61');
%zres = xlsread('FrameRateVSNumPixels.xlsx','C2:C61');
%numoffreq = xlsread('FrameRateVSNumPixels.xlsx','D2:D61');
xres = 150;
zres = 25;
numoffreq = 2720; 
averageframerate = 5;
framerate = zeros(1,100000);
slowholder = zeros(1,100000);
xresholder = zeros(1,100000);
zresholder = zeros(1,100000);
freqresholder = zeros(1,100000);
g = gpuDevice;

counter = 0;
while xres <= 500
    counter = counter +1;
    
    if (numoffreq < 2800) && (averageframerate > 1) && (averageframerate ~= 0) && ((numoffreq+ 100) < 2800)
        numoffreq = numoffreq +100;
    elseif (numoffreq < 2800) && (averageframerate < 1) && (averageframerate > 0.5) && (exist('numoffreqholder','var') == 0) && (averageframerate ~= 0)
        numoffreqholder = linspace(numoffreq,2800,10);
        numoffreq = numoffreqholder(2);
        counter2 = 3;
    elseif (numoffreq < 2800) && (averageframerate < 1) && (averageframerate > 0.5) && (exist('numoffreqholder','var') == 1) && (counter2 < 11) && (averageframerate ~= 0)
        numoffreq = numoffreqholder(counter2);
        counter2 = counter2 + 1;
    else
        numoffreq = 20;
        zres = zres + 25;
        clear numoffreqholder
        counter2 = 3;
    end
    
    if zres == 500
        zres = 25;
        xres = xres + 25;
    end
        
        xresholder(counter) = xres;
        zresholder(counter) = zres;
        freqresholder(counter) = numoffreq;
        
    
    
    
    tic
    
    % Establish
    
    fc = 5;
    
    %folder_name = uigetdir('C:\Users\palcicha\Documents\Data_to_run');
    BASEFolder = 'L:\Shared\LabUsers\Niusha\20170727\Trial 1\BD 317.80 (1)\';
    VSXFolder ='L:\Shared\LabUsers\Niusha\20170727\';
    
    
    DataFile = 'PData_28156';
    VerasonicsFile = 'L7-4_PCI_w128RyLns_vAW1c';
    nfreq = [];
    
    %Choose range of desired frequencies MHz
    
    
    [YY,MM,DD,hh,mm,ss] = datevec(now);
    
    % Load Data and Verasonics files
    %load([BASEFolder DataFile]);
     RData = single(RData);  % RData = single(RData);
    %load([VSXFolder VerasonicsFile]);
    LoadingTime = toc;
    
    tic
    RDataWidth = 1:Resource.Parameters.numRcvChannels;
    RDataLength = 1:(2*(Receive(1).endDepth - Receive(1).startDepth)*4);
    RData = RData(RDataLength,RDataWidth,:);
    
    %Code to convert digitized signal to volts, assuming a flat TGC
    if ~exist('RcvProfile','var')
        RcvProfile.LnaGain = 18;
        RcvProfile.PgaGain = 24;
    else
        if ~isfield(RcvProfile, 'LnaGain')   %If Low Noise Amplifier gain is undefined, value is set to 18dB (default)
            RcvProfile.LnaGain = 18;
        end
        if ~isfield(RcvProfile, 'PgaGain')   %If Gain Buffer Amplifier gain is undefined, value is set to 24dB (default)
            RcvProfile.PgaGain = 24;
        end
    end
    TG = RcvProfile.LnaGain + RcvProfile.PgaGain + (40/1023*mean(TGC(1).CntrlPts) - 40) - 0.6; %Compute total gain in dB
    RData = single(RData).*(1/16384).*(1/10^(TG/20));    %Conversion to volts
    
    
    % TX/RX/Material Parameters
    %A = TW(1).Parameters(1);                   % lambda/2 in number of clock cycles
    SampPerWave = Receive(1).samplesPerWave;    % Receive.samplesperWave... The sample rate for stored samples is 1, 2, or 4 samples per wavelength of the specified transducer frequency, Trans.frequency. The maximum bandwidth (4 samps/wavelength) is 200% of Fc, the transducer center frequency. The maximum value of Receive.endDepth - Receive.startDepth is 1024 wavelengths.
    %Fs = 180/(2*A)*SampPerWave;                % Fs in MHz
    Fs = Trans.frequency.*SampPerWave;          % Fs in MHz
    SoS = Resource.Parameters.speedOfSound/1000;    % Speed of sound mm/us
    NFrameStart = 1;                            % Frame to start analysis
    %NFramesAnalyze = 1;
    NFramesAnalyze = 5;% number of frames to analyze
    if strcmp(Trans.units,'mm')
        Lambda_mm = 1;  % if data is already in units of mm, then set scaling factor to unity
    else
        Lambda_mm = SoS./Trans.frequency;           % if data is no in mm, use this scaler to convert wavelength in mm
    end
    
    % Transducer Parameters
    ElementPitch = Trans.spacing*Lambda_mm;         % in mm
    Aperture = Trans.ElementPos(:,1)*Lambda_mm;     % Element positions in the array in mm
    NChannels = length(Aperture);
    zarray = 0;         % range location of array (should always be zero with a linear array)
    
    % PCI Dimensions in mm
    % x = (0:size(RData,2)-1)*ElementPitch;   % Lateral (x) pixel locations to beamform, in mm
    % x = x-x(round(length(x)/2));            % Set the center of the image to have a value of x=0;
    % x = Aperture;                           % Lateral (x) pixel locations to beamform, in mm
    
    x = linspace(0,15,xres); % in mm
    sizex = length(x);
    
    z = linspace(25,35,zres); % Range  (z)  pixel locations to beamform, in mm
    sizez = length(z);
    

    
    
    % Mean Subtract to remove DC component on each channel
    RData_DC = mean(RData,1);
    RData = RData-repmat(RData_DC,[size(RData,1) 1 1]);
    NFrames = size(RData,3);    % Number of frames
    
    % Zeropad the data so that when it is shifted it does not loop around.
    % zpad = 500; % This value is set imperically, should set it based on max shifts to be performed when computing 'tof'
    % RData = [zeros(zpad,NChannels,NFrames); RData; zeros(zpad,NChannels,NFrames);
    
    
    % Timing and frequency info
    dt = 1/Fs;                  % in us
    T = (size(RData,1)-1)*dt;   % in us
    Time = (0:dt:T)';
    % f = linspace(0,Fs,size(RData,1));
    if isempty(nfreq)
        nfreq = size(RData,1);
    end
    fcpu = linspace(-Fs/2,Fs/2,nfreq); %linspace(-Fs/2,Fs/2,size(RData,1));
    f = gpuArray(single(fcpu));
    f = f(1:numoffreq);
    %f1 is the lower end of the frequency spectrum you want to analyze
    %f2 is the upper end of the frequency spectrum you want to analyze
    
    f1 = -20;
    f2 = 20;
    
    fselecti = find(f > f1 & f < f2);
    
    fprintf('The length of x is %i\n',sizex)
    fprintf('The length of z is %i\n',sizez)
    fprintf('The length of f is %i\n',length(f))
    
    
    % Data Parameter size
    % EndData = find(squeeze(sum(abs(RData(:,:,1)),2))==0);   % end of recorded data (note the buffer is larger than this)
    % Nrange=EndData(1)-1;        % Number of range points recorded (receive buffer is larger)
    Nlines = size(RData,2);     % Number of data channels (i.e. elements)
    % Nlines = length(Aperture);
    
    
    
    % Distance from point source with position (x(i),z(k)) along the x axis to
    % element w and corresponding time of flight (with bulk shifts removed).
    % Additionally compute apodization matrix based on transmit angle
    d = single(zeros(length(z),length(x),length(Aperture)));
    tof = d;    costheta = d;
    for w = 1:Nlines            % Index over channels
        for i = 1:length(x)     % Index over lateral pixel index
            for k = 1:length(z) % Index over range pixel index
                d(k,i,w) = sqrt((z(k)-zarray)^2+( x(i)-Aperture(w))^2); % in mm
                tof(k,i,w) = d(k,i,w)./SoS;     % in us
                %tof(k,i,w) = tof(k,i,w) - z(k)./SoS;    % remove bulk shift
                costheta(k,i,w) = z(k)./d(k,i,w); % compute angle between each element and each image point, use for apodization matrix
            end
        end
    end
    
    % Perform PCI beamforming
    AvePowerSpectraAllFrames = gpuArray(single(zeros([length(x) length(z) length(fselecti)])));
    
    DeclarationTime = toc;
    
    
    
    
    %Creates gpuArrays for data to be processed
    tofgpu = gpuArray(single(tof)); fgpu = gpuArray(f(fselecti));
    RData = gpuArray(single(RData));
    RFFT = fftshift(fft(RData(:,:,:),[],1),1);
    
    RFFTgpu = RFFT(fselecti,:,1:NFramesAnalyze);
    
    %AvePowerSpectraAllFrames2gpu = AvePowerSpectraAllFramesgpu;
    
    FrameComputationTimegpu = zeros(1,length(NFrameStart:NFramesAnalyze));
    
    f1 = f(fselecti);
    f = f1;
    %Function to run PCI Beam Forming
    [AvePowerSpectraAllFrames_g,averageframerate,large] = Analyze_FrameV3_20170707_test(x,z,RFFTgpu,NFrameStart,NFramesAnalyze,tofgpu,fgpu,f);
    framerate(counter) = averageframerate;
    slowholder(counter) = large;
    %xlswrite('FrameRateVSNumPixels.xlsx',averageframerate,'Sheet1',hzrange)
    %xlswrite('FrameRateVSNumPixels.xlsx',large,'Sheet1',slowrange)
    
    %[AvePowerSpectraAllFrames_g,FocalWaveformWin] = Analyze_FrameV5(x,z,RFFTgpu,NFrameStart,NFramesAnalyze,tofgpu,fgpu,f);
    %clear RData RFFT tofgpu tof d costheta AvePowerSpectraAllFrames AvePowerSpectraAllFramesgpu f1  g RFFTgpu
    
    fgpu = gather(fgpu);
    f = gather(f);
    fselecti = gather(fselecti);
    
    
    window = '_chebwin_300';
    
    
    
end
% save([BASEFolder DataFile(1:end) '_PCI_' num2str(YY) num2str(MM) num2str(DD) num2str(hh) num2str(mm) num2str(round(ss)) 'gpu' window '.mat'])
% disp('Done Saving')
% 
% 
% 
% FileName = DataFile;
% FileNameEnd = length(FileName);
% FileName = FileName(1:FileNameEnd);
% FileNameUlt = strcat(FileName,'_Ultraharmonic',window);
% FileNameFun = strcat(FileName,'_Fundamental',window);
% FileNameHarm = strcat(FileName,'_Harmonic',window);
% FileNameIH = strcat(FileName,'_Inharmonic',window);
% FileNameSpec = strcat(FileName,'_FocSpectrum',window);
% FileNameWave = strcat(FileName,'_FocWaveform',window);
% FileNameWaveInd = strcat(FileName,'_FocWaveformInd',window);
% 
% 
% DataFile = strrep(FileName,'_',' ');
% 
% %Ultraharmonic
% fsub = closest(f,3*fc/2); PCI_Sub = sum(AvePowerSpectraAllFrames_g(:,:,fsub),3);
% PCI_Subnorm = PCI_Sub./max(PCI_Sub(:));
% h_spec = figure(1); imagesc(x,z,10*log10(squeeze(PCI_Subnorm))'); colorbar; colormap(hot);
% caxis([-30 0]); axis equal;
% AxesLabelsNTitles(gca,'Lateral (mm)','Range (mm)',sprintf('Ultraharmonic PCI'),16,'black')
% %print(fullfile(BASEFolder,FileNameUlt),'-dpng');
% savefig(h_spec,fullfile(BASEFolder,FileNameUlt));
% 
% 
% % Fundamental
% ffund = closest(f,fc); PCI_Fund = sum(AvePowerSpectraAllFrames_g(:,:,ffund),3);
% PCI_Fundnorm = PCI_Fund./max(PCI_Fund(:));
% h_spec = figure(2); imagesc(x,z,10*log10(squeeze(PCI_Fundnorm))'); colorbar; colormap(hot);
% caxis([-30 0]); axis equal;
% AxesLabelsNTitles(gca,'Lateral (mm)','Range (mm)',sprintf('Fundamental PCI'),16,'black')
% %print(fullfile(BASEFolder,FileNameFun),'-dpng');
% savefig(h_spec,fullfile(BASEFolder,FileNameFun));
% 
% % second harmonic
% fharm = closest(f,2*fc); PCI_Harm = sum(AvePowerSpectraAllFrames_g(:,:,fharm),3);
% PCI_Harmnorm = PCI_Harm./max(PCI_Harm(:));
% h_spec = figure(3); imagesc(x,z,10*log10(squeeze(PCI_Harmnorm))'); colorbar; colormap(hot);
% caxis([-30 0]); axis equal;
% AxesLabelsNTitles(gca,'Lateral (mm)','Range (mm)',sprintf('2nd Harm PCI'),16,'black')
% %print(fullfile(BASEFolder,FileNameHarm),'-dpng');
% savefig(h_spec,fullfile(BASEFolder,FileNameHarm));
% 
% % inharmonic
% fIH = closest(f,6.9); PCI_IH = sum(AvePowerSpectraAllFrames_g(:,:,fIH),3);
% PCI_IHnorm = PCI_IH./max(PCI_IH(:));
% h_spec = figure(4); imagesc(x,z,10*log10(squeeze(PCI_IHnorm))'); colorbar; colormap(hot);
% caxis([-30 0]); axis equal;
% AxesLabelsNTitles(gca,'Lateral (mm)','Range (mm)',sprintf('Inharmonic PCI'),16,'black')
% %print(fullfile(BASEFolder,FileNameIH),'-dpng');
% savefig(h_spec,fullfile(BASEFolder,FileNameIH));
% 
% clear h_spec
% 
% h_spec = figure(5); plot(f,10*log10(FocalPowerSpectra)); xlim([0 14]);
% AxesLabelsNTitles(gca,'Frequency (MHz)','Power Spectra (dB)',sprintf('PCI Max Loc Spectra',DataFile),16,'black')
% %print(fullfile(BASEFolder,FileNameSpec),'-dpng');
% savefig(h_spec,fullfile(BASEFolder,FileNameSpec));
% 
% clear h_spec
% 
% h_spec = figure(6); plot(Time,FocalWaveform, Time,Window, Time, FocalWaveformWindowed); %xlim([0 14]);
% AxesLabelsNTitles(gca,'Time (us)','Amplitude',sprintf('PCI Max Loc Waveforms'),16,'black')
% %print(fullfile(BASEFolder,FileNameSpec),'-dpng');
% savefig(h_spec,fullfile(BASEFolder,FileNameWave));
% 
% clear h_spec
% 
% h_spec = figure(7); plot(1:length(Time),FocalWaveform, 1:length(Time),Window, 1:length(Time), FocalWaveformWindowed); %xlim([0 14]);
% AxesLabelsNTitles(gca,'Time (index)','Amplitude',sprintf('PCI Max Loc Waveforms'),16,'black')
% %print(fullfile(BASEFolder,FileNameSpec),'-dpng');
% savefig(h_spec,fullfile(BASEFolder,FileNameWaveInd));
