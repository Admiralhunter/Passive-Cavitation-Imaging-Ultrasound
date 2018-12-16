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
%20180321: Auto Window placement had been added by utilizing kevin's code.
%Cumsum is used to find the middle of the signal and then a desired window
%length is placed centered at the index.


clc; clear;close all
g = gpuDevice;

tic

% Establish

fc = 5;

StartFolder = 'D:\20161220 AW-20180409T012856Z-001\';
VerasonicsFolder = 'D:\20161220 AW-20180409T012856Z-001\20161220 AW\';
SaveFolder = 'D:\20161220 AW-20180409T012856Z-001\20161220 AW\Analyzed_WindowLength240';
VerasonicsFile = 'L7-4_PCI_w128RyLns_v06_78125_NewMainVH6';





D = dir([StartFolder]);
Folders = D([D.isdir]);

nfreq = [];

counter = 0;
Folders = D([D.isdir]);
for ii = 1:length(Folders)
    if ii > 2
        
        DirFile = dir(strcat(StartFolder,Folders(ii).name));
        for df= 3:length(DirFile)%1:length(FileNames)
            counter = counter +1;
            BaseFolder = strcat(StartFolder,Folders(ii).name,'\',DirFile(df).name,'\');
            Files = dir(BaseFolder);
            for ff = 1:length(Files)
                if length(strfind(Files(ff).name,'PData')) > 0 && length(strfind(Files(ff).name,'B')) == 0
                    %Selects next data file from folder
                    FileName = strcat(BaseFolder,Files(ff).name);
                    % clear all
                    clearvars -except df StartFolder Folders DirFile Files BaseFolder counter NumofFiles FileName VerasonicsFolder VerasonicsFile ii PulseDurCyc PulseFreq PulseArrival WindowFlag BASEFolder DataFolder NFrameStart NFramesAnalyze SaveIndividFrames nfreq fc SaveFolder ff
                    
                    
                    DataFile = FileName;   % 5 MPa Deg: RData_66477; 5 MPa ADV: RData_71667
                    
                    
                    %Choose range of desired frequencies MHz
                    
                    
                    [YY,MM,DD,hh,mm,ss] = datevec(now);
                    
                    % Load Data and Verasonics files
                    load([DataFile]);
                    RData = double(RData);  % RData = single(RData);
                    load([VerasonicsFolder VerasonicsFile]);
                    LoadingTime = toc;
                    fprintf('Data loading took %08.4fs \n',LoadingTime)
                    
                    tic
                    RDataWidth = 1:Resource.Parameters.numRcvChannels;
                    RDataLength = 1:(2*(Receive(1).endDepth - Receive(1).startDepth)*4);
                    RData = squeeze(RData);
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
                    NFramesAnalyze = length(RData(1,1,:));                         % number of frames to analyze
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

                    x = linspace(-10,10,30);  % in mm
                    sizex = length(x);

                    z = linspace(25,35,50); % Range  (z)  pixel locations to beamform, in mm
                    sizez = length(z);
                    
                    fprintf('The length of x is %i\n',sizex)
                    fprintf('The length of z is %i\n',sizez)
                    
                    
                    % Mean Subtract to remove DC component on each channel
                    RData_DC = mean(RData,1);
                    RData = RData-repmat(RData_DC,[size(RData,1) 1 1]);
                    NFrames = size(RData,3);    % Number of frames
                    
                    % Zeropad the data so that when it is shifted it does not loop around.
                    % zpad = 500; % This value is set imperically, should set it based on max shifts to be performed when computing 'tof'
                    % RData = [zeros(zpad,NChannels,NFrames); RData; zeros(zpad,NChannels,NFrames)];
                    
                    
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
                    
                    %f1 is the lower end of the frequency spectrum you want to analyze
                    %f2 is the upper end of the frequency spectrum you want to analyze
                    f1 = -20;
                    f2 = 20;
                    fselecti = find(f > f1 & f < f2);
                    
                    
                    
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
                    AvePowerSpectraAllFrames = gpuArray(double(zeros([length(x) length(z) length(fselecti)])));
                    
                    DeclarationTime = toc;
                    fprintf('Initial declarations took %08.4fs\n',DeclarationTime)
      
                    
                    %Transfer data into gpuArray's for processing
                    tofgpu = gpuArray(single(tof)); fgpu = gpuArray(f(fselecti));
                    RFFT = fftshift(fft(RData(:,:,:),[],1),1);
                    RFFTgpu = single(gpuArray(RFFT(fselecti,:,:)));
                    
                    AvePowerSpectraAllFramesgpu = zeros([length(x) length(z) length(fselecti)],'gpuArray');

                    
                    FrameComputationTimegpu = zeros(1,length(NFrameStart:NFramesAnalyze));
                    
                    f1 = f(fselecti);
                    f = f1;
                    
                    %Function to Run PCI Beam Forming
                    [AvePowerSpectraAllFrames_g,FocalPowerSpectra,FocalWaveform,Window,FocalWaveformWindowed] = Analyze_FrameV3_20171221(x,z,RFFTgpu,NFrameStart,NFramesAnalyze,tofgpu,fgpu,f);
                    clear RData RFFT tofgpu tof d costheta AvePowerSpectraAllFrames AvePowerSpectraAllFramesgpu f1  g RFFTgpu
                    
                    fgpu = gather(fgpu);
                    f = gather(f);
                    fselecti = gather(fselecti);

                    window = '_chebwin_300';
                    matfile = fullfile(([SaveFolder '\' Folders(ii).name '_' DirFile(df).name '_' Files(ff).name  '_PCI_' num2str(YY) num2str(MM) num2str(DD) num2str(hh) num2str(mm) num2str(round(ss)) 'gpu' window '.mat']));
                    save(matfile)
                    disp('Done Saving')
                    
                    
                    
                    FileName = DataFile;
                    FileNameEnd = length(FileName)-4;
                    FileName = strcat(SaveFolder,'\',Folders(ii).name,'_',DirFile(df).name, '_', Files(ff).name);
                    FileNameUlt = strcat(FileName,'_Ultraharmonic',window,'.fig');
                    FileNameFun = strcat(FileName,'_Fundamental',window,'.fig');
                    FileNameHarm = strcat(FileName,'_Harmonic',window,'.fig');
                    FileNameIH = strcat(FileName,'_Inharmonic',window,'.fig');
                    FileNameSpec = strcat(FileName,'_FocSpectrum',window,'.fig');
                    FileNameWave = strcat(FileName,'_FocWaveform',window,'.fig');
                    FileNameWaveInd = strcat(FileName,'_FocWaveformInd',window,'.fig');
                    
                    
                    DataFile = strrep(FileName,'_',' ');
                    
                    %Ultraharmonic
                    fsub = closest(f,3*fc/2); PCI_Sub = sum(AvePowerSpectraAllFrames_g(:,:,fsub),3);
                    PCI_Subnorm = PCI_Sub./max(PCI_Sub(:));
                    h_spec = figure(1); imagesc(x,z,10*log10(squeeze(PCI_Subnorm))'); colorbar; colormap(hot);
                    caxis([-30 0]); axis equal;
                    AxesLabelsNTitles(gca,'Lateral (mm)','Range (mm)',sprintf('Ultraharmonic'),16,'black')
                    %print(fullfile(BASEFolder,FileNameUlt),'-dpng');
                    %savefig(h_spec,FileNameUlt);
                    
                    
                    % Fundamental
                    ffund = closest(f,fc); PCI_Fund = sum(AvePowerSpectraAllFrames_g(:,:,ffund),3);
                    PCI_Fundnorm = PCI_Fund./max(PCI_Fund(:));
                    h_spec = figure(2); imagesc(x,z,10*log10(squeeze(PCI_Fundnorm))'); colorbar; colormap(hot);
                    caxis([-30 0]); axis equal;
                    AxesLabelsNTitles(gca,'Lateral (mm)','Range (mm)',sprintf('Fundamental'),16,'black')
                    %print(fullfile(BASEFolder,FileNameFun),'-dpng');
                    %savefig(h_spec,FileNameFun);
                    
                    % second harmonic
                    fharm = closest(f,2*fc); PCI_Harm = sum(AvePowerSpectraAllFrames_g(:,:,fharm),3);
                    PCI_Harmnorm = PCI_Harm./max(PCI_Harm(:));
                    h_spec = figure(3); imagesc(x,z,10*log10(squeeze(PCI_Harmnorm))'); colorbar; colormap(hot);
                    caxis([-30 0]); axis equal;
                    AxesLabelsNTitles(gca,'Lateral (mm)','Range (mm)',sprintf('2nd Harm'),16,'black')
                    %print(fullfile(BASEFolder,FileNameHarm),'-dpng');
                    %savefig(h_spec,FileNameHarm);
                    
                    % inharmonic
                    fIH = closest(f,6.9); PCI_IH = sum(AvePowerSpectraAllFrames_g(:,:,fIH),3);
                    PCI_IHnorm = PCI_IH./max(PCI_IH(:));
                    h_spec = figure(4); imagesc(x,z,10*log10(squeeze(PCI_IHnorm))'); colorbar; colormap(hot);
                    caxis([-30 0]); axis equal;
                    AxesLabelsNTitles(gca,'Lateral (mm)','Range (mm)',sprintf('Inharmonic'),16,'black')
                    %print(fullfile(BASEFolder,FileNameIH),'-dpng');
                    %savefig(h_spec,FileNameIH);
                    
                    h_spec = figure(5); plot(f,10*log10(FocalPowerSpectra)); xlim([0 14]);
                    AxesLabelsNTitles(gca,'Frequency (MHz)','Power Spectra (dB)',sprintf('PCI Max Loc Spectra'),16,'black')
                    %print(fullfile(BASEFolder,FileNameSpec),'-dpng');
                    %savefig(h_spec,FileNameSpec);
                    
                    h_spec = figure(6); plot(Time,FocalWaveform, Time,Window, Time, FocalWaveformWindowed); %xlim([0 14]);
                    AxesLabelsNTitles(gca,'Time (us)','Amplitude',sprintf('PCI Max Loc Waveforms'),16,'black')
                    %print(fullfile(BASEFolder,FileNameSpec),'-dpng');
                    %savefig(h_spec,FileNameWave);
                    
                    h_spec = figure(7); plot(1:length(Time),FocalWaveform, 1:length(Time),Window, 1:length(Time), FocalWaveformWindowed); %xlim([0 14]);
                    AxesLabelsNTitles(gca,'Time (index)','Amplitude',sprintf('PCI Max Loc Waveforms'),16,'black')
                    %print(fullfile(BASEFolder,FileNameSpec),'-dpng');
                    %savefig(h_spec,FileNameWaveInd);
                end
            end
        end
    end
end
