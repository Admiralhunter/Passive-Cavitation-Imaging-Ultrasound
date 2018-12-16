function [ AvePowerSpectraAllFrames_g, FocalPowerSpectra,FocalWaveform,Window,FocalWaveWin] = Analyze_FrameV3_20170707(x,z,RFFTgpu,NFrameStart,NFramesAnalyze,tofgpu,fgpu,f)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Ccode is optimized and allows any sized array of Frequency, X range, and Z
%range that is smaller than VRAM of the GPU. 
%In particular windowing is applied at the focal point to supress side
%lobes. Windowing is incorrect for any other pixel location
%This function loops through a set of frames from a data file.
%
%NNOTE: Higher Z resolutions deplete vRAM alot faster than X resolution.
%This is because the for loop runs through X dimensions and so just runs
%more times, not increasing the storage much. The Z dimension will increase
%the matrices for calculations and thus requires alot more memory.
%
%NNOTE: GTX 1000 series is pascal series and should be run using MATLAB
%2017a. Previous versions do utilize the gpu, but not as efficiently. Also
%GTX series is only efficient for singles data type not doubles.
%DO NOT USE DOUBLES FOR COMPUTATION, IS ROUGHLY 1/32 THE SPEED OF SINGLES.
g = gpuDevice;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%preallocation of a few arrays
AvePowerSpectraAllFramesgpu = gpuArray(single(zeros(length(x),length(z),length(fgpu))));
%RFFT_shifted = gpuArray(single(zeros(size(RFFTgpu,1),1,'single')));
FocalPowerSpectra = single(zeros(size(f)));
Fs = 2*max(f); dt = 1/Fs; Time = (1:length(f))*dt;
spacing = single(size(tofgpu,1)*size(tofgpu,2));

fgpu= transpose(fgpu);
%if array is small enough, compute with fgpu outside of loop (vastly
%improves speed, but requires more memory. Value is currently set to 3GB's,
%if larger VRAM is used this value can be increased
if numel(tofgpu)*numel(fgpu)*4 < 3e+09
    large =0;
    y = exp(fgpu*-2i*pi*-1*reshape(tofgpu,1,size(tofgpu,1)*size(tofgpu,2)*size(tofgpu,3),1));
else
    large = 1;
    y = -2i*pi*-1*reshape(tofgpu,1,size(tofgpu,1)*size(tofgpu,2)*size(tofgpu,3),1);
end
totaltime = 0;
xmm = x;
for n = NFrameStart:NFramesAnalyze

    tic
    SingleFrameOfData = squeeze(RFFTgpu(:,:,n));
    
   
    
    %h is the step index for individual sensors for each pixel (128
    %sensors for each pixel.
    h =1:spacing:length(y);

    %Loops through h to get single column associated with pixel from h and
    %stored in a vector
    for i = 2:length(h)+1
       holder = (h(i-1):h(i-1) + length(z)-1);
        xx((i-2)*length(z) +1:(i-1)*length(z)) = holder;
    end
    sx = single(zeros(1,length(xx)));
    
    
    %Rearranges matrix xx to get individual pixels from each 128 sensors to
    %be next to eachother in the row
    for i = 1:length(holder)    
        sx(1,(i-1)*length(h) +1 :(i)*length(h)) = xx(1, i:length(z):length(xx));
    end
    x =single(sx -size(tofgpu,1)) ;
    %Loops through all 128 sensors for each column
    for i = 1:size(tofgpu,2)
        x =x +size(tofgpu,1);

        %shifts for entire column of data and then reshapes into a new
        %matrix that seperates each pixel column into a seperate sheet.
        %Switch statement determines if enough VRAM exists to hold matrix
        %of f*tofgpu. If not this has to be calculated on smaller arrays,
        %which is less efficient
        switch large
            case 0
                holder = y(:,x);
            case 1
                holder = exp(fgpu*y(:,x));
        end

        holdershape = reshape(holder,length(fgpu),128,[]);
        
        %%%%%%%
        %THIS IS IMPLICIT EXPANSION
        %MATLAB 2016b AND HIGHER IS REQUIRED
        shifted = SingleFrameOfData.*holdershape;
        %shifted = SingleFrameOfData.*reshape(exp(fgpu*y(:,x)),length(fgpu),128,[]);
        %%%%%%%
        shifted3d = reshape(shifted,length(fgpu),128,[]);

      
        %Each column is summed via the row which creates z amount of pixels with f
        %amount of frequencies.
        
        RFFT_shifted3d = sum(shifted3d,2);
        RFFT_shifted = permute(reshape(RFFT_shifted3d,length(fgpu),size(tofgpu,1),1),[2 1]);
        
        
        WindowLength = 240; % must be an even value
        WaveformIndicies = 1700:2200;   % Set these values manually to a region that includes the main pulse but not any reverberations
        

        %Location for window
        Lateral = .3448;
        Range = 29.08;
        
        %index range
        %Lower = 1830;
        %Upper = Lower + 240;

        
        %Computes window
       if i==closest(xmm,Lateral)
            RFFT_shifted_g = gather(RFFT_shifted);
            FocalWaveform = fftshift(ifft(ifftshift(RFFT_shifted_g(closest(z,Range),:))));
            % Find center point of waveform
            MidPointInd = closest(cumsum(abs(FocalWaveform(WaveformIndicies)))./max(cumsum(abs(FocalWaveform(WaveformIndicies)))),0.5)+WaveformIndicies(1);
            Lower = MidPointInd - WindowLength/2;
            Upper = MidPointInd + WindowLength/2-1;           
            Window = [zeros(Lower,1); chebwin(Upper-Lower,300); zeros(length(f)-Upper,1)]';
            FocalWaveWin = FocalWaveform.*Window/mean(chebwin(Upper-Lower,300).^2);
            FocalPowerSpectra = FocalPowerSpectra + abs(fftshift(fft(FocalWaveWin))).^2;
       end
        

        %Adds the shifted column into a matrix for final processing
        container = abs(RFFT_shifted).^2;
        container = reshape(container,1,[],length(fgpu));
        AvePowerSpectraAllFramesgpu(i,:,:) = AvePowerSpectraAllFramesgpu(i,:,:) + container;
           
            
       
        

        
        

    end
    
    
    wait(g)
    timer = toc;
    totaltime = totaltime + timer;
    Hz = 1/timer;
    fprintf('Frame %i was completed in %0.3fs which corresponds to a frame rate of %0.1fHz\n',n,timer,Hz)
end
%averages power spectra and then gathers data from GPU
AvePowerSpectraAllFrames = AvePowerSpectraAllFramesgpu./n;
FocalPowerSpectra = FocalPowerSpectra./n;
AvePowerSpectraAllFrames_g  = gather(AvePowerSpectraAllFrames);
fprintf('%0.3fs \n',totaltime)
g.AvailableMemory
end

