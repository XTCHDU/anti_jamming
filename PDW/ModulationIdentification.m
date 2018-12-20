function [ModType] = ModulationIdentification(vSignalInput)
%ModulationIdentification---Identify the modulation type of input signal 
%INPUTs
%fSampling
%   sampling frequency. unit: Hz
%vSignalInput
%   vSignalInput means a vector of signal 

%OUTPUTs

%ModType
%   modulation type:
%       1----Monopulse
%       2----BPSK
%       3----LFM
%       4----NLFM
%       5----fail to identify

%   extract the signal pulse
sigPulse = vSignalInput;
numPoints = length(sigPulse);

if numPoints>= 512
    numFFT = numPoints;     %FFT points
else
    numFFT = 512;
end


%   estimate the signal bandwidth
sigPulseFFT = abs(fft(sigPulse,numFFT));%
numPoints = length(sigPulseFFT);
[maxFFT,maxIndex] = max(sigPulseFFT);
%   find the left bound of 3dB bandwidth

valueCoef = sqrt(1/4);  %sqrt(1/8) bandwidth
valueThre = valueCoef*maxFFT; 

for kLine = 1:maxIndex
    if ((sigPulseFFT(kLine) > valueThre))
        indexBoundLeft = kLine; %The left bound index
        break;
    end
end
%   find the right bound of 3dB bandwidth
for kLine = numPoints:-1:maxIndex
    if ((sigPulseFFT(kLine) > valueThre))
        indexBoundRight = kLine; %The right bound index
        break;
    end
end
valueBandWidth = indexBoundRight - indexBoundLeft;


sigPulseSquared = sigPulse.^2;
%   estimate the 3dB signal bandwidth
sigPulseSquaredFFT = abs(fft(sigPulseSquared,numFFT));%
[maxSqFFT,maxSqIndex] = max(sigPulseSquaredFFT);
%   find the left bound of 3dB bandwidth
valueSqThre = valueCoef*maxSqFFT; 
for kLine = 1:maxSqIndex
    if ((sigPulseSquaredFFT(kLine) > valueSqThre))
        indexSqBoundLeft = kLine; %The left bound index
        break;
    end
end
%   find the right bound of 3dB bandwidth
for kLine = numPoints:-1:maxSqIndex
    if ((sigPulseSquaredFFT(kLine) > valueSqThre))
        indexSqBoundRight = kLine; %The right bound index
        break;
    end
end

valueSqBandWidth = indexSqBoundRight - indexSqBoundLeft; %original method


if abs(valueBandWidth - valueSqBandWidth)<2
    ModType = 1;    % Identification of monopulse
elseif valueBandWidth - valueSqBandWidth >= 3%BPSKºÍNLFM»á»ì
    ModType = 2;     % Identification of BPSK
elseif valueBandWidth - valueSqBandWidth <= -5
    ModType = 3;    %   Identification of LFM
else
    ModType = 5;    %   fail to Identify 
end

if ModType == 3
    numDelay = 20;
    if numPoints - numDelay<1
        numDelay = 1;
    end
        
    numEff = numPoints - numDelay;
    

    sigPrevious = vSignalInput(1:numEff);  %The previous part
    sigDelayed = vSignalInput(numDelay+1:numPoints);  %The delayed part

    % conjuate-multiplication operation to convert signal into monopulse
    sigConjMult = sigDelayed.*conj(sigPrevious);    %Result of Conj-Mult

    sigConFFT = fftshift(abs(fft(sigConjMult)));
    [maxConFFT,maxConIndex] = max(sigConFFT);
    valueCoef1 = sqrt(0.7);
    valueConThre = valueCoef1*maxConFFT;
    
    for kLine = 1:maxConIndex
        if ((sigConFFT(kLine) > valueConThre))
            indexConLeft = kLine; %The left bound index
            break;
        end
    end
    %   find the right bound of 3dB bandwidth
    for kLine = length(sigConFFT):-1:maxConIndex
        if ((sigConFFT(kLine) > valueConThre))
            indexConRight = kLine; %The right bound index
            break;
        end
    end
    valueConBandWidth = indexConRight - indexConLeft;
    
    %ÂË²¨
    sigConFFTTemp = fftshift(fft(sigConjMult));
    sigConFFTTemp(1:indexConLeft - 20) = 0;
    sigConFFTTemp(indexConRight + 20:end) = 0;
    sigConjMultTemp = ifft(fftshift(sigConFFTTemp));
    sigConjMult = sigConjMultTemp;


    sigConSquared = sigConjMult.^2;
    %   estimate the 3dB signal bandwidth
    sigConSquaredFFT = fftshift(abs(fft(sigConSquared,numFFT)));%
    [maxConSqFFT,maxConSqIndex] = max(sigConSquaredFFT);
    %   find the left bound of 3dB bandwidth
    valueCoef2 = sqrt(0.5);
    valueConSqThre = valueCoef2*maxConSqFFT; 
    for kLine = 1:maxConSqIndex
        if ((sigConSquaredFFT(kLine) > valueConSqThre))
            indexSqConLeft = kLine; %The left bound index
            break;
        end
    end
    %   find the right bound of 3dB bandwidth
    for kLine = numPoints:-1:maxConSqIndex
        if ((sigConSquaredFFT(kLine) > valueConSqThre))
            indexSqConRight = kLine; %The right bound index
            break;
        end
    end
    
    valueConSqBandWidth = indexSqConRight - indexSqConLeft; %original method
    if abs(valueConBandWidth - valueConSqBandWidth)<1
        ModType = 3;    % Identification of LFM
    else
        ModType = 4;    % Identification of NLFM
    end
end
    
    