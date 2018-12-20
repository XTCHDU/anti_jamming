function valueCodeRate = CodeRateEst2(fSampling,vSignalInput)
%ChirpRateEstimation--Chirp rate estimation of a certain LFM ignal
%INPUTs
%fSampling
%   sampling frequency
%vSignalInput
%   vSignalInput means a vector of signal
%valueCodeRate
%   value of code rate

%Pre-processing
numSig = length(vSignalInput);  %number of sampling points
numDelay = 10;   %number of delay points
numEff = numSig - numDelay; %Effective number of points in operation


%Realization of the delay-conjugate-multiplication algorithm
%Two parts of multiplication operation
sigPrevious = vSignalInput(1:numEff);  %The previous part
sigDelayed = vSignalInput(numDelay+1:numSig);  %The delayed part

% conjugate-multiplication operation
sigConjMult = conj(sigDelayed).*sigPrevious;    %Result of Conj-Mult

spectrumSig = abs(fft(sigConjMult));
spectrumSig(1:10) = 0;

locMax = 0;
for kIndex = 2:length(spectrumSig)/2
    if spectrumSig(kIndex) > spectrumSig(kIndex - 1)&&spectrumSig(kIndex) > spectrumSig(kIndex + 1)
        locMax = kIndex;
        break;
    end
end
vFreq = linspace(0,fSampling,numEff);
% valueCodeRate = (locMax-1)/numEff*fSampling;%
valueCodeRate = vFreq(locMax);
figure;
plot(spectrumSig/max(spectrumSig))%vFreq,
xlabel('frequency/Hz')
ylabel('Normalized Amplitude')
figure;
plot(vFreq,spectrumSig/max(spectrumSig))%
