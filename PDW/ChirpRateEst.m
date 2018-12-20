function valueChirpRate = ChirpRateEst(fSampling,vSignalInput)
%ChirpRateEstimation--Chirp rate estimation of a certain LFM signal
%INPUTs
%fSampling
%   sampling frequency
%vSignalInput
%   vSignalInput means a vector of signal
%OUTPUT
%valueChirpRate
%   value of chirp rate

%Pre-processing
numSig = length(vSignalInput);  %number of sampling points
numDelay = floor(0.4*numSig);   %number of delay points
numEff = numSig - numDelay; %Effective number of points in operation

%Realization of the delay-conjugate-multiplication algorithm
%Two parts of multiplication operation
sigPrevious = vSignalInput(1:numEff);  %The previous part
sigDelayed = vSignalInput(numDelay+1:numSig);  %The delayed part

% conjuate-multiplication operation to convert signal into monopulse
sigConjMult = sigDelayed.*conj(sigPrevious);    %Result of Conj-Mult
valueFreq = FreqEstFFT(fSampling,sigConjMult);   %Frequency estimation
valueChirpRate = valueFreq/numDelay*fSampling;   %The final result