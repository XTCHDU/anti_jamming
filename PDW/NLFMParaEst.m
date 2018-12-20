function a2 = NLFMParaEst(fSampling,vSignalInput)
%NLFMParaEst --- Intrapulse parameter estimation for NLFM signal
%INPUTs
%fSampling
%   sampling frequency
%vSig
%   signal input
%OUTPUTs:
%values of estimated NLFM signal parameters

%Pre-processing
numSig = length(vSignalInput);  %number of sampling points
numDelay = floor(numSig/3);   %number of delay points
numEff = numSig - numDelay; %Effective number of points in operation

%Realization of the delay-conjugate-multiplication algorithm
%Two parts of multiplication operation
sigPrevious = vSignalInput(1:numEff);  %The previous part
sigDelayed = vSignalInput(numDelay + 1:numSig);  %The delayed part

% DPT
sigConjMult = sigDelayed.*conj(sigPrevious);    %Result of Conj-Mult
numPoints = length(sigConjMult);
sigConjMultFFT = fft(sigConjMult);

[maxFFT,maxIndex] = max(abs(sigConjMultFFT));
valueCoef = 0.8;  %sqrt(0.5) sqrt(1/8) bandwidth
valueThre = valueCoef*maxFFT; 

for kLine = 1:maxIndex
    if abs(sigConjMultFFT(kLine)) > valueThre
        indexConLeft = kLine; %The left bound index
        break;
    end
end
%   find the right bound of 3dB bandwidth
for kLine = numPoints:-1:maxIndex
    if abs(sigConjMultFFT(kLine)) > valueThre
        indexConRight = kLine; %The right bound index
        break;
    end
end
sigConjMultFFT(1:indexConLeft - 20) = 0;
sigConjMultFFT(indexConRight + 20:end) = 0;
sigConjMult = ifft(sigConjMultFFT);


figure;
plot(abs(sigConjMultFFT))



%DPT again
numSigAg = length(sigConjMult);  %number of sampling points
numDelayAg = floor(numSig/3);   %number of delay points
numEffAg = numSigAg - numDelayAg; %Effective number of points in operation
sigPreviousAg = sigConjMult(1:numEffAg);  %The previous part
sigDelayedAg = sigConjMult(numDelayAg + 1:numSigAg);  %The delayed part
sigConjMultAg = sigDelayedAg.*conj(sigPreviousAg);    %Result of Conj-Mult

valueFreq = FreqEstFFTMore(fSampling,sigConjMultAg);   %Frequency estimation
a2 = valueFreq/(3*numDelay*numDelayAg)*(fSampling^2);   %The final result


