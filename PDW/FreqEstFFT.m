function valueFreq = FreqEstFFT(fSampling,vSig)
%FreqEstiamtionFFT --- Intrapulse parameter estimation
%INPUTs
%fSampling
%   sampling frequency
%vSig
%   signal input
%OUTPUTs:
%value of estimated frequency

numPoints = 1024;%length(vSig)
vFreq = linspace(-fSampling/2,fSampling/2,numPoints);
vSpectrumSig = fftshift(abs(fft(vSig,numPoints)));
[~,maxIndex] = max(vSpectrumSig);
valueFreq = vFreq(maxIndex);