function valueFreq = FreqEstFFTMore(fSampling,vSig)
%FreqEstiamtionFFT --- frequency estimation with high accuracy
%INPUTs
%fSampling
%   sampling frequency
%vSig
%   signal input
%OUTPUTs:
%value of estimated frequency

numPoints = 2048;
vFreq = linspace(-fSampling/2,fSampling/2,numPoints);
vSpectrumSig = fftshift(abs(fft(vSig,numPoints)));
[~,maxIndex] = max(vSpectrumSig);
valueFreq = vFreq(maxIndex);