%fs = 4*35e6;
function [ lfmwav , E2] = LFM(fs) 
sLFM = phased.LinearFMWaveform('SampleRate',fs,...
    'SweepBandwidth',2e6,...
    'PulseWidth',1e-3,'PRF',1e3);

lfmwav =step( sLFM);
nsamp = size(lfmwav,1);
t = [0:(nsamp-1)]/fs;

% figure(2);
% plot(t*1000,real(lfmwav));
% xlabel('Time (millisec)')
% ylabel('Amplitude')
grid
end
