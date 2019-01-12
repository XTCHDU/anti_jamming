%This simulation is used for radar signal (MP, BPSK, LFM and NLFM) detection and parameter estimation
close all
clear all
warning off
fADC = 10e6; %sampling rate
fCar = 4e6; %carrier frequency
Ttotal = 7.5e-6; %total time duration
T = 5.5e-6; %signal duration
A = 1; %signal amplitude
num_jam = 2;
t = 0:1/fADC:T - 1/fADC; %signal time

tTotal = 0:1/fADC:Ttotal - 1/fADC; %time axis

Nspace = round((Ttotal - T)*fADC/2); %points of time without signal

%Generate MP signal
sigMP = A*exp(j*2*pi*fCar*t);    %Generate the MP signal
sigMP = [zeros(1,Nspace),sigMP, zeros(1,Nspace)];

%Generate BPSK signal (11 bit Barker code)
Len = length(t);
Ncode = 11;
for kIndex = 1:Len
    if kIndex<Len/Ncode
        sBPSK(kIndex)=1;
    elseif kIndex<Len/Ncode*2
        sBPSK(kIndex)=1;
    elseif kIndex<Len/Ncode*3
        sBPSK(kIndex)=1;
    elseif kIndex<Len/Ncode*4
        sBPSK(kIndex)=-1;
    elseif kIndex<Len/Ncode*5
        sBPSK(kIndex)=-1;
    elseif kIndex<Len/Ncode*6
        sBPSK(kIndex)=-1;
    elseif kIndex<Len/Ncode*7
        sBPSK(kIndex)=1;
    elseif kIndex<Len/Ncode*8
        sBPSK(kIndex)=-1;
    elseif kIndex<Len/Ncode*9
        sBPSK(kIndex)=-1;
    elseif kIndex<Len/Ncode*10
        sBPSK(kIndex)=1;
    else
        sBPSK(kIndex)=-1;
    end
end
sMP = A*exp(j*2*pi*fCar*t);    %Generate the MP signal
sigBPSK = sMP.*sBPSK;
sigBPSK = [zeros(1,Nspace),sigBPSK, zeros(1,Nspace)];

%Generate LFM signal
B = 4e6;
chirpRate =B/Ttotal;
sigLFM = A*exp(j*2*pi*fCar*t + j*pi*chirpRate*t.^2);
sigLFM = [zeros(1,Nspace),sigLFM, zeros(1,Nspace)];

%Generate NLFM signal
a1 = 1e11;
a2 = 5e16;
sigNLFM = A*exp(j*2*pi*fCar*t + j*pi*a1*t.^2 + j*pi*a2*t.^3);
sigNLFM = [zeros(1,Nspace),sigNLFM, zeros(1,Nspace)];

%Input signal
sig = sigLFM;%change signal type for different input
figure;
plot(tTotal*1e6,sig)
xlabel('t/us')
ylabel('Amplitude')

sig = awgn(sig,1000,'measured');
figure;
plot(tTotal*1e6,sig)
xlabel('t/us')
ylabel('Amplitude')
title('无干扰回波时域波形')
%%%%%%%%%%%加干扰
res = zeros(1,1000);
count = 0;
sig_sum = zeros(120,length(sig)*5);
norm_f = zeros(num_jam,length(sig));
Amp = zeros(num_jam,1);
for JNR = 30
    receiveSignal = sig;
    jam_total = zeros(1,length(sig));
    for jam_idx = 1:num_jam
        jam_type = 2;%1 噪声调幅干扰；2 噪声调频干扰
        jam_freq = 12e6;
        w_signal = sum(abs(sig).^2)/length(sig);
        w_jamming = sqrt(10^(JNR/10)*w_signal)*jam_idx;
        B = 3e6;
        [Amp(jam_idx),norm_f(jam_idx,:),jam]= jamming(length(tTotal),jam_type,fADC,w_jamming,B,T,jam_freq);%生成干扰
        receiveSignal = receiveSignal + jam;
        jam_total = jam_total + jam;
        figure;
        plot(tTotal*1e6,jam)
        xlabel('t/us')
        ylabel('Amplitude')
        title(sprintf('干扰时域波形%d',jam_idx))
    end
    jam = jam_total;
    figure;
    plot(tTotal*1e6,jam_total)
    xlabel('t/us')
    ylabel('Amplitude')
    title('干扰时域波形')
    figure;
    plot(tTotal*1e6,receiveSignal)
    xlabel('t/us')
    ylabel('Amplitude')
    title('干扰+回波时域波形')
    
    
    %%%%%%%%%%%%
    fAxis = linspace(-fADC/2,fADC/2,length(tTotal));
    figure;
    plot(fAxis,fftshift(abs(fft(sig))))%
    xlabel('频率/hz')
    ylabel('幅度')
    title('无干扰回波频谱图')
    figure;
    plot(fAxis,fftshift(abs(fft(jam))))
    xlabel('频率/hz')
    ylabel('幅度')
    title('干扰频谱')
    figure;
    plot(fAxis,fftshift(abs(fft(receiveSignal))))%
    xlabel('频率/hz')
    ylabel('幅度')
    title('干扰+回波频谱图')
    
    %%%%%%%%%%%%
    re_sig = dataProcess(receiveSignal,length(tTotal),jam_type,fADC,norm_f,Amp);
    figure;
    plot(real(re_sig));
    xlabel('t/us')
    ylabel('Amplitude')
    title('复原信号时域图')
    figure;
    plot(fAxis,fftshift(abs(fft(re_sig))))
    xlabel('频率/hz')
    ylabel('幅度')
    title('复原信号频谱图')
    w_re_j = sum(abs(sig-re_sig).^2)/length(re_sig);
end





