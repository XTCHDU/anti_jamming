%This simulation is used for radar signal (MP, BPSK, LFM and NLFM) detection and parameter estimation
close all
clear all

fADC = 100e6; %sampling rate
fCar = 5e6; %carrier frequency
Ttotal = 7.5e-6; %total time duration
T = 5.5e-6; %signal duration
A = 1; %signal amplitude

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
chirpRate = 1e12;
sigLFM = A*exp(j*2*pi*fCar*t + j*pi*chirpRate*t.^2);
sigLFM = [zeros(1,Nspace),sigLFM, zeros(1,Nspace)];

%Generate NLFM signal
a1 = 1e12;
a2 = 5e17;
sigNLFM = A*exp(j*2*pi*fCar*t + j*pi*a1*t.^2 + j*pi*a2*t.^3);
sigNLFM = [zeros(1,Nspace),sigNLFM, zeros(1,Nspace)];

%Input signal
sig = sigLFM;%change signal type for different input
figure;
plot(tTotal*1e6,sig)
xlabel('t/us')
ylabel('Amplitude')


sig = awgn(sig,10,'measured');
figure;
plot(tTotal*1e6,sig)
xlabel('t/us')
ylabel('Amplitude')
title('�޸��Żز�ʱ����')
%%%%%%%%%%%�Ӹ���
res = zeros(1,1000);
count = 0;
for JNR = -10
     count = count + 1;
     x_lab(count) = JNR;

    for ROUND = 1:1
        close all;
        %JNR = 10;%����JNR
        jam_type = 1;%1 �����������ţ�2 ������Ƶ����
        jam_freq = 13e6; %type==1ʱ�����ŵ�Ƶ������ ע�����ź��м�����ź���5e6-10e6֮�䣬�ź������12e6����
        w_signal = sum(abs(sig).^2)/length(sig);
        w_jamming = sqrt(10^(JNR/10)*w_signal);
        B = 3e6;
        [mag1,jam]= jamming(length(tTotal),jam_type,fADC,w_jamming,B,T,jam_freq);%���ɸ���
        mag1 = mod(mag1,2*pi);
        receiveSignal = sig + jam;
        figure;
        plot(tTotal*1e6,jam)
        xlabel('t/us')
        ylabel('Amplitude')
        title('����ʱ����')
        figure;
        plot(tTotal*1e6,receiveSignal)
        xlabel('t/us')
        ylabel('Amplitude')
        title('����+�ز�ʱ����')
        
        
        %%%%%%%%%%%%
        fAxis = linspace(-fADC/2,fADC/2,length(tTotal));
        figure;
        plot(fAxis,fftshift(abs(fft(sig))))%
        xlabel('Ƶ��/hz')
        ylabel('����')
        title('�޸��Żز�Ƶ��ͼ')
        figure;
        plot(fAxis,fftshift(abs(fft(jam))))
        xlabel('Ƶ��/hz')
        ylabel('����')
        title('����Ƶ��')
        figure;
        plot(fAxis,fftshift(abs(fft(receiveSignal))))%
        xlabel('Ƶ��/hz')
        ylabel('����')
        title('����+�ز�Ƶ��ͼ')
        
        %%%%%%%%%%%%
        re_sig = dataProcess(receiveSignal,length(tTotal),jam_type,fADC);
        figure;
        plot(real(re_sig));
        xlabel('t/us')
        ylabel('Amplitude')
        title('��ԭ�ź�ʱ��ͼ')
        figure;
        plot(fAxis,fftshift(abs(fft(re_sig))))
        xlabel('Ƶ��/hz')
        ylabel('����')
        title('��ԭ�ź�Ƶ��ͼ')
        w_re_j = sum(abs(sig-re_sig).^2)/length(re_sig);
        res(count) = res(count)+10*log10(w_re_j/w_signal);
        
        
    end
end
 res = res./ROUND;
 figure;
 plot(x_lab,res(1:length(x_lab)))
xlabel('����JNR��dB��')
ylabel('���ƺ�JNR��dB��')
%%%%%%%%%%%


%Signal detection
% Flag = Detection(sig);
%
% %TOA and PW estimation
% if Flag == 1
%     [numTOA,numTOE,sigData] = TOA(sig);
% else
%     numTOA = 0;
%     numTOE = 0;
% end
% estTOA = tTotal(numTOA);
% estTOE = tTotal(numTOE);
% estPW = estTOE - estTOA + 1/fADC;
%
% tTemp = 0:1/fADC:(length(sigData)-1)/fADC;
%
% figure;
% plot(tTemp*1e6,sigData)
% xlabel('t/us')
% ylabel('Amplitude')
% %Signal modulation identification
% sigType = ModulationIdentification(sigData);
%
% %Signal parameter estimatiom
% if sigType == 1
%     fCarEst = FreqEstFFT(fADC,sigData);
% elseif sigType == 2
%     codeRateEst = CodeRateEst(fADC,sigData);
%     sigDataSq = sigData.^2;
%     fCarEst = FreqEstFFT(fADC,sigDataSq)/2;
% elseif sigType == 3
%     chirpRateEst = ChirpRateEst(fADC,sigData);
%     sigDemod = sigData.*exp(-j*pi*chirpRateEst*tTemp.^2);
%     fCarEst = FreqEstFFT(fADC,sigDemod);
% elseif sigType == 4
%     a2Est = NLFMParaEst(fADC,sigData);
%     sigDemod = sigData.*exp(-j*pi*a2Est*tTemp.^3);
%     a1Est =  ChirpRateEst(fADC,sigDemod);
%     sigDemod2 = sigData.*exp(-j*pi*a1Est*tTemp.^2 - j*pi*a2Est*tTemp.^3);
%     fCarEst = FreqEstFFT(fADC,sigDemod2);
% else
%     fCarEst = 0;
% end
%
%
%
