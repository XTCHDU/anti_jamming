clear all;close all;

%����LFM�ź�
JSR = 60; %���ű�
j = sqrt(-1); 
fs = 20e6; %����Ƶ��
B = 2e6; %����
f0 = 4e6; %����Ƶ��
T = 50e-6;
k = B/T;
N = T*fs; %��������
t = linspace(0,T,N);
lfm = exp(j*2*pi*f0*t+j*pi*k*t.^2); %LFM�ź�

figure(1);

subplot(2,1,1);
plot(t,lfm); %ʵ��
title('lfm�ź�ʱ����');
lfm_fft=fft(lfm);
f=(0:N-1)*fs/N;
subplot(2,1,2);
plot(f,abs(lfm_fft));
title('lfm�ź�Ƶ����');

%���������ź�
u0=1;
noise_power_dB = 10; 
fj=4e6; %������Ƶ
Tr=1000/fs;     
t1=0:1/fs:Tr-1/fs; 
N=length(t1); 
u=wgn(1,N,noise_power_dB); %��˹������
df1=fs/N; %���
n=0:N/2;
wp=1e6;
ws=2e6;  
rp=1; 
rs=60;  
[n1,wn1]=buttord(wp/(fs/2),ws/(fs/2),rp,rs);  
[b,a]=butter(n1,wn1);  
u1=filter(b,a,u);   %�õ���������
p1=(1/N)*sum(u1.^2);   %����������ƽ������
y=(u0+u1).*exp(j*2*pi*fj*t1)
%�������������ź�  �ѵ���

figure(2);

subplot(2,2,1);
plot(t1,u1);
title('�������Ʋ���');

subplot(2,2,2);
u1_fft=fft(u1);
plot(10*log10(abs(u1_fft(n+1)*2/N)));
title('��������������');
rand('state',0);

subplot(2,2,3);
plot(t1,real(y));%ʵ��
title('������������ʱ����');

subplot(2,2,4);
y_fft=fft(y,N);
mag=abs(y_fft);
plot(f*1e-6,mag);
xlabel('MHZ');
title('������������Ƶ����');

pj=(1/N)*sum(y.^2); %������������ƽ������
plfm=(1/N)*sum(lfm.^2); %LFM�ź�ƽ������

temp = 10^(JSR/10)*plfm;
Jamming = sqrt(temp/pj) .* y;%����JSR�����µĸ����ź�
% Jamming = (temp + u1).*exp(j*2*pi*fj*t1)

pJamming = (1/N)*sum(Jamming.^2);
JSR1= 10*log10(pJamming/plfm);

xn=lfm + Jamming;%�״���ջ��ź�

figure(3);

xn_fft=fft(xn);
mag3=abs(xn_fft);
plot(f,mag3);
xlabel('MHZ');
title('�������ź�Ƶ����');
Jamming_fft=abs(fft(Jamming));
[M I] = max(Jamming_fft);
estmated_f = (I-1)*fs/N;%�����ź���Ƶ����ֵ

%��� ��Ƶ����ֵestmated_f
y3 = (lfm+Jamming) .* exp(-j*2*pi*estmated_f*t1);  %y(t)
y4 = lfm.* exp(-j*2*pi*estmated_f*t1);  %y1(t)
y5 = Jamming.* exp(-j*2*pi*estmated_f*t1);%y2(t)
y6 = y4 + y5;
F_y3 = abs(fft(y3)); %xnƵ��ͼ Y(t)
F_y4 = abs(fft(y4)); %lfmƵ��ͼ Y1(t)
F_y5 = abs(fft(y5)); %JammingƵ��ͼ Y2(t)
F_y6 = F_y4 + F_y5; %Y1(t)+Y2(t)

figure(4);
plot(f*1e-6,F_y6);
title('Y(f)Ƶ��ͼ');

%Ƶ�����

Y_L=F_y6(2:501);
Y_R=F_y6(501:1000);
 
% Y_L=F_y3(2:N/2);
% Y_R=F_y3(N/2+1:N-1);

Y_R1=conj(fliplr(Y_R));
% Y_R1=Y_R(end:-1:1);

figure(5)
subplot(3,1,1);
plot(abs(Y_L));
xlabel('MHZ');
title('Y_L(f)Ƶ��ͼ');
subplot(3,1,2);
plot(abs(Y_R1));
xlabel('MHZ');
title('Y_R1(f)����Ƶ��ͼ');

% Y_L=Y_L(2:end);
% Y_L=[Y_L,zeros(1,1)];

subplot(3,1,3);
plot(abs(Y_L-Y_R1));
title('���Ŷ������Ƶ��ͼ');

figure(6)
subplot(2,1,1);
plot(abs(lfm_fft));
subplot(2,1,2);
plot(abs(Y_L-Y_R1));
axis([0 1000 0 120]);