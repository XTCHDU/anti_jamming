clear all;close all;
%生成LFM信号
JSR = 20;
j=sqrt(-1);
fs=20e6;   %采样频率
B=2e6;  %带宽
f0=4e6;   %中心载频
T=50e-6;  %脉冲宽度
k=B/T;  
N=T*fs;   %采样点
t=linspace(0,T,N);
lfm=exp(j*2*pi*f0*t+j*pi*k*t.^2);
ht=exp(-j*2*pi*f0*t-j*pi*k*t.^2);%匹配滤波器
figure(1);
subplot(2,1,1);
plot(t,real(lfm));
title('lfm信号时域波形');
lfm_fft=fft(lfm);
% f=linspace(-fs/2,fs/2,N);
f=(0:N-1)*fs/N;
subplot(2,1,2)
plot(f,abs(lfm_fft));
title('lfm信号频域波形');

%噪声调频
fj=4.5e6;   %干扰载频
bj=10e6;                   
Tr=1000/fs; 
t1=0:1/fs:Tr-1/fs;
N=length(t1);
 uj=60;
 mf=0.6; 
 wpp=6;     
u=wgn(1,N,wpp);               
wp=10e6; %调制噪声带宽                       
ws=13e6; 
rp=1; rs=60; 
[Nn,wn]=buttord(wp/(30e6/2),ws/(30e6/2),rp,rs); 
[b,a]=butter(Nn,wn);           
u1=filter(b,a,u);              
ss(1)=0;                       
for i=1:N-1                    
ss(i+1)=u1(i)+ss(i);
end
i=1:N-1;  %ss=cumsum([0,u1(i)]);
ss=ss*Tr/N;                        
y0=uj*exp(j*(2*pi*fj*t1+2*pi*mf*bj*ss+2*pi*rand(1,1)));  

 %计算干信比
 pj=(1/N)*sum(y0.^2);   %噪声调频干扰的平均功率=u0^2/2+p1^2/2
plfm=(1/N)*sum(lfm.^2);   %LFM信号平均功率
 temp = 10^(JSR/10)*plfm;
 y0 = sqrt(temp/pj) .* y0;
% pJamming = (1/N)*sum(Jamming.^2);
% JSR1= 10*log10(pJamming/plfm);

 figure(7)
 subplot(2,1,1)
 y1=fft(y0,N);
 mag2=abs(y1);
 plot(f,mag2);
 xlabel('HZ');
 title('噪声调频干扰频域波形');
 subplot(2,1,2)
 plot(f,mag2 + abs(lfm_fft));
 xlabel('HZ');
 title('噪声调频干扰加雷达信号频域波形');


%fmwav0 雷达信号  fmwav1 常规距离假目标欺骗干扰  noiseFM()噪声调频干扰
J1 = y0;
J0 =  lfm ;
s = J0 ./ J1 ;
x = J0 + J1 ;
figure(2)
subplot(2,2,1)
plot(real(x));
s1 = log(x);
%for i=1:size(x,2)
 %   s1(i)= log(J1(i)) + log(a(i) + s(i));
%end
Rx = real(s1) ;
Ix = imag(s1) ;
A = exp(mean(Rx));
s2 = x - A  * exp(j .* Ix) ;
figure(2)
subplot(2,1,1)
plot(real(y0));
title('调频噪声时域');
subplot(2,1,2)
plot(real(s2));
title('干扰抑制后信号时域');

s21 = fliplr(s2);
s22 = conj(s21);
n = conv(lfm , s22);   %估计脉压
fmwav1=fliplr(lfm);
fmwav2=conj(fmwav1);
m=conv(lfm ,fmwav2);  % 真实脉压
x1 = fliplr(x);
x2 = conj(x1);
n1 = conv(lfm , x2);   %接受信号脉压
figure(3)
subplot(2,1,1)
plot(mapminmax(real(m)));
title('雷达信号脉压');
subplot(2,1,2)
plot(mapminmax(real(n1)));
title('雷达信号加干扰的脉压');
figure(4)
subplot(2,1,1)
plot(mapminmax(real(m)));
title('雷达信号脉压');
subplot(2,1,2)
plot(mapminmax(real(n)));
title('抑制干扰后的脉压');




% %调幅干扰
% u0=1;
% noise_power_dB = 10; 
% fj=4.5e6;%干扰载频
% %fs=4*fj; %采样频率
% Tr=1000/fs;     
% t1=0:1/fs:Tr-1/fs; 
% N=length(t1); 
% u=wgn(1,N,noise_power_dB); %噪声
% df1=fs/N; %间隔
% n=0:N/2;
% 
% wp=3e6; %调制噪声带宽10M
% ws=5e6;  
% rp=1; 
% rs=60;  
% [n1,wn1]=buttord(wp/(fs/2),ws/(fs/2),rp,rs);  
% [b,a]=butter(n1,wn1);  
% u1=filter(b,a,u);   %得到带限噪声
% p1=(1/N)*sum(u1.^2);   %带限噪声的平均功率
% y=(u0+u1).*exp(j*2*pi*fj*t1+2*pi*rand(1,1));%噪声调幅干扰信号  已调波
% 
% % figure(2);
% % subplot(2,2,1);
% % plot(t1,u1);
% % title('噪声调制波形');
% % %axis([0,0.05e-4,-2,2]);
% % 
% % subplot(2,2,2);
% % j1=fft(u1);
% % plot(f,10*log10(abs(j1(n+1)*2/N)));
% % title('调制噪声功率谱');
% % rand('state',0);
% % 
% % subplot(2,2,3);
% % plot(t1,y);
% % title('噪声调幅干扰时域波形');
% % %axis([0,50e-6,-10,10]);
% % 
% % subplot(2,2,4);

% 
% %干扰下的信号
% y=y(:,1:1000);
% xn=lfm+y;%雷达接收机信号
% 
%
% 
% 
% 
% 
% figure(2);
% y2=fft(xn);
% mag3=abs(y2);
% plot(f,mag3);
% xlabel('MHZ');
% 
% [M I] = max(abs(y2));
% 
% estmated_f = (I-1)*fs/N;
% 
% y3 = xn .* exp(-j*2*pi*estmated_f*t);
% F_y3 = fft(y3);
% 
% figure(3);
% plot(f*1e-6,abs(F_y3));
% xlabel('MHZ');








%x1=sawtooth(2*pi*fj*t);
% figure
% subplot(3,1,1);
% plot(t,x1);
% title('合成锯齿波信号,初相为0');
% axis([0 1e-5 -1 1]);





% %互相关求干扰初相
% r=xcorr(x1,s1);
% subplot(3,1,3);
% plot(r);
% title('互相关局部图');
% axis([150 200 -100 100]);
% grid



