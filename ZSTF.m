%function y=noiseAM(u0,N,wpp); 
u0=1;
noise_power_dB = 10; 
fj=35e6;
fs=4*fj; 
Tr=520e-6;     
t1=0:1/fs:3*Tr-1/fs; 
N=length(t1); 
u=wgn(1,N,noise_power_dB); 
df1=fs/N;
n=0:N/2;
f=n*df1;   
wp=10e6; 
ws=14e6;  
rp=1; 
rs=60;  
[n1,wn1]=buttord(wp/(fs/2),ws/(fs/2),rp,rs);  
[b,a]=butter(n1,wn1);  
u1=filter(b,a,u);   %得到带限噪声

figure,plot(abs(fftshift(fft(u1))));


p1=(1/N)*sum(u1.^2);   %带限噪声的平均功率



figure;
y=(u0+u1).*cos(2*pi*fj*t1+2*pi*rand(1,1));     
p=(1/N)*sum(y.^2);   %噪声调幅干扰的平均功率=u0^2/2+p1^2/2
subplot(2,1,1), plot(t1,y),title('噪声调幅干扰时域波形'); axis([0,0.05e-4,-2,2])
subplot(2,1,2), J=fft(y);plot(f,10*log10(abs(J(n+1))));
title('已调波功率谱');
figure;
% subplot(2,2,1),plot(t1,u1),title('噪声调制波形'); axis([0,0.05e-4,-2,2])
% subplot(2,2,2), j2=fft(u1);plot(f,10*log10(abs(j2(n+1)*2/N)));
% title('调制噪声功率谱');
% rand('state', 0); 

