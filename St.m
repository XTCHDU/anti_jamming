clear all 
close all
JSR = 10;

%gan rao
fj=1e5;
bj=3e4;                   
Tr=520e-6; 
fs=10e5;
t1=0:1/fs:3*Tr-1/fs;
N=length(t1);
 uj=1;
 mf=0.6; 
 wpp=6;     
u=wgn(1,N,wpp);               
wp=10e6;                        
ws=13e6; 
rp=1; rs=60; 
[Nn,wn]=buttord(wp/(30e6/2),ws/(30e6/2),rp,rs); 
[b,a]=butter(Nn,wn);           
u1=filter(b,a,u);              
ss(1)=0;                       
for i=1:N-1                    
ss(i+1)=u1(i)+ss(i);
end
%   i=1:N-1; ss=cumsum([0,u1(i)]);
ss=ss*Tr/N;                        
y0=uj*cos(2*pi*fj*t1+2*pi*mf*bj*ss+2*pi*rand(1,1));  

%LFM 
ffs=10e5;
[lfmwav0] = LFM(ffs) ;

%gan  xin bi 
pj=(1/N)*sum(y0.^2);  
plfm=(1/N)*sum(lfmwav0.^2);  
tem= 10^(JSR/10)*plfm;
y0 = sqrt(tem/pj).*y0;
pJamming = (1/N)*sum(y0.^2);
JSR1= 10*log10(pJamming/plfm);

%qu dian
num = min(length(lfmwav0),length(y0)) ;
la1 = linspace(1 , length(lfmwav0) , num) ;
lfmwav = lfmwav0(round(la1)) ;
la2 = linspace(1 , length(y0) , num) ;
y =y0(round(la2)) ;

%%
% plot(t1(round(la2)),y,'r');    
% %  axis([0,0.05e-4,-6,6]);
% plot(t1(round(la2)),real(lfmwav),'b');  
% % axis([0,0.05e-4,-6,6]); 
% hold off
%%
%xinhao + ganrao
lfmwav=lfmwav';
xt=lfmwav+y;  
figure(4)
plot(real(xt));

% axis([0,1.5e-3,-6,6]); t1(round(la2)),
% slo=[zeros(1,400)  lfmwav];
% slo1=conj(slo);
% yt=xt.*slo1;
% r=fftfli(yt);
% figure(5);
% plot(fftshift(abs(r)));
% figure(6);
% lowpass(yt,250e3,fs);
% temp=lowpass(yt,250e3,fs);
% srt1=temp;
% figure(7)
% plot(real(srt1));
 % srt1=real(srt1);
 
%emd 
[imf , f2c , c2f] = emd(real(xt)) ;
emd_visu(real(xt),t1(round(la2)),imf,1) ;
 temp=imf(7,:)+imf(6,:)+imf(5,:)+imf(4,:)+imf(3,:) +imf(2,:);
% temp=imf(2,:);
figure(5);
plot(real(temp));
title('')

%mai ya
lfmwav1=conj(fliplr(lfmwav));
w=conv(lfmwav1,temp);
figure(6)
plot(real(w))

w1=conv(lfmwav1,lfmwav);
figure(7)
plot(real(w1))



