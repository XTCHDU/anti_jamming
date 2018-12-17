close all
clear all
T = 5e-6;
fADC = 50e6;
t = 0:1/fADC:T - 1/fADC;
f1 = 5e6;
f2 = 15e6;

x1 = exp(j*2*pi*f1*t);
x1 = [zeros(1,100),x1,zeros(1,100)];
x2 = exp(j*2*pi*f2*t);
x2 = [zeros(1,100),x2,zeros(1,100)];
x = x1 + x2;

figure;
plot(real(x))

xfft = fft(x);
figure;
plot(abs(xfft))

xfft(105:165) = 0;
figure;
plot(abs(xfft))

y = ifft(xfft);
figure;
plot(real(y))