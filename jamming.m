function [mag,y] = jamming(N,type,fADC,Amp,B,T,jam_freq)
mag = zeros(1,N);
switch type
    case 1
        y = (normrnd(0,Amp/4,1,N)+Amp/2).*exp(1j*(2*pi*jam_freq*(0:N-1)/fADC));
    case 2
        n = normrnd(0,2,1,N);
        for i = 2:length(n)
            n(i) = n(i)+n(i-1);
        end
        y = Amp*exp(1j*(2*pi*5e6*(0:N-1)/fADC+2*pi*B/T.*n));
        mag = (2*pi*5e6*(0:N-1)/fADC+2*pi*B/T.*n);
end
