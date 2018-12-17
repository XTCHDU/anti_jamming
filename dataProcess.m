function y = dataProcess(sig,N,type,fADC)
switch type
    case 2
        lg_sig = log(sig);
        Rx = real(lg_sig);
        Ix = imag(lg_sig);
        base = Ix(1);
        for i = 2:N
            if Ix(i)>base+pi
                Ix(i) = Ix(i)-2*pi;
            else
                if Ix(i)<base-pi
                    Ix(i) = Ix(i)+2*pi;
                end
            end
        end
        A_ = exp(mean(Rx));
        J = A_*exp(1j*(Ix));
        y = sig - J;
        
    case 1
        lg_sig = sig;
        FT_sig = fft(sig);
        FT_sig(90:110) = 0;
        y = ifft(FT_sig);
        return
        FFT_N = 2^(ceil(log2(length(sig)))+1);
        FFT_N = 750;
        FT_log_sig = fftshift(abs(fft(lg_sig,FFT_N)));
        fAxis = linspace(-fADC/2,fADC/2,FFT_N);
        base = FFT_N/2+1;
        FT_right = FT_log_sig(base:end);
        [~,index] = max(FT_right);
        W_ = fAxis(base+index-1);
        W_2 = W_;
        threshold = 1.5e6;
        while (abs(W_2-W_)<threshold)
            FT_right(index) = 0;
            FT_log_sig(base+index-1-10:base+index-1+10) = 0;
            [~,index] = max(FT_right);
            W_2 = fAxis(base+index-1);
        end
        y = ifft(ifftshift(FT_log_sig),length(sig));
        return 
        %W_ = 13e6;
        for i = 1:length(sig)
            re_j(i) = abs(sig(i)).*cos(2*pi*W_*(0:N-1)/fADC)%exp(1j*2*pi*W_*(i-1)/fADC);
        end
        y = sig-re_j;
        if W_2<W_
            Wstop=2*(W_-threshold)/fADC;
            [b,a]=butter(30,Wstop,'low');
            y=filter(b,a,y);
        else
            Wpass = 2*(W_+threshold)/fADC;
            [b,a] = butter(5,Wpass,'high');
            y = filter(b,a,y);
        end        
end


