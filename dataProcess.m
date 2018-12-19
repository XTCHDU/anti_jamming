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
        %return 
        W_ = 13e6;
        W_2 = 5e6;
        for i = 1:length(sig)
            re_j(i) = abs(sig(i))*exp(1j*2*pi*W_*(i-1)/fADC);
        end
        y = sig-re_j;
        if W_2<W_
            Wstop=2*(W_-threshold)/fADC;
            wp=Wstop;ws=2*W_/fADC;rp=0.1;rs=60;   %DF指标（低通滤波器的通、阻带边界频）
            [N,wp]=ellipord(wp,ws,rp,rs); %调用ellipord计算椭圆DF阶数N和通带截止频率wp
            [b,a]=ellip(N,rp,rs,wp);

            y=filter(b,a,y);
        else
            Wpass = 2*(W_+threshold)/fADC;
            wp = Wpass;ws = 2*W_/fADC;rp=0.1;rs=60;   %DF指标（低通滤波器的通、阻带边界频）
            [N,wp]=ellipord(wp,ws,rp,rs);    %调用ellipord计算椭圆DF阶数N和通带截止频率wp
            [b,a]=ellip(N,rp,rs,wp,'high'); %调用ellip计算椭圆带通DF系统函数系数向量B和A
            y = filter(b,a,y);
        end            
        y = mapminmax(real(y),-1,1);
end


