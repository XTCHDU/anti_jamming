function y = dataProcess(sig,N,type,fADC,norm_f,Amp)
switch type
    case 2
        Rx = real(sig);
        Ix = imag(sig);
        abs_x = abs(sig).^2;
        Ex = mean(abs_x);
        Varx = var(abs_x);
        x1 = sqrt((2*Ex+sqrt(4*Ex^2-8*Varx))/4)
        x2 = sqrt((2*Ex-sqrt(4*Ex^2-8*Varx))/4)
        J = x1*exp(1j*norm_f(2,:))+x2*exp(1j*norm_f(1,:));
        y = sig-J;
    case 1
        lg_sig = sig;
        FT_sig = fft(sig);
        [~,index] = max(abs(FT_sig));
        FT_sig(index:index) = 10;
        y = ifft(FT_sig);
        return
        FFT_N = 2^(ceil(log2(length(sig)))+1);
        %FFT_N = length(sig);
        FT_log_sig = fftshift(abs(fft(lg_sig,FFT_N)));
        fAxis = linspace(-fADC/2,fADC/2,FFT_N);
        base = FFT_N/2+1;
        FT_right = FT_log_sig(base:end);
        [~,index] = max(FT_right);
        W_ = fAxis(base+index-1);
        indexW_ = base+index-1;
        W_2 = W_;
        threshold = 1e6;
        while (abs(W_2-W_)<threshold)
            FT_right(index) = 0;
            FT_right(base+index-1-2:base+index-1+2) = 0;
            [~,index] = max(FT_right);
            W_2 = fAxis(base+index-1);
        end
        %return 
        W_=12e6;
        for i = 1:length(sig)
            re_j(i) = abs(sig(i))*exp(1j*2*pi*W_*(i-1)/fADC);
        end
        y = sig-re_j;
        fft_re_sig = fft(y,FFT_N);
        %y = solution(fft_re_sig,indexW_,jam);
        %return 
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


