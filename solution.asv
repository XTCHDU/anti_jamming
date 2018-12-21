function [y] = solution( sig, index,jam )
s = abs(fft(jam,length(sig)));
y = sig;
base = length(sig)/2;
y_left = y(1:base);
y_right = y(base+1:end);
y_final = y_left-fliplr(conj(y_right));
end

