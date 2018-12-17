function y = get_sum(N,X)
y = zeros(1,N);
    for i = 1:length(X)
        y = y+jamming(N,X(i));
    end
    figure;
    plot(y);
    title('get sum')