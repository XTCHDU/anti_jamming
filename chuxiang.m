%%
%%% �ο���Դ
%%% http://www.doc88.com/p-7856386736794.html
%%
clear all
% load('k.mat') ;
k = 1 : 3155 ;
% M = 110 ;
M = 631 ;
N = 128 ;
cx = [] ;
for kk = 1 : 5
    Y = [] ;
    for i = (M * (kk - 1) + 1) : M * kk
        load(['C:\Users\94591\Documents\MATLAB\����matlab����\IF_' , num2str(k(i)) , '.mat']);
        a = sigout ;
        x = hilbert(a);
        am = abs(a + 1i * x) ;%����
        d = ceil(length(a) / N) ; %���ݼ��
        nn = 5 ; %ƽ��ָ��
        l = N * d ;
        tau = l - length(a) ;
        a = [a, zeros(1 , tau)] ; 
        y = reshape(a , d , N) ;
        y = max(y) ;
        [l1 l2] = size(y) ;
        y = reshape(y , 1 , l1 * l2) ;
        x = linspace(0 , length(a) , N) ;
        y1 = smooth(y , nn) ;
        Y = [Y ; y1'] ;
    end
    Y = 1 / sqrt(N) * Y ;
    R1 = Y * conj(Y') ;
    [u lamda] = eig(R1) ;
    lam = diag(lamda) ;  %ȡ���Խ�Ԫ��
    [lamda1 ind] = max(lam)  ; %lamda1Ϊ�������ֵ��ihd��Ӧ�������ֵ����λ��
    u1 = u(:,ind(1)) ;  %�������ֵ��Ӧ����������
    A = sqrt(M) * u1 ;
    cx = [cx ; A] ;
end