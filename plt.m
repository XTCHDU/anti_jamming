lfm = load('jsr_lfm.mat');
lfm = lfm.res;
nlfm = load('jsr_nlfm.mat');
nlfm = nlfm.res;
mp = load('jsr_mp.mat');
mp = mp.res;
bpsk = load('jsr_bpsk.mat');
bpsk = bpsk.res;
x = 0:5:50;
plot(x,lfm(1:length(x)),'-*');
hold on;
plot(x,nlfm(1:length(x)),'-*');
hold on;
mp(1) = -10;
plot(x,mp(1:length(x)),'-*');
hold on;
plot(x,bpsk(1:length(x)),'-*');
xlabel('抑制前JSR（dB）')
ylabel('抑制后JSR（dB）')
grid on;
legend('LFM','NLFM','MP','BPSK')