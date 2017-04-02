close all
clear
clc

T=2000;
T1 = 1; %% define time of interval, 1 seconds
t = [0:T-1]/T; %% define time
t = t*T1; %% define time in seconds
f_1 = 10 ;
f_2 = 200 ;
f_3 = 980 ;

pie = 3.1416 ;
f1 = (sqrt(2/T))*cos(2*pi*f_1*t) ;
f2 = (sqrt(2/T))*cos(2*pi*f_2*t) ;
f3 = (sqrt(2/T))*cos(2*pi*f_3*t) ;

g1 = sin(2*pi*5*t) ;
g = ones(1,T);
g(T/10:T) = 0 ;
g = g1.*g ;
N = 60000 ;% number of symbols
M = 8; % constellation size
k = log2(M); % bits per symbol
% defining the real and imaginary PAM constellation
alphai = [0-(nthroot(M,3))/2 0+(nthroot(M,3)/2)];
alphaj = [0-(nthroot(M,3))/2 0+(nthroot(M,3)/2)];
alphak = [0-(nthroot(M,3))/2 0+(nthroot(M,3)/2)];
k_8QAM = 1/sqrt(3);

Eb_N0_dB = [-5:15]; % multiple Es/N0 values
Es_N0_dB = Eb_N0_dB + 10*log10(k);

ref = [0:k-1];
map = bitxor(ref,floor(ref/2));
[tt ind] = sort(map);

ipBitReshape = reshape(rand(1,N*k,1)>0.5,k,N).';
s = cat(3,(2*(ipBitReshape(:,k-2)) - ones(N,1) )*f1,2*(ipBitReshape(:,[k-1])) - ones(N,1)*f2, (2*(ipBitReshape(:,[k])) -ones(N,1) )*f3 * k_8QAM); % normalization of transmit power to one
z = s(:,:,1) + 10^(-Es_N0_dB(1)/20)*1/sqrt(2)* transpose(randn(1,N))* f1 + s(:,:,2) + 10^(-Es_N0_dB(1)/20)*1/sqrt(2)*transpose(randn(1,N))*f2 + s(:,:,3) + 10^(-Es_N0_dB(1)/20)*1/sqrt(2)*transpose(randn(1,N))*f3; % additive white gaussian noise
avg = (1/N)*sum(z,1);

z_conv = conv(avg,g);
length_conv = length(z_conv);
disp(length_conv) ;

p = abs(fft(z_conv))/(length_conv/2); %% absolute value of the fft
p = p(1:length_conv/2).^2 ;%% take the power of positve freq. half
freq = [0:length_conv/2-1]/T1; %% find the corresponding frequency in Hz
semilogy(freq,p); % plot on semilog scale
axis([0 1000 0 1]); % zoom in
xlabel('Frequency, Hz')
ylabel('Power spectral density')
title('Power spectral density of implemented system on a semi-log plot using a wave shaping function');