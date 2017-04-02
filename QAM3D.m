close all
clear
clc

T=2000; %array size of the carrier wave.
T1 = 1; % define time of interval, 1 seconds
t = 0:T-1/T; % define time
t = t*T1; % define time in seconds

f_1 = 60 ;
f_2 = 61 ;
f_3 = 62 ;
f1 = (sqrt(2/T))*cos(2*pi*f_1*t) ;
f2 = (sqrt(2/T))*cos(2*pi*f_2*t) ;
f3 = (sqrt(2/T))*cos(2*pi*f_3*t) ;

N = 10^3 ;% number of symbols
M = 8; % constellation size
k = log2(M); % bits per symbol

% defining the constellation for PAM constellation
alphai = [0-(nthroot(M,3))/2 0+(nthroot(M,3)/2)];
alphaj = [0-(nthroot(M,3))/2 0+(nthroot(M,3)/2)];
alphak = [0-(nthroot(M,3))/2 0+(nthroot(M,3)/2)];
k_8QAM = 1/sqrt(3);

Eb_N0_dB = -5:15; % multiple Es/N0 values
Es_N0_dB = Eb_N0_dB + 10*log10(k);
ref = 0:k-1;
map = bitxor(ref,floor(ref/2));
[tt, ind] = sort(map);

for ii = 1:length(Eb_N0_dB)
    ipBitReshape = reshape(rand(1,N*k,1)>0.5,k,N).';
    s = cat(3,(2*(ipBitReshape(:,[k-2])) - ones(N,1) )*f1,(2*(ipBitReshape(:,[k-1])) - ones(N,1) )*f2, (2*(ipBitReshape(:,[k])) - ones(N,1) )*f3) * k_8QAM;
    % Adding noise in each of the channels
    z = s(:,:,1) + 10^(-Es_N0_dB(ii)/20)*1/sqrt(2)* transpose(randn(1,N))* f1 + s(:,:,2) + 10^(-Es_N0_dB(ii)/20)*1/sqrt(2)*transpose(randn(1,N))*f2 + s(:,:,3) + 10^(-Es_N0_dB(ii)/20)*1/sqrt(2)*transpose(randn(1,N))*f3; % additive white gaussian noise
    
    % demodulation
    % rounding to the nearest alphabet after multiplying the received signal with each of the basis vectors and integrating over period T.
    ipHati = 2*floor((z(:,:)/k_8QAM) * transpose(f1)/2)+1;
    ipHati(ipHati>max(alphai)) = max(alphai);
    ipHati(ipHati<min(alphai)) = min(alphai);
    ipHatj = 2*floor((z(:,:)/k_8QAM) * transpose(f2)/2)+1;
    ipHatj(ipHatj>max(alphaj)) = max(alphaj);
    ipHatj(ipHatj<min(alphaj)) = min(alphaj);
    ipHatk = 2*floor((z(:,:)/k_8QAM) * transpose(f3)/2)+1;
    ipHatk(ipHatk>max(alphak)) = max(alphak);
    ipHatk(ipHatk<min(alphak)) = min(alphak);
    
    % Constellation to Decimal conversion
    ipBinHati = ((ipHati+1)/2); % LUT based
    ipBinHatj = ((ipHatj+1)/2); % LUT based
    ipBinHatk = ((ipHatk+1)/2); % LUT based
    
    %Counting the number of errors present.
    nBitErr(ii) = size(find(ipBitReshape(:,k-2)- ipBinHati),1) + size(find(ipBitReshape(:,k-1)- ipBinHatj),1) + size(find(ipBitReshape(:,k)- ipBinHatk),1);
end

%plotting simulation and theoretical bit errors
simBer = nBitErr/(N*k);
theoryBer = 0.5*erfc(sqrt((10.^(Eb_N0_dB/10)))) ;

figure
semilogy(Eb_N0_dB,theoryBer,'bs-','LineWidth',2);
hold on
semilogy(Eb_N0_dB,simBer,'mx-','LineWidth',2);
axis([-5 15 10^-5 1])
grid on
legend('theory', 'simulation');
xlabel('Eb/No, dB')
ylabel('Bit Error Rate')
title('Bit error probability curve for 8-QAM modulation');
