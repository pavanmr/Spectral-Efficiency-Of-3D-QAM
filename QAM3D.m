close all
clear
clc

N = 10^5;       % number of symbols
M = 8;          % constellation size
k = log2(M);    % bits per symbol
T=1000;         %samples
tb=1;           %bit period

t=(0:T-1)/T;
t=t*tb;
f1=150;
f2=151;
f3=152;

%orthomormal basis functions
f1_t=sqrt(2/T)*cos(2*pi*f1*t);
f2_t=sqrt(2/T)*cos(2*pi*f2*t);
f3_t=sqrt(2/T)*cos(2*pi*f3*t);

k_8QAM = 1/sqrt(3);

Eb_N0_dB  = [-5:15]; % multiple Es/N0 values
Es_N0_dB  = Eb_N0_dB + 10*log10(k);

for ii = 1:length(Eb_N0_dB)
    
    % symbol generation
    % ------------------
    ipBit = rand(1,N*k,1)>0.5; % random 1's and 0's
    ipBitReshape = reshape(ipBit,k,N).';
    bin2DecMatrix = 2*ipBitReshape-ones(N,k);
    % conversion from binary to decimal
    %----------------------------------
    s1=bin2DecMatrix(:,1)*f1_t;
    s2=bin2DecMatrix(:,2)*f2_t;
    s3=bin2DecMatrix(:,3)*f3_t;
    
    s1_t=s1*k_8QAM;
    s2_t=s2*k_8QAM;
    s3_t=s3*k_8QAM;
    
    %noise generation
    %----------------
    n1=1/sqrt(2)*randn(1,N);
    n2=1/sqrt(2)*randn(1,N);
    n3=1/sqrt(2)*randn(1,N);
    
    n1_t=transpose(n1)*f1_t;
    n2_t=transpose(n2)*f2_t;
    n3_t=transpose(n3)*f3_t;
    
    y1=s1_t+ 10^(-Es_N0_dB(ii)/20)*n1_t;
    y2=s2_t+10^(-Es_N0_dB(ii)/20)*n2_t;
    y3=s3_t+10^(-Es_N0_dB(ii)/20)*n3_t;
    y = y1 + y2 + y3;
    
    %demodulation
    %------------
    y1_t=(y*transpose(f1_t))/k_8QAM;
    y2_t=(y*transpose(f2_t))/k_8QAM;
    y3_t=(y*transpose(f3_t))/k_8QAM;
    
    
    % Constellation to Decimal conversion
    %------------------------------------
    ipHat1 = [2*floor(y1_t/2)]+1;
    ipHat1(ipHat1>max(bin2DecMatrix(:,1)))= max(bin2DecMatrix(:,2));
    ipHat1(ipHat1<min(bin2DecMatrix(:,1)))=  min(bin2DecMatrix(:,1));
  
    ipHat2 = 2*floor(y2_t/2)+1;
    ipHat2(ipHat2>max(bin2DecMatrix(:,2))) = max(bin2DecMatrix(:,2));
    ipHat2(ipHat2<min(bin2DecMatrix(:,2))) = min(bin2DecMatrix(:,2));
    
    ipHat3 = 2*floor(y3_t/2)+1;
    ipHat3(ipHat3>max(bin2DecMatrix(:,3))) = max(bin2DecMatrix(:,3));
    ipHat3(ipHat3<min(bin2DecMatrix(:,3))) = min(bin2DecMatrix(:,3));
    
    % converting to binary string
    ipHat1=(ipHat1+1)/2;
    ipHat2=(ipHat2+1)/2;
    ipHat3=(ipHat3+1)/2;
    
    % counting errors
    nBitErr(ii)=size(find(ipBitReshape(:,1)-ipHat1),1)+size(find(ipBitReshape(:,2)-ipHat2),1)+size(find(ipBitReshape(:,3)- ipHat3),1);
end

% Simulation error
simBer = nBitErr/(N*k);

% Theoretical error
Berr_theory = 0.5*erfc(sqrt(10.^(Eb_N0_dB/10)));

EbNo_lin=10.^(Eb_N0_dB/10);

figure
semilogy(Eb_N0_dB,Berr_theory,'bs-','LineWidth',2);
hold on
semilogy(Eb_N0_dB,simBer,'mx-','LineWidth',2);
axis([0 15 10^-5 1])
grid on
legend('theory', 'simulation');
xlabel('Eb/No, dB')
ylabel('Bit Error Rate')
title('Bit error probability curve for 8-QAM modulation')
