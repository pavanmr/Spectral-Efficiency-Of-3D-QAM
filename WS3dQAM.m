close all
clear
clc

N = 100000;     % number of symbols
M = 8;          % constellation size
k = log2(M);    % bits per symbol
T=2000;         %samples
tb=1;           %bit period
t=(0:T-1)/T;
t=t*tb;

f1=120;
f2=240;
f3=480;

%orthomormal basis functions
f1_t=sqrt(2/T)*cos(2*pi*f1*t);
f2_t=sqrt(2/T)*cos(2*pi*f2*t);
f3_t=sqrt(2/T)*cos(2*pi*f3*t);

k_8QAM = 1/sqrt(3);

% symbol generation
% ------------------
ipBit = rand(1,N*k,1)>0.5; % random 1's and 0's
ipBitReshape = reshape(ipBit,k,N).';
bin2DecMatrix = 2*ipBitReshape-ones(N,k); % conversion from binary to decimal
s1=bin2DecMatrix(:,1)*f1_t;
s2=bin2DecMatrix(:,2)*f2_t;
s3=bin2DecMatrix(:,3)*f3_t;

s1_t=sqrt(T/2)*s1*k_8QAM;
s2_t=sqrt(T/2)*s2*k_8QAM;
s3_t=sqrt(T/2)*s3*k_8QAM;
s_t=s1_t+s2_t+s3_t;
g_t=mean(s_t);

% generating wave shaping function
gt=ones(1,T);
gt(T/2:T)= 0;
h = gt(T/2:T);
raisedcosine=sin(2*pi*5*t);
g1=raisedcosine.*gt ;

%convolving
conv_t=conv(g1,g_t);
len_t=length(conv_t);

%power calculation
p = abs(fft(conv_t))/(len_t/2); %% absolute value of the fft
p = p(1:len_t/2).^2; %% take the power of positve freq. half
freq = (0:len_t/2-1)/tb; %% find the corresponding frequency in Hz
%plot of power
semilogy(freq,p); %% plot on semilog scale
axis([0 500 0 1]); %% zoom in
