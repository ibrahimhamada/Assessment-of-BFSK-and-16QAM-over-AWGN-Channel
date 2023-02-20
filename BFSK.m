clc;
clear;

%% 1) BFSK 
%% Define transmitted signal (BFSK)
fs = 64; %Sampling Frequency
df = 32; %Frequency Separation
Tb=1; %Bit duration = 1 msec
N = 100; %Number of bits
X_input= randi([0, 1],1,N);  %Binary signal  
M = 2;
m = log2(M);
ns=2;  %Number of samples per symbol
X_digit=[]; 
nb=10000;  %Number of points between two symbols (it's used to convert the symbols into continuous digital signal)
for i=1:N
    if X_input(i)==1
       x_temp=ones(1,nb);
    else
        x_temp=zeros(1,nb);
    end
    X_digit=[X_digit x_temp];
end
t_sig = Tb/nb : Tb/nb : N*Tb; %Time vector of continuous digital signal
%Plotting the input message signal
figure();
plot(t_sig,X_digit, 'LineWidth',2,'Color','black'); grid on; xlim([0 Tb*N]); ylim([-0.5 1.5]);
xlabel('Time(sec)'); ylabel('Amplitude(volt)'); title('Input Message Signal');

%% BFSK Modulation (BFSK)
BFSKMOD = fskmod(X_input, M, df, ns, fs); %BFSK Modulation

%% Constellation Diagram (BFSK)
scatterplot(BFSKMOD); grid on; 
title('Constellation Diagram of Transmitted BFSK Symbols');

%% Noise in the Communication Channel (BFSK)
SNR = 15+10*log10(m);
Y = awgn(BFSKMOD, SNR, 'measured');  %Adds white Gaussian noise to the Modulated signal

%% Plotting the Noisy Signal
scatterplot(Y); grid on; 
title('Constellation Diagram of Received Noisy BFSK Symbols');

%% BFSK Demodulation (BFSK)
X_demod=fskdemod(Y, M, df, ns, fs); %Demodulation 

%Convert the symbols into continuous digital signal
X_dem_sig = [];
for i=1:length(X_demod)
    if X_demod(i)==1
       x_temp_dem=ones(1,nb);
    else
        x_temp_dem=zeros(1,nb);
    end
    X_dem_sig=[X_dem_sig x_temp_dem];
end
t_sig_dem = Tb/nb : Tb/nb : length(X_demod)*Tb; %Time vector of continuous digital signal
figure();
plot(t_sig_dem,X_dem_sig, 'LineWidth',2,'Color','blue'); grid on; xlim([0 Tb*length(X_demod)]); ylim([-0.5 1.5]);
xlabel('Time(sec)'); ylabel('Amplitude(volt)'); title('BPSK Demodulated Signal');

%% PSD
power_Spec = spectrum.welch;
PSD_BFSK=psd(power_Spec,BFSKMOD); 
figure(); plot(PSD_BFSK); title('PSD of the Transmitted BFSK Signal'); grid on; 

%% Theoretical BER (BFSK)
EbNo=15;
BER_Theo=zeros(1,EbNo/0.1);
j=1;
for i=0:0.1:EbNo 
    BER_Theo(j)=0.25*erfc(sqrt(i/2));
    j=j+1;
end
Eb_No=0:0.1:15;
figure(); plot(Eb_No,BER_Theo); grid on; 
ylabel('BER'); xlabel('Eb/N0 (dB)'); title('Theoretical BER for BFSK');

%% Measured BER (BFSK)
k=1;
l=10000;
BER_Measured = zeros(1,EbNo/0.1);
for i=0:0.1:EbNo
    rand_bits = randi([0 1],l,1);
    BFSKmodul=fskmod(rand_bits, M, df, ns, fs);
    yy=awgn(BFSKmodul,i+10*log10(m),'measured');
    demod=fskdemod(yy,M,df,ns,fs);
    [n ,e_ratio]=biterr(demod,rand_bits);
    BER_Measured(k)=e_ratio;
    k=k+1;
end
figure(); plot(Eb_No,BER_Measured); grid on; 
ylabel('BER'); xlabel('Eb/N0 (dB)'); title('Measured BER for BFSK');