%%%%%%%******************* Vikram Singh *******************%%%%%
%%%%%%%------------------    16104100     ----------------%%%%%%
%%%%%%%$$$$$$$ Multimedia Wireless and Networks Lab $$$$$$%%%%%%
%%%%%%%+++++ Indian Innstitute of Technology, Kanpur +++++%%%%%%


clear all;
close all;


%% System Performance Parameters
Nt=4;
Nr=64;
L=10;
Cp=L-1;
Nc=12;
Nrays=1;
M=16;
N_fft=512;
SNRd_dB=10;
SNRp_dB=(0:2:20)';
SNRd=10.^(SNRd_dB*0.1);
SNRp=10.^(SNRp_dB*0.1);
montecarloiterations = 10;
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Auxilliary variable definiition just for speed up

%% Transmitter Side Parameters
X_f     =   zeros(Nt,N_fft);  
X_mod	=	zeros(Nt,N_fft);
x_n     =	zeros(Nt,N_fft);
x_cp	=	zeros(Nt,N_fft+L+1);
%%

%% L- path (Nr,Nt) Rayleigh faded MIMO Wireless Channel
H       =	zeros(Nr,Nt*L);
%%

%% Receiver Side Parameters
Tx_buf  =	zeros(Nt,L);
Rx_buf	=   zeros(Nr,N_fft);
xn_bar  =   zeros(Nt,N_fft);
Xmod_bar=   zeros(Nt,N_fft);
Xf_bar  =   zeros(Nt,N_fft);
Xk_mod  =   zeros(Nt,N_fft);
Xk_bar  =   zeros(Nt,N_fft);
H_k     =   zeros(Nt,Nr);
%%


MSE = zeros(length(SNRp),1);
SER = zeros(length(SNRp),1);
for snr=1:length(SNRp)
    for mci = 1:montecarloiterations
        %Tx Side Processing
        X_f   =   randi([0 M-1],Nt, N_fft);                   % Symbol generated in Frequency domain
        X_mod =   qammod(X_f,M);                              % M-QAM Modulation
        x_n   =   ifft(X_mod.',N_fft).';                      % Taking for time domain transmission
        x_cp  =   [x_n(:,N_fft-Cp+1:N_fft) x_n zeros(Nt,1)];              % Adding Cyclic prefix to counter
                                                                % Multipath Interference
        %Propagation over Channel
        H     =   sqrt(.5)*(randn(Nr,Nt*L)+1j*randn(Nr,Nt*L));%L-path (Nr,Nt) rayleigh faded MIMO Channel

        Tx_buf(:,1:L) = fliplr(x_cp(:,1:L));                  % Tx buffer will transmit Nt bits every instant
                                                              % and maintain the record of path L transmitted vectors
        for pkt = 1:N_fft
            y=H*reshape(Tx_buf,[],1)+sqrt(.5/SNRp(snr))*(randn(Nr,1)+1j*randn(Nr,1));
            Tx_buf(:,2:L)=Tx_buf(:,1:L-1);
            Tx_buf(:,1)=x_cp(:,L+pkt);
            Rx_buf(:,pkt) = y;
        end
        %H_hat = Channel_Estimation(y,X_pilot);
        Y_k=fft(Rx_buf.').';

        %Receive Side Processing
        xn_bar   = pinv(H)*Rx_buf;               % Estimate of x_n ? Remember 1st Nt are relevent rest are redundent vector
        Xmod_bar = fft(xn_bar(1:Nt,:).').';      % Estimate of X_mod
        Xf_bar   = qamdemod(Xmod_bar,M);         % Estimate of X_f

        
%         Performance Measures
        SER(snr) = SER(snr)+biterr(Xf_bar,X_f);	% Symbol Error Rate at SNR=SNRp
        MSE(snr) = MSE(snr)+sum(sum(abs(x_n-xn_bar(1:Nt,:)).^2));

    end
end
No_of_Symbols_transmitted = Nt*N_fft*montecarloiterations;
MSE = MSE/No_of_Symbols_transmitted;
SER = SER/No_of_Symbols_transmitted;

figure
semilogy(SNRp_dB, MSE,'b--o');
title('MSE vs SNR');
xlabel('Mean Square Error in H');
ylabel('Symbol Error Rate')
% legend('cos(x)','cos(2x)')
figure
semilogy(SNRp_dB, SER,'r--o');
title('SER vs SNR')
xlabel('Signal to Noise Ratio')
ylabel('Symbol Error Rate')
% legend('cos(x)','cos(2x)')