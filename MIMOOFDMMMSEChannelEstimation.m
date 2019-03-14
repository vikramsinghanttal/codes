%%%%%%%******************* Vikram Singh *******************%%%%%
%%%%%%%------------------    16104100     ----------------%%%%%%
%%%%%%%$$$$$$$ Multimedia Wireless and Networks Lab $$$$$$%%%%%%
%%%%%%%+++++ Indian Innstitute of Technology, Kanpur +++++%%%%%%


clear all;
close all;


%% System Performance Parameters
Nt=4;
Nr=8;
L =10;
Cp=L-1;
Nc=12;
Nrays=1;
M=16;
N_fft=1024;
SNRd_dB=10;
SNRp_dB=(0:4:40)';
% SNRp_dB = 10;
SNRd=10.^(SNRd_dB*0.1);
SNRp=10.^(SNRp_dB*0.1);
montecarloiterations = 10;
% For CE
Np=2*max(Nr,Nt);
Codebook_Size=1024;
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Auxilliary variable definiition just for speed up

%% Transmitter Side Parameters
Xk_unmod     =   zeros(Nt,N_fft);  
X_k	=	zeros(Nt,N_fft);
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
X_k     =   zeros(Nt,N_fft);
Xk_bar  =   zeros(Nt,N_fft);
H_k_cap =   zeros(Nt,Nr);
%For channel Estimation
Y_check = zeros(0,0);
X_check = zeros(0,0);
Rx_buf_Pilot = zeros(Nr,Np*N_fft);
Tx_buf_Pilot = zeros(Nt,Np*N_fft);
H_aggregate  = zeros(Nr,Nt*N_fft);
Haggregate_cap = zeros(Nr,Nt*N_fft);
%%


MSE = zeros(length(SNRp),1);
SER = zeros(length(SNRp),1);
for snr=1:length(SNRp)
    for mci = 1:montecarloiterations        
%         Channel is constant for one coherence period
%         Pilots and data are transmitted during this coherence period
        H     =   sqrt(.5)*(randn(Nr,Nt*L)+1j*randn(Nr,Nt*L));
        
        % For estimating Channel
        for n_p = 1:Np
            %Tx Side Processing
            Xk_unmod   =   randi([0 M-1],Nt, N_fft);                   % Symbol generated in Frequency domain
            X_k =   qammod(Xk_unmod,M);                              % M-QAM Modulation
            x_n   =   ifft(X_k.',N_fft).';                      % Taking for time domain transmission
            x_cp  =   [x_n(:,N_fft-Cp+1:N_fft) x_n zeros(Nt,1)];              % Adding Cyclic prefix to counter
                                                                    % Multipath Interference
            %Propagation over Channel
%             H     =   sqrt(.5)*(randn(Nr,Nt*L)+1j*randn(Nr,Nt*L));%L-path (Nr,Nt) rayleigh faded MIMO Channel

            Tx_buf(:,1:L) = fliplr(x_cp(:,1:L));                  % Tx buffer will transmit Nt bits every instant
                                                                  % and maintain the record of path L transmitted vectors
            for pkt = 1:N_fft
                y=H*reshape(Tx_buf,[],1)+sqrt(.5/SNRp(snr))*(randn(Nr,1)+1j*randn(Nr,1));
                Tx_buf(:,2:L)=Tx_buf(:,1:L-1);
                Tx_buf(:,1)=x_cp(:,L+pkt);
                Rx_buf(:,pkt) = y;
            end
            Y_k=fft(Rx_buf.').';
            Rx_buf_Pilot(:,Np*(0:1:N_fft-1)+n_p)=Y_k;
            Tx_buf_Pilot(:,Np*(0:1:N_fft-1)+n_p)=X_k;
        end
        
%       Receive Side Processing
        for k = 1:N_fft  
            H_k=FFT_of_Channel(H,Nt,L,N_fft,k);
            H_aggregate(:,(k-1)*Nt+1:k*Nt) = H_k;
            H_k_cap=MMSEChannelEstimation(Rx_buf_Pilot(:,(k-1)*Np+1:k*Np), Tx_buf_Pilot(:,(k-1)*Np+1:k*Np), SNRp(snr));
            Haggregate_cap(:,(k-1)*Nt+1:k*Nt) = H_k_cap;
        end

        MSE(snr) = MSE(snr)+sum(sum(abs(H_aggregate-Haggregate_cap).^2));

    end
end
Scaling_factor = Nt*Nr*N_fft*montecarloiterations;
MSE = MSE/Scaling_factor;
% SER = SER/Scaling_factor;

figure
semilogy(SNRp_dB, MSE,'b--o');
title('MSE vs SNR');
ylabel('Mean Square Error in H per path');
xlabel('SNR(dB)')
% legend('cos(x)','cos(2x)')
% figure
% semilogy(SNRp_dB, SER,'r--o');
% title('SER vs SNR')
% xlabel('Signal to Noise Ratio')
% ylabel('Symbol Error Rate')
% legend('cos(x)','cos(2x)')