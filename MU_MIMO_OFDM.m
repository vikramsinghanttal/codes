%%%%%%%******************* Vikram Singh *******************%%%%%
%%%%%%%------------------    16104100     ----------------%%%%%%
%%%%%%%$$$$$$$ Multimedia Wireless and Networks Lab $$$$$$%%%%%%
%%%%%%%+++++ Indian Innstitute of Technology, Kanpur +++++%%%%%%


clear all;
%/usr/local/MATLAB/MATLAB_Production_Server/R2015a/bin

%% System Performance Parameters
Nr=2;
K =4;
Nt=Nr*K;

if Nt < K*Nr
 disp('Error : Nt < K*Nr :: Not possible for system to accomodate this much users');
 return;
end

L =10;
Cp=L-1;
Nc=12;
Nrays=1;
M=16;
N_fft=32;
SNRd_dB=10;
SNRp_dB=(0:4:40)';
% SNRp_dB = 10;
SNRd=10.^(SNRd_dB*0.1);
SNRp=10.^(SNRp_dB*0.1);
montecarloiterations = 1;
% For Channel Estimation
Np=2*max(Nr,Nt);
Codebook_Size=1024;
%For Data Txn
Nd=8*Np;
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Auxilliary variable definiition just for speed up

%% Transmitter Side Parameters
Xk_unmod     =   zeros(K,Nt,N_fft);  
X_k	=	zeros(K,Nt,N_fft);
x_n     =	zeros(K,Nt,N_fft);
x_cp	=	zeros(K,Nt,N_fft+L+1);
%%

%% L- path (Nr,Nt) Rayleigh faded MIMO Wireless Channel
H       =	zeros(K,Nr,Nt*L);
Hk_aggregate  = zeros(K,Nr,Nt*N_fft);
%%

%% Receiver Side Parameters
Tx_buf  =	zeros(K,Nt,L);
Rx_buf	=   zeros(K,Nr,N_fft);
Xk_bar  =   zeros(K,Nt,N_fft);
H_k_cap =   zeros(K,Nt,Nr);
%For channel Estimation
Rx_buf_Pilot = zeros(K,Nr,Np*N_fft);
Tx_buf_Pilot = zeros(NK,t,Np*N_fft);
Hkaggregate_cap = zeros(K,Nr,Nt*N_fft);
%%


MSE = zeros(K,length(SNRp),1);
No_of_Bits_in_error = zeros(K,length(SNRp),1);
No_of_Bits_Txed = zeros(K,1);
SER = zeros(K,length(SNRp),1);

for snr=1:length(SNRp)
    for mci = 1:montecarloiterations        
%% &&&&&&&&&&&&&&&&&&&&& mci-th Coherence Time &&&&&&&&&&&&&&&&&&&&&&&&&&&&
%% &&&&&&& Channel stays nearly constant for one coherence period &&&&&&&&&

%         Pilots and data are transmitted during this coherence period
        H     =   sqrt(.5)*(randn(K,Nr,Nt*L)+1j*randn(K,Nr,Nt*L));
        
%% ++++++++++++++++++++ For Channel Estimation +++++++++++++++++++++++++++
        for k = 1: K
            for n_p = 1:Np
    %% --------------------Tx Side Processing --------------------------------
    %             Xk_unmod   =   randi([0 M-1],Nt, N_fft);                   % Symbol generated in Frequency domain
    %             X_k =   qammod(Xk_unmod,M);                              % M-QAM Modulation
                X_k(k,:,:)   =   sqrt(.5)*(randn(Nt,N_fft)+1j*randn(Nt,N_fft));
                x_n(k,:,:)   =   ifft(X_k(k,:,:).',N_fft).';                      % Taking for time domain transmission
                x_cp(k,:,:)  =   [x_n(k,:,N_fft-Cp+1:N_fft) x_n(k,:,:) zeros(Nt,1)];              % Adding Cyclic prefix to counter
                                                                        % Multipath Interference

                Tx_buf(k,:,1:L) = fliplr(x_cp(k,:,1:L));                  % Tx buffer will transmit Nt bits every instant
                                                                      % and maintain the record of path L transmitted vectors

    %% ------------------ Propagation over Channel-----------------------------
                for pkt = 1:N_fft
                    y(k,:)=H*reshape(Tx_buf,[],1)+sqrt(.5/SNRp(snr))*(randn(Nr,1)+1j*randn(Nr,1));
                    Tx_buf(k,:,:)=[x_cp(:,L+pkt) Tx_buf(:,1:L-1)];
                    Rx_buf(k,:,pkt) = y;
                end
                Y_k(k,:)=fft(Rx_buf(k,:,:).').';
                Rx_buf_Pilot(k,:,Np*(0:1:N_fft-1)+n_p)=Y_k(k,:);
                Tx_buf_Pilot(k,:,Np*(0:1:N_fft-1)+n_p)=X_k(k,:,:);
            end
        end
%% ------------------ Receive Side Processing ----------------------------
        for k = 1:N_fft  
            H_k=FFT_of_Channel(H,Nt,L,N_fft,k);
            Hk_aggregate(:,(k-1)*Nt+1:k*Nt) = H_k;
            H_k_cap=MMSEChannelEstimation(Rx_buf_Pilot(:,(k-1)*Np+1:k*Np), Tx_buf_Pilot(:,(k-1)*Np+1:k*Np), SNRp(snr));
            Hkaggregate_cap(:,(k-1)*Nt+1:k*Nt) = H_k_cap;
        end

        MSE(snr) = MSE(snr)+sum(sum(abs(Hk_aggregate-Hkaggregate_cap).^2));

%% +++++++++++++++++++++++ For Data Transmission +++++++++++++++++++++++++++
        
        for n_d = 1:Nd
            %Tx Side Processing
            Xk_unmod =   randi([0 M-1],Nt, N_fft);                  % Symbol generated in Frequency domain
            X_k      =   qammod(Xk_unmod,M);                        % M-QAM Modulation
            x_n      =   ifft(X_k.',N_fft).';                       % Taking for time domain transmission
            x_cp     =   [x_n(:,N_fft-Cp+1:N_fft) x_n zeros(Nt,1)]; % Adding Cyclic prefix to counter
                                                                    % Multipath Interference

            Tx_buf(:,1:L) = fliplr(x_cp(:,1:L));                    % Tx buffer will transmit Nt bits every instant
                                                                    % and maintain the record of path L transmitted vectors
            for pkt = 1:N_fft
                y=H*reshape(Tx_buf,[],1)+sqrt(.5/SNRp(snr))*(randn(Nr,1)+1j*randn(Nr,1));
                Tx_buf(:,2:L)=Tx_buf(:,1:L-1);
                Tx_buf(:,1)=x_cp(:,L+pkt);
                Rx_buf(:,pkt) = y;
            end
            
            Y_k=fft(Rx_buf.').';
            
            for k = 1:N_fft  
                Xk_bar(:,k)=pinv(Hkaggregate_cap(:,(k-1)*Nt+1:k*Nt))*Y_k(:,k);
            end

            Xk_unmod_bar   = qamdemod(Xk_bar,M);         % Estimate of X_f

            % Performance Measures
            No_of_Bits_in_error(snr) = No_of_Bits_in_error(snr)+biterr(Xk_unmod_bar,Xk_unmod);	% Symbol Error Rate at Pilot SNR = SNRp(snr)

        end
    end
    snr
end
Scaling_factor = Nt*Nr*N_fft*montecarloiterations;
MSE = MSE/Scaling_factor;
No_of_Bits_Txed=Nt*N_fft*log2(M)*montecarloiterations*Nd;
SER = No_of_Bits_in_error./No_of_Bits_Txed;

figure
semilogy(SNRp_dB, MSE,'b--o');
title(['MSE vs SNR (N_{fft} =',num2str(N_fft),')']);
ylabel('Mean Square Error in H per tap');
xlabel('SNR(dB)')
figure
semilogy(SNRp_dB, SER,'r--o');
title(['SER vs SNR (N_{fft} =',num2str(N_fft),')'])
xlabel('Signal to Noise Ratio')
ylabel('Symbol Error Rate');