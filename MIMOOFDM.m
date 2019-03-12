%MIMO-OFDM
clear all;
close all;
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
mse = zeros(length(SNRp),1);
ser = zeros(length(SNRp),1);
montecarloiterations = 10;
for snr=1:length(SNRp)
    for mci= 1:montecarloiterations
    X_f         =   randi([0 M-1],Nt, N_fft);
    X_mod       =   qammod(X_f,M);
%     norm_x        =   sqrt(sum(abs(X_mod).^2));
%     Xmod_norm   =   bsxfun(@rdivide,X_mod,norm_x);
    x_n = ifft(X_mod.',N_fft).';

    x_cp = [x_n(:,N_fft-Cp+1:N_fft) x_n];

    H = sqrt(.5)*(randn(Nr,Nt*L)+1j*randn(Nr,Nt*L));

    xTx_buff = zeros(Nt,L);
    xTx_buff(:,1:L)= fliplr(x_cp(:,1:L));
    x_cp = [x_cp zeros(Nt,1)];
    Rx_buff = zeros(Nr,N_fft);
        for pkt = 1:N_fft
            y=H*reshape(xTx_buff,[],1)+sqrt(.5/SNRp(snr))*(randn(Nr,1)+1j*randn(Nr,1));
            xTx_buff(:,2:L)=xTx_buff(:,1:L-1);
            xTx_buff(:,1)=x_cp(:,L+pkt);
            Rx_buff(:,pkt) = y;
        end

    xn_bar = pinv(H)*Rx_buff;
    Xfmod_bar = fft(xn_bar(1:Nt,:).').';
    Xf_bar = qamdemod(Xfmod_bar,M);
    ser(snr)= biterr(Xf_bar,X_f);
    mse(snr) = mse(snr)+sum(sum(abs(x_n-xn_bar(1:Nt,:)).^2));

    end
end
mse = mse/(Nt*N_fft*montecarloiterations);
ser = ser/(Nt*N_fft*montecarloiterations);
figure
semilogy(SNRp_dB, mse,'b--o');
title('MSE vs SNR');
xlabel('Mean Square Error in H');
ylabel('Symbol Error Rate')
% legend('cos(x)','cos(2x)')
figure
semilogy(SNRp_dB, ser,'r--o');
title('SER vs SNR')
xlabel('Signal to Noise Ratio')
ylabel('Symbol Error Rate')
% legend('cos(x)','cos(2x)')