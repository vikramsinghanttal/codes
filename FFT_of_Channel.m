function H_k = FFT_of_Channel(H,Nt,L,N_fft,k)
H_k=zeros(size(H,1),Nt);
    for l=0:L-1
       H_k=H_k+H(:,1+l*Nt:(l+1)*Nt)*exp(-1j*(2*pi*k*l)/N_fft);
    end
end