function H_k = MMSEChannelEstimation(Y, X, SNR,Nt)
H_k=Y*X'/((1/SNR)*eye(Nt)+X*X');
end