function H_k = MMSEChannelEstimation(Y, X, SNR)
H_k=Y*X'/((1/SNR)*eye(size(X,1))+X*X');
end