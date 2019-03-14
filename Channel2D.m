%function H = Channel2D(d,lambda,AoA,AoD)
clear all;
close;
d      = 0.005;
lambda = 0.01;
Nt     = 4;
Nr     = 8;
Nc     = 10;
AoA    = 180*rand(Nc,1)-90;
AoD    = 180*rand(Nc,1)-90;

Tx_streering_vector = steervec((0:d:(Nt-1)*d),AoD');
Rx_streering_vector = steervec((0:d:(Nr-1)*d),AoA'); 


alpha = sqrt(0.5)*(randn(Nc,1)+1j*randn(Nc,1));

H = Rx_streering_vector*diag(alpha)*Tx_streering_vector';
H_imp=zeros(Nr,Nt);
for i=1:Nc
    H_imp = H_imp+alpha(i)*Rx_streering_vector(:,i)*Tx_streering_vector(:,i)';
end

%Analysis of Channel
% H_k = fft2(H);
% H_k_bar = dftmtx(Nr)*H*dftmtx(Nt)';
% diff=H_k-H_k_bar;
% power_H=abs(H_k);
% [U,S,V]=svd(H)
% [Uk,Sk,Vk]=svd(H_k)
H_inv=pinv(H)
% end