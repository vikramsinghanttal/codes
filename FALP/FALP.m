clear all;
N=4;
Kc=10;
H=channel(N,Kc);
X1=ifft2(H);
% X2=dftmtx(N)'*H*dftmtx(N);
abs(X1)/sum(sum(abs(X1)))*100