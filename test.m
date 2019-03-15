Nt = 48;
Nr = 64;
L = 100;
k=4;
N_fft = 128;
H_k2=zeros(Nr,Nt);
H_k1=zeros(Nr,Nt);
mci=1000;

tic
for i=1:mci
H   =sqrt(.5)*(randn(Nr,Nt*L)+1j*randn(Nr,Nt*L));
H_k1=(reshape(H.*repelem(exp((-1j*2*pi*k/N_fft)*(0:1:L-1)),Nt),[Nr*L,Nt]));
end
toc

%{
tic
for i=1:mci
H     =   sqrt(.5)*(randn(Nr,Nt*L)+1j*randn(Nr,Nt*L));
for l=0:L-1
       H_k2=H_k2+H(:,1+l*Nt:(l+1)*Nt)*exp(w*l);
end
end
toc
%}