%function H = Channel3D(d,lambda,AoA,AoD,ZoA,ZoD)
clear all;
close;
d      = 0.005;
lambda = 0.01;
Nt     = 4;
Nr     = 4;
Nc     = 10;
AoA    = 180*rand(Nc,1)-90;
AoD    = 180*rand(Nc,1)-90;
ZoA    = 180*rand(Nc,1)-90;
ZoD    = 180*rand(Nc,1)-90;

Tx_streering_vector_AoD = steervec((0:d:(Nt-1)*d),AoD');
Rx_streering_vector_AoA = steervec((0:d:(Nr-1)*d),AoA'); 
Tx_streering_vector_ZoD = steervec((0:d:(Nt-1)*d),ZoD');
Rx_streering_vector_ZoA = steervec((0:d:(Nr-1)*d),ZoA');

Tx_streering_vector = zeros(Nt*Nt,Nc);
Rx_streering_vector = zeros(Nr*Nr,Nc);

alpha = sqrt(.5)*(randn(Nc,1)+1j*randn(Nc,1));

for i=1:length(Nc)
    Tx_streering_vector(:,i)=kron(Tx_streering_vector_AoD(:,i),Tx_streering_vector_ZoD(:,i));
    Rx_streering_vector(:,i)=kron(Rx_streering_vector_AoA(:,i),Rx_streering_vector_ZoA(:,i));
end

H = Tx_streering_vector*diag(alpha)*Rx_streering_vector';

% end