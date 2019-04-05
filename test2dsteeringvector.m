clear all;
Kc=10; 
N=4;
AoD = -90+rand(Kc,1)*90;
ZoD = -90+rand(Kc,1)*90;
gama= sqrt(0.5)*(randn(Kc,1)+1j*randn(Kc,1));
pi_sin_AoD=pi*sin(AoD);
we  = pi_sin_AoD.*sin(ZoD);
wa  = pi_sin_AoD.*cos(ZoD);
Rx_streering_vector = exp(1j*we*(0:1:N-1)).';
Tx_streering_vector = exp(1j*wa*(0:1:N-1)).';
H   = Rx_streering_vector*diag(gama)*Tx_streering_vector.';

% both are same
% H_imp=zeros(N,N);
% for i=1:Kc
%     H_imp = H_imp+gama(i)*Rx_streering_vector(:,i)*Tx_streering_vector(:,i).';
% end

pos = (0:1:N-1);
steering_vec = steervec(pos,[AoD,ZoD]');
% H_m = 