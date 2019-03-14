%function H = Channel(d,lambda,AoA,AoD,ZoA,ZoD)
d      = 0.005;
lambda = 0.01;
Nt     = 4;
Nr     = 4;
Nc     = 10;
AoA    = randi([-90,90],Nc,1);
AoD    = randi([-90,90],Nc,1);
ZoA    = randi([-90,90],Nc,1);
ZoD    = randi([-90,90],Nc,1);
Tx_streering_vector = steervec(0:d:length(AoD),[AoD ZoD]);
Rx_streering_vector = steervec(0:d:length(AoA),[AoA ZoA]);   

% end