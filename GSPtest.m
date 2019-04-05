clear all;
clc;

n=3;%% order of the filter
% f=[linspace(0.5e9,1.5e9,2000)]; %% range of frequency 0.5GHz to 1.5GHz
f=1e9;
f0=1e9; %% resonant frequency 1GHz
BW=1e8; %% f2-f1
Qu=inf;
lambda= (f0/BW)*((1/Qu)+(f/f0)-(f0/f) );
%% g values
g = zeros(3,1);
g(1)=1.5963;
g(2)=1.0967;
g(3)=1.5963;
%% R,I and M matrix
%{
R= [1 0 0 0 0;0 0 0 0 0;0 0 0 0 0;0 0 0 0 0;0 0 0 0 1];
Qu=inf; %% quality factor for lossless 
I=eye(n+2);
I(1,1)=0;
I(n+2,n+2)=0;
M=[0            1/sqrt(g1)      0            0              0;
   1/sqrt(g1)   0           1/sqrt(g1*g2)       0           0;
   0            1/sqrt(g1*g2)   0           1/sqrt(g1*g2)   0;
   0            0           1/sqrt(g1*g2)   0      1/sqrt(g3);
   0            0               0           1/sqrt(g3)      0]; 
%%M=[0 1/sqrt(g1*g2) 0; 1/sqrt(g1*g2) 0 1/sqrt(g2*g3);0 1/sqrt(g2*g3) 0];
%%A=R+  ((1i*2*pi*f0)*U)+1i*M;
S21=zeros(1,2000);
S11=zeros(1,2000);
lambda= (f0/BW)*((1/Qu)+(f/f0)-(f0/f) );
A= lambda*I- 1i*R +M;
S21_X=-2*1i*inv(A);
S11_X=1+2*1i*inv(A);
S11=S11_X(1,1);
S21=S21_X(n+2,1);
%}
% S11_bar = 1-2*1i*g1*(lambda*lambda*g1*g2-1-1i*lambda*lambda*lambda*g1*g2*g3)/(lambda*g1*g2*(lambda*lambda*g1*g3+1)-1i*(lambda*lambda*g1*g1*g2-g1+g3-lambda*lambda*g1*g2*g3))
% S11

mu=0.01;
No_of_iterations = 1000;
error = zeros(No_of_iterations,1);
%{
for iter = 1:No_of_iterations
    
    Numerator   = -g(1)*(lambda*lambda*g(1)*g(2)-1-1i*lambda*lambda*lambda*g(1)*g(2)*g(3));
    Denomenator = (lambda*g(1)*g(2)*(lambda*lambda*g(1)*g(3)+1)-1i*(lambda*lambda*g(1)*g(1)*g(2)-g(1)+g(3)-lambda*lambda*g(1)*g(2)*g(3)));
    
    error(iter) = abs(Numerator/Denomenator);
    
    grad_g1     = 2*g(1)*g(2)*lambda*lambda*Denomenator-lambda*g(2)*(2*g(1)*g(3)*lambda*lambda+1)*Numerator-1i*(2*g(1)*g(2)*g(3)*lambda*lambda*lambda*Denomenator+(2*g(1)*g(2)*lambda*lambda-1-lambda*lambda*g(1)*g(3)));
    grad_g2     = g(1)*lambda*lambda*Denomenator-(lambda*lambda*lambda*g(1)*g(1)*g(3)+lambda)*Numerator-1i*(g(1)*g(2)*lambda*lambda*lambda*Denomenator-(g(1)*g(1)*lambda*lambda-lambda*lambda*g(1)*g(3))*Numerator);
    grad_g3     = -Numerator*g(1)*g(1)*g(2)-1i*(g(1)*g(2)*lambda*lambda*lambda*Denomenator-(1-g(1)*g(2)*lambda*lambda)*Numerator);
    
    g = g + mu/(Denomenator*Denomenator)*[grad_g1 ;grad_g2;grad_g3];
end
%}

for iter = 1:No_of_iterations
Nr = lambda*g(1)*g(2)-lambda*lambda*lambda*g(1)*g(1)*g(2)*g(3);
Dr = lambda*g(1)*g(2)+lambda*lambda*lambda*g(1)*g(1)*g(2)*g(3);
Ni = g(3) + g(1) - lambda*lambda*g(1)*g(1)*g(2) - lambda*lambda*g(1)*g(2)*g(3);
Di = g(3) - g(1) + lambda*lambda*g(1)*g(1)*g(2) - lambda*lambda*g(1)*g(2)*g(3);

N  = (Nr*Nr + Ni*Ni);
D  = (Dr*Dr + Di*Di);

error(iter)     = N/D;

del_Nr_g1 = lambda*g(2) - 2*lambda*lambda*lambda*g(1)*g(2)*g(3);
del_Ni_g1 = 1 - 2*lambda*lambda*g(1)*g(2)-lambda*lambda*g(2)*g(3);
del_Dr_g1 = lambda*g(2) + 2*lambda*lambda*lambda*g(1)*g(2)*g(3);
del_Di_g1 = -1 + 2*lambda*lambda*g(1)*g(2)-lambda*lambda*g(2)*g(3);


del_Nr_g2 = lambda*g(1) - lambda*lambda*lambda*g(1)*g(1)*g(3);
del_Ni_g2 = -lambda*lambda*g(1)*g(1)-lambda*lambda*g(1)*g(3);
del_Dr_g2 = lambda*g(1) + lambda*lambda*lambda*g(1)*g(1)*g(3);
del_Di_g2 = lambda*lambda*g(1)*g(1)-lambda*lambda*g(1)*g(3);


del_Nr_g3 = - lambda*lambda*lambda*g(1)*g(1)*g(2);
del_Ni_g3 = 1-lambda*lambda*g(1)*g(2);
del_Dr_g3 = lambda*lambda*lambda*g(1)*g(1)*g(2);
del_Di_g3 = 1-lambda*lambda*g(1)*g(2);

del_mse_g1 = 2*((Nr*del_Nr_g1+Ni*del_Ni_g1)*D-N*(Dr*del_Dr_g1+Di*del_Di_g1))/(D*D);
del_mse_g2 = 2*((Nr*del_Nr_g2+Ni*del_Ni_g2)*D-N*(Dr*del_Dr_g2+Di*del_Di_g2))/(D*D);
del_mse_g3 = 2*((Nr*del_Nr_g3+Ni*del_Ni_g3)*D-N*(Dr*del_Dr_g3+Di*del_Di_g3))/(D*D);

grad_mse = [del_mse_g1;del_mse_g2;del_mse_g3];

g = g - mu*grad_mse;
end
plot(1:1:No_of_iterations,20*log(sqrt(error)));