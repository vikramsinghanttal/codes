clear all;
clc;

n=3;%% order of the filter
%f=[linspace(0.5e9,1.5e9,2000)]; %% range of frequency 0.5GHz to 1.5GHz
f =0.5*1e9;
f0=1e8; %% resonant frequency 1GHz
BW=1e8; %% f2-f1
Qu=inf;
lambda= (f0/BW)*((1/Qu)+(f/f0)-(f0/f) );
%% g values
g1=1.5963;
g2=1.0967;
g3=1.5963;
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
    cvx_begin
        variable g(3)
        minimize(abs(1-2*1i*g(1)*(lambda*lambda*g(1)*g(2)-1-1i*lambda*lambda*lambda*g(1)*g(2)*g(3))/(lambda*g(1)*g(2)*(lambda*lambda*g(1)*g(3)+1)-1i*(lambda*lambda*g(1)*g(1)*g(2)-g(1)+g(3)-lambda*lambda*g(1)*g(2)*g(3)))))
    cvx_end