clear all;
clc;
n=3;                            %% order of the filter
f=linspace(0.5e9,1.5e9,2000);   %% range of frequency 0.5GHz to 1.5GHz
f_length = length(f);
f0=10^9;                        %% resonant frequency 1GHz
BW=10^8;                        %% f2-f1

%% g values
g_opt_inf = [1.5963; 1.0967; 1.5963];
g1=1.5963;
g2=1.0967;
g3=1.5963;
g_opt_Qu = zeros(3,1);

g_inf_f = zeros(3,f_length);
g_Qu_f  = zeros(3,f_length);

%% R,I and M matrix
R= [1 0 0 0 0;
    0 0 0 0 0;
    0 0 0 0 0;
    0 0 0 0 0;
    0 0 0 0 1];
QuV=[inf,1000];                 %% quality factor for lossless 
I=eye(n+2);
I(1,1)=0;
I(n+2,n+2)=0;
S21=zeros(f_length,1);
S11_inf_opt =   zeros(f_length,1);
S11_inf_g   =   zeros(f_length,1);
S11_Qu_opt  =   zeros(f_length,1);
S11_Qu_g    =   zeros(f_length,1);
S11_inf_opt_inv=zeros(f_length,1);
mu          =   0.01;
No_of_iterations = 1000;
error       =   zeros(f_length,1);

Qu = 1000;

for fop=1:f_length
    fi = f(fop);
    %%A=R+  ((1i*2*pi*(f(1,i)/f0))*U)+1i*M;
    lambda= (f0/BW)*((1/Qu)+(fi/f0)-(f0/fi));
    M=[0            1/sqrt(g1)      0               0               0         ;
      1/sqrt(g1)   0               1/sqrt(g1*g2)   0               0         ;
      0            1/sqrt(g1*g2)   0               1/sqrt(g1*g2)   0         ;
      0            0               1/sqrt(g1*g2)   0               1/sqrt(g3);
      0            0               0               1/sqrt(g3)      0        ];
       
    A= lambda*I - 1i*R + M;
        
%   S21_X=-2*1i*sqrt(1*1)*inv(A);
    inv_A = inv(A);
    S11_inf_g(fop)= 20*log10(abs(1+2*1i*inv_A(1,1)));
        
    g = [g1;g2;g3];
    N=0;
    D=0;
    for iter = 1:No_of_iterations
        Nr = lambda*g(1)*g(2)-lambda*lambda*lambda*g(1)*g(1)*g(2)*g(3);
        Dr = lambda*g(1)*g(2)+lambda*lambda*lambda*g(1)*g(1)*g(2)*g(3);
        Ni = g(3) + g(1) - lambda*lambda*g(1)*g(1)*g(2) - lambda*lambda*g(1)*g(2)*g(3);
        Di = g(3) - g(1) + lambda*lambda*g(1)*g(1)*g(2) - lambda*lambda*g(1)*g(2)*g(3);

        N  = (Nr*Nr + Ni*Ni);
        D  = (Dr*Dr + Di*Di);

    %     error(iter)     = N/D;

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
    M=[0            1/sqrt(g(1))      0               0               0         ;
      1/sqrt(g(1))   0               1/sqrt(g(1)*g(2))   0               0         ;
      0            1/sqrt(g(1)*g(2))   0               1/sqrt(g(1)*g(2))   0         ;
      0            0               1/sqrt(g(1)*g(2))   0               1/sqrt(g(3));
      0            0               0               1/sqrt(g(3))      0        ];
       
    A= lambda*I - 1i*R + M;
        
%   S21_X=-2*1i*sqrt(1*1)*inv(A);
    inv_A = inv(A);
    fop
    S11_inf_opt_inv(fop)= 20*log10(abs(1+2*1i*inv_A(1,1)));
    
    g_inf_f(:,fop) = g;
    S11_inf_opt(fop) = 20*log10(sqrt(N/D));
end
figure
plot(f,S11_inf_opt,'-r');
hold on
plot(f,S11_inf_opt_inv,'--y');
plot(f,S11_inf_g,':g');
