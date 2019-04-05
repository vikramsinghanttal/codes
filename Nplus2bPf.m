clear all;
clc;
n=3;                            %% order of the filter
f=linspace(0.5e9,1.5e9,2000);     %% range of frequency 0.5GHz to 1.5GHz
f0=10^9;                         %% resonant frequency 1GHz
BW=10^8;                         %% f2-f1

%% g values
g1=1.5963;
g2=1.0967;
g3=1.5963;

%% R,I and M matrix
R= [1 0 0 0 0;0 0 0 0 0;0 0 0 0 0;0 0 0 0 0;0 0 0 0 1];
QuV=[inf,1000]; %% quality factor for lossless 
I=eye(n+2);
I(1,1)=0;
I(n+2,n+2)=0;
M=[0 1/sqrt(g1) 0 0 0;1/sqrt(g1) 0 1/sqrt(g1*g2) 0 0;0 1/sqrt(g1*g2) 0 1/sqrt(g1*g2) 0;0 0 1/sqrt(g1*g2) 0 1/sqrt(g3);0 0 0 1/sqrt(g3) 0]; 
%%M=[0 1/sqrt(g1*g2) 0; 1/sqrt(g1*g2) 0 1/sqrt(g2*g3);0 1/sqrt(g2*g3) 0];
%%A=R+  ((1i*2*pi*f0)*U)+1i*M;
S21=zeros(1,2000);
S11=zeros(1,2000);

markers= ['>','<','h','^','p'];
linestyles = [':','-','--','-',':'];
color = ['m','c','y','g','k'];
txt = ['Ideal Waveguide : S11','Ideal Waveguide : S21','Practical Waveguide : S11','Practical Waveguide : S21'];
inc=1;
for Q=1:2
    Qu =QuV(Q);
    for i=1:2000
        %%A=R+  ((1i*2*pi*(f(1,i)/f0))*U)+1i*M;
        lambda= (f0/BW)*  (   (1/Qu)+   (f(1,i)/f0)-     (f0/f(1,i))     );
        %%lambda= (f0/BW)*((f(1,i)/f0)  - (f0/f(1,i))); %%bpf pass
        %%lambda= f0/f(1,i); %% high pass
        A= lambda*I- 1i*R +M;
        S21_X=-2*1i*sqrt(1*1)*inv(A);
        S11_X=1+2*1i*1*inv(A);
        S11(1,i)=S11_X(1,1);
        S21(1,i)=S21_X(n+2,1);
    end
    S21_m=20*log10(abs(S21));
    hold on
    p =plot(f,S21_m);
    p.LineWidth = 1;
    p.Color = color(inc);
    p.LineStyle = linestyles(inc);
    p.Marker    = markers(inc);
    p.MarkerSize = 5;
    S11_m=20*log10(abs(S11));
    hold on
    
    inc = inc+1;
    p =plot(f,S11_m);
    p.LineWidth = 1;
    p.Color = color(inc);
    p.LineStyle = linestyles(inc);
    p.Marker    = markers(inc);
    p.MarkerSize = 2;
end
xlabel('frequency');
legend('Ideal Waveguide : S11','Ideal Waveguide : S21','Practical Waveguide : S11','Practical Waveguide : S21');