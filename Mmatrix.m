N=5;
x=.4;
y=0.45;
lambda = .25;
M = zeros(N,N);
R = zeros(N,N);
v=[0,1,1,2,2,3,3,4]*(N)+[2,1,3,2,4,3,5,4];
M(v)=[x,x,x*y,x*y,x*y,x*y,x,x];
R([1,25])=1;
Z = M-1j*R;
A=lambda*eye(N)+Z;
invZ=inv(Z);
invA=inv(A);