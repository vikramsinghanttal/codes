Nc=10; 
Nt=8;
Ang_departure = -90+rand(2,Nc)*90;
Stv2D=steervec(.5*(0:1:Nc-1),Ang_departure);
w_e=pi*sin(Ang_departure(1,:))*sin(Ang_departure(2,:));
w_a=pi*sin(Ang_departure(1,:))*cos(Ang_departure(2,:));
ind=0:1:Nt-1;
a_N=