function [t1,y1,t2,y2]=osmomodel

%% initial steady process
global glu
initial=[1,50,0.1,0.1,0,0,0,1,0];

glu=0.1;
[t01,y01]=ode23s(@func,[0,1000],initial);
steady1=y01(end,:); index=find(t01>-10); t0=t01(index); y0=y01(index,:);
glu=2;
[t02,y02]=ode23s(@func,[0,1000],initial);steady2=y02(end,:); 
%% stimulus
glu=0.1;
[t1,y1]=ode23s(@func,[0:0.1:6000],steady1);
glu=2;
[t2,y2]=ode23s(@func,[0:0.1:6000],steady2);

function dy=func(t,y)
global glu 
dy=zeros(9,1);
bGB=300;

k1=5;k2=0.25;k3=100;a2=1;a3=1;a4=0.05/50;kp1=100;kp2=100;
a5=1;a5_1=0.01/50;ks=100;  a6=0.1*10; a7=0.1*10; kh=0.05;a7_1=0.001;
a8=0.5; a9=0.05; c0=50; a10=0.1;kh2=0.001;kr=0.2*10;    
kcl=0;deltcmax=10;
a11=0.3;a1=1/1000;
if t>1800
    kcl=80;
end

deltc=max(kcl-y(4)-5,0);deltc2=-min(kcl-y(4),0);
dy(1)=a1*y(2)/(k1+y(2))*(1/(1+(deltc/deltcmax)^2));
dy(2)=a2*glu/(k2+glu)-bGB*a1*y(2)/(k1+y(2))-a3*y(2)/(k3+y(2))*(a4+y(5)/(kp1+y(5)))-dy(1)*y(2); 
dy(3)=2*a3*y(2)/(k3+y(2))*(a4+y(5)/(kp1+y(5)))-a5*y(3)/(ks+y(3))*(a5_1+y(9)/(y(9)+kp2))-dy(1)*y(3);
dy(4)=a5*y(3)/(ks+y(3))*(a5_1+y(9)/(y(9)+kp2))-dy(1)*y(4)-kr*y(8)*deltc2;
dy(5)=a6*y(7)/(y(7)+kh)-dy(1)*y(5);
dy(6)=a7*y(7)/(y(7)+kh)-dy(1)*y(6)-a7_1*y(6);
dy(7)=a8*deltc/(c0+deltc)-(a8*deltc/(c0+deltc)+a9)*y(7);
dy(8)=a10*(kh2/(y(7)+kh2))-a10*y(8);
dy(9)=a11*y(6)-dy(1)*2*y(9);



