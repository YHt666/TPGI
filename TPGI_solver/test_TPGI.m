close all; clear; clc

para=fun_para;
DH=para.DH;
qlim=para.qlim;
Uk=para.Uk;
robot=para.robot;

a=para.a;
d=para.d;
a2=a(2); a3=a(3); a4=a(4); a5=a(5);
a8=a(8); d1=d(1); d7=d(7); d11=d(11);

n=size(qlim,1);
q_rand=qlim(:,1)+(qlim(:,2)-qlim(:,1)).*rand(n,1);
q0_rand=qlim(:,1)+(qlim(:,2)-qlim(:,1)).*rand(n,1);
rad2deg([q_rand,q0_rand])

Te=robot.fkine(Uk*q_rand).T;
[q,exitflag,t,k]=TPGI(Te,q0_rand,para);

rad2deg(q)
exitflag
t
k
Tc=robot.fkine(Uk*q).T;
delta_T=Tc-Te
