function para=fun_para
% Parameters of manipulator.
% DH parameters
a=[0 ,  1000 , 80 , 650 , 80 , 0 ,  0 , 100 , 0 , 0 , 0]/1000;
d=[300 , 0 ,   0 ,  0 ,   0 , 0 , 300 , 0 , 0 , 0 , 150]/1000;
alpha=deg2rad([90,0,0,0,0,-90,90,0,-90,0,0]);
theta0=deg2rad([0,90,0,0,0,-90,0,90,-90,0,0]);
DH=[theta0;d;a;alpha]';

L(11)=Link();
for i=1:11
    L(i)=Link(DH(i,:),'standard');
    L(i).offset=theta0(i);
end
robot=SerialLink(L,'name','CDM');

Uk=[1  0   0    0    0    0   0  0   0   0   0
    0  1   0    0    0    0   0  0   0   0   0
    0  0   1    1    0    0   0  0   0   0   0
    0  0   0    0    1    1   0  0   0   0   0
    0  0   0    0    0    0   1  0   0  -1   0
    0  0   0    0    0    0   0  1   1   0   0
    0  0   0    0    0    0   0  0   0   0   1 ]'; 

Ut=[1  0   0    0    0    0   0  0   0   0   0
    0  1   0    0    0    0   0  0   0   0   0
    0  0   1    0    0    0   0  0   0   0   0
    0  0   0    0    1    0   0  0   0   0   0
    0  0   0    0    0    0   1  0   0   0   0
    0  0   0    0    0    0   0  1   0   0   0
    0  0   0    0    0    0   0  0   0   0   1 ];  

% Joint limits
qlim(1,:)=deg2rad([-180 180]); %q1
qlim(2,:)=deg2rad([-90 90]); %q2
qlim(3,:)=deg2rad([0 80]); %q3
qlim(4,:)=deg2rad([-80 0]); %q5
qlim(5,:)=deg2rad([-180 180]); %q7
qlim(6,:)=deg2rad([0 45]); %q8
qlim(7,:)=deg2rad([-180 180]); %q11

para.a=a;
para.d=d;
para.alpha=alpha;
para.theta0=theta0;
para.DH=DH;
para.robot=robot;
para.Uk=Uk;
para.Ut=Ut;
para.qlim=qlim;
