% Parameters of planar 4-DOF manipulator.
function para=fun_para_4R
% DH parameters
a=[500,80,400,80,200,80,200]/1000;
d=[0,0,0,0,0,0,0]/1000;
alpha=deg2rad([0,0,0,0,0,0,0]);
theta0=deg2rad([0,0,0,0,0,0,0]);
DH=[theta0;d;a;alpha]';

L(7)=Link();
for i=1:7
    L(i)=Link(DH(i,:),'standard');
    L(i).offset=theta0(i);
end
robot=SerialLink(L,'name','CDM-4R');

Uk=[
    1 0 0 0
    0 1 0 0
    0 1 0 0
    0 0 1 0
    0 0 1 0
    0 0 0 1
    0 0 0 1
]; 

Ut=[
    1 0 0 0 0 0 0
    0 1 0 0 0 0 0
    0 0 0 1 0 0 0
    0 0 0 0 0 1 0
];

% Joint limits
qlim(1,:)=deg2rad([-90 90]);
qlim(2,:)=deg2rad([0 80]);
qlim(3,:)=deg2rad([-80 0]);
qlim(4,:)=deg2rad([0 80]);

para.a=a;
para.d=d;
para.alpha=alpha;
para.theta0=theta0;
para.DH=DH;
para.robot=robot;
para.Uk=Uk;
para.Ut=Ut;
para.qlim=qlim;
