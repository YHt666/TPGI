% Parameters of spatial 6-DOF manipulator.
function para=fun_para_6R
% DH parameters
a=[0,300,80,0,0,90,0,0,0]/1000;
d=[120,0,0,0,220,0,0,0,150]/1000;
alpha=deg2rad([90,0,0,-90,90,0,-90,0,0]);
theta0=deg2rad([-90,90,0,-90,0,90,-90,0,0]);
DH=[theta0;d;a;alpha]';

L(9)=Link();
for i=1:9
    L(i)=Link(DH(i,:),'standard');
    L(i).offset=theta0(i);
end
robot=SerialLink(L,'name','CDM-6R');

Uk=[
    1 0 0 0 0 0
    0 1 0 0 0 0
    0 0 1 0 0 0
    0 0 1 0 0 0
    0 0 0 1 0 0
    0 0 0 0 1 0
    0 0 0 0 1 0
    0 0 0 -1 0 0
    0 0 0 0 0 1
]; 

Ut=[
    1 0 0 0 0 0 0 0 0
    0 1 0 0 0 0 0 0 0
    0 0 1 0 0 0 0 0 0
    0 0 0 0 1 0 0 0 0
    0 0 0 0 0 1 0 0 0
    0 0 0 0 0 0 0 0 1
];

% Joint limits
qlim(1,:)=deg2rad([-180 180]);
qlim(2,:)=deg2rad([-90 90]);
qlim(3,:)=deg2rad([0 80]);
qlim(4,:)=deg2rad([-180 180]);
qlim(5,:)=deg2rad([0 45]);
qlim(6,:)=deg2rad([-180 180]);

para.a=a;
para.d=d;
para.alpha=alpha;
para.theta0=theta0;
para.DH=DH;
para.robot=robot;
para.Uk=Uk;
para.Ut=Ut;
para.qlim=qlim;
