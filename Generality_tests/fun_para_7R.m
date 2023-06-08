% Parameters of spatial 7-DOF manipulator.
function para=fun_para_7R
% DH parameters
a=[0,0,0,80,0,0,100,0,0,0]/1000;
d=[0,0,300,0,0,200,0,0,0,100]/1000;
alpha=deg2rad([90,90,90,0,90,90,0,90,0,0]);
theta0=deg2rad([0,90,90,90,90,180,90,90,0,0]);
DH=[theta0;d;a;alpha]';

L(10)=Link();
for i=1:10
    L(i)=Link(DH(i,:),'standard');
    L(i).offset=theta0(i);
end
robot=SerialLink(L,'name','CDM-7R');

Uk=[
    1 0 0 0 0 0 0
    0 1 0 0 0 0 0
    0 0 1 0 0 0 0
    0 0 0 1 0 0 0
    0 0 0 1 0 0 0
    0 0 0 0 1 0 0
    0 0 0 0 0 1 0
    0 0 0 0 0 1 0
    0 0 0 0 -1 0 0
    0 0 0 0 0 0 1
]; 

Ut=[
    1 0 0 0 0 0 0 0 0 0
    0 1 0 0 0 0 0 0 0 0
    0 0 1 0 0 0 0 0 0 0
    0 0 0 1 0 0 0 0 0 0
    0 0 0 0 0 1 0 0 0 0
    0 0 0 0 0 0 1 0 0 0
    0 0 0 0 0 0 0 0 0 1
];

% Joint limits
qlim(1,:)=deg2rad([-180 180]);
qlim(2,:)=deg2rad([0 90]);
qlim(3,:)=deg2rad([-180 180]);
qlim(4,:)=deg2rad([0 80]);
qlim(5,:)=deg2rad([-180 180]);
qlim(6,:)=deg2rad([0 45]);
qlim(7,:)=deg2rad([-180 180]);

para.a=a;
para.d=d;
para.alpha=alpha;
para.theta0=theta0;
para.DH=DH;
para.robot=robot;
para.Uk=Uk;
para.Ut=Ut;
para.qlim=qlim;
