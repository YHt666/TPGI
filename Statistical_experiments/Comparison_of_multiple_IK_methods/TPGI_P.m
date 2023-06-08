function [q,exitflag,t,k]=TPGI_P(Td,q0,para)
% Plain version of TPGI without random disturbance measure.
% Input
% Td: desired end pose;
% q0: initial guess;
% para: parameters of manipulator.
% Output
% q: IK solution;
% exitflag: exit flag whose value is 0 (failure) or 1 (success);
% t: runtime;
% k: number of iterations.

t_fabrik=tic;
exitflag=0;
max_time=50e-3;     % Max runtime
k_max=1000;     % Max number of iterations
eps_po=1e-6;    % Error threshold
step_tol=1e-7;  % Threshold of joint angle variation

qlim=para.qlim;
pt=Td(1:3,4);
et=Td(1:3,3);

d=para.d; a=para.a;
d1=d(1); a2=a(2); a3=a(3); a4=a(4);
a5=a(5); d7=d(7); a8=a(8); d11=d(11);
qslim=qlim;
qslim(3,:)=2*qslim(3,:);
qslim(4,:)=2*qslim(4,:);
qslim(6,:)=2*qslim(6,:);
q=q0(1:6);
q1=q(1); q2=q(2); q3=q(3);
q5=q(4); q7=q(5); q8=q(6);

k=0;
t=999;
N=zeros(3,5);
N(:,1)=[0;0;d1];
x0=[1;0;0]; y0=[0;1;0]; z0=[0;0;1];
l(1)=d1; l(2)=a2+a3/(2*cos(q3));
l(3)=a4+a3/(2*cos(q3))+a5/(2*cos(q5));
l(4)=d7+a5/(2*cos(q5))+a8/(2*cos(q8));
l(5)=d11+a8/(2*cos(q8));
x2=[(-1).*cos(q1).*sin(q2),(-1).*sin(q1).*sin(q2),cos(q2)]';
x4=[(-1).*cos(q1).*sin(q2+2.*q3),(-1).*sin(q1).*sin(q2+2.*q3),cos(q2+ ...
    2.*q3)]';
z6=[(-1).*cos(q1).*sin(q2+2.*(q3+q5)),(-1).*sin(q1).*sin(q2+2.*(q3+ ...
    q5)),cos(q2+2.*(q3+q5))]';
z10=[sin(q1).*sin(q7).*sin(2.*q8)+(-1).*cos(q1).*(cos(2.*q8).*sin(q2+ ...
    2.*(q3+q5))+cos(q2+2.*(q3+q5)).*cos(q7).*sin(2.*q8)),(-1).*cos(2.* ...
    q8).*sin(q1).*sin(q2+2.*(q3+q5))+(-1).*(cos(q2+2.*(q3+q5)).*cos( ...
    q7).*sin(q1)+cos(q1).*sin(q7)).*sin(2.*q8),cos(q2+2.*(q3+q5)).* ...
    cos(2.*q8)+(-1).*cos(q7).*sin(q2+2.*(q3+q5)).*sin(2.*q8)]';
e(:,1)=z0; e(:,2)=x2; e(:,3)=x4;
e(:,4)=z6; e(:,5)=z10;
for i=2:5
    N(:,i)=N(:,i-1)+e(:,i)*l(i);
end

while toc(t_fabrik)<max_time/2
    %% If converge
    dp=N(:,5)-pt;
    do=acos(dot(e(:,5),et));
    epo=[dp;do];
    t=toc(t_fabrik);
    if norm(epo)<eps_po && t<max_time
        % Solve q11
        x11=Td(1:3,1);
        x10=[(1/4).*((-2).*cos(q1+(-1).*q7).*cos(q7).*cos(2.*q8)+(-1).*((-2)+ ...
            cos(2.*q8)).*sin(q1).*sin(2.*q7)+2.*cos(q1).*((1+2.*cos(q2+2.*(q3+ ...
            q5))).*cos(q7).^2.*cos(2.*q8)+2.*cos(q2+2.*(q3+q5)).*sin(q7).^2+( ...
            -2).*cos(q7).*sin(q2+2.*(q3+q5)).*sin(2.*q8))),(-1).*cos(q1).*sin( ...
            2.*q7).*sin(q8).^2+sin(q1).*(cos(q2+2.*(q3+q5)).*(cos(q7).^2.*cos( ...
            2.*q8)+sin(q7).^2)+(-1).*cos(q7).*sin(q2+2.*(q3+q5)).*sin(2.*q8)), ...
            sin(q2+2.*(q3+q5)).*(cos(q7).^2.*cos(2.*q8)+sin(q7).^2)+cos(q2+2.* ...
            (q3+q5)).*cos(q7).*sin(2.*q8)]';
        q11=atan2(dot(cross(x10,x11),z10),dot(x10,x11));
        q=[q1;q2;q3;q5;q7;q8;q11];
        exitflag=1;
        return
    end
    
    k=k+1;
    
    if k>k_max
        k=k-1;
        break
    end
    
    %% Forward reaching phase
    N1(:,5)=pt; e(:,5)=et;
    N1(:,4)=N1(:,5)-l(5)*e(:,5);
    np=cross(N(:,1),N1(:,4));
    np=np/norm(np);
    Nh3=N(:,3)-dot(N(:,3),np)*np;
    e(:,4)=N1(:,4)-Nh3;
    e(:,4)=e(:,4)/norm(e(:,4));
    qs6=acos(dot(e(:,4),e(:,5)));
    if qs6>qslim(6,2)
        qs6=qslim(6,2);
        vn4=cross(e(:,5),e(:,4));
        e(:,4)=rot(vn4,qs6,e(:,5));
    end
    N1(:,3)=N1(:,4)-l(4)*e(:,4);
    Nh2=N(:,2)-dot(N(:,2),np)*np;
    e(:,3)=N1(:,3)-Nh2;
    e(:,3)=e(:,3)/norm(e(:,3));
    qs4=atan2(dot(cross(e(:,3),e(:,4)),np),dot(e(:,3),e(:,4)));
    if qs4>qslim(4,2)
        qs4=qslim(4,2);
        vn3=cross(e(:,4),e(:,3));
        e(:,3)=rot(vn3,qs4,e(:,4));
    elseif qs4<qslim(4,1)
        qs4=qslim(4,1);
        vn3=cross(e(:,4),e(:,3));
        e(:,3)=rot(vn3,qs4,e(:,4));
    end
    N1(:,2)=N1(:,3)-l(3)*e(:,3);
    
    %% Backward reaching phase
    N2(:,1)=N(:,1);
    e(:,2)=N1(:,2)-N2(:,1);
    e(:,2)=e(:,2)/norm(e(:,2));
    qs2=atan2(dot(cross(e(:,1),e(:,2)),np),dot(e(:,1),e(:,2)));
    if qs2>qslim(2,2)
        qs2=qslim(2,2);
        vn1=cross(e(:,1),e(:,2));
        e(:,2)=rot(vn1,qs2,e(:,1));
    elseif qs2<qslim(2,1)
        qs2=qslim(2,1);
        vn1=cross(e(:,1),e(:,2));
        e(:,2)=rot(vn1,qs2,e(:,1));
    end
    N2(:,2)=N2(:,1)+l(2)*e(:,2);
    e(:,3)=N1(:,3)-N2(:,2);
    e(:,3)=e(:,3)/norm(e(:,3));
    qs3=atan2(dot(cross(e(:,2),e(:,3)),np),dot(e(:,2),e(:,3)));
    if qs3>qslim(3,2)
        qs3=qslim(3,2);
        vn2=cross(e(:,2),e(:,3));
        e(:,3)=rot(vn2,qs3,e(:,2));
    elseif qs3<qslim(3,1)
        qs3=qslim(3,1);
        vn2=cross(e(:,2),e(:,3));
        e(:,3)=rot(vn2,qs3,e(:,2));
    end
    N2(:,3)=N2(:,2)+l(3)*e(:,3);
    e(:,4)=N1(:,4)-N2(:,3);
    e(:,4)=e(:,4)/norm(e(:,4));
    qs4=atan2(dot(cross(e(:,3),e(:,4)),np),dot(e(:,3),e(:,4)));
    if qs4>qslim(4,2)
        qs4=qslim(4,2);
        vn3=cross(e(:,3),e(:,4));
        e(:,4)=rot(vn3,qs4,e(:,3));
    elseif qs4<qslim(4,1)
        qs4=qslim(4,1);
        vn3=cross(e(:,3),e(:,4));
        e(:,4)=rot(vn3,qs4,e(:,3));
    end
    N2(:,4)=N2(:,3)+l(4)*e(:,4);
    e(:,5)=N1(:,5)-N2(:,4);
    e(:,5)=e(:,5)/norm(e(:,5));
    qs6=acos(dot(e(:,4),e(:,5)));
    if qs6>qslim(6,2)
        qs6=qslim(6,2);
        vn4=cross(e(:,4),e(:,5));
        e(:,5)=rot(vn4,qs6,e(:,4));
    end
    N2(:,5)=N2(:,4)+l(5)*e(:,5);
    
    %% State update phase
    q1=atan2(dot(x0,np),dot(-y0,np));
    q2=qs2; q3=qs3/2; q5=qs4/2;
    x6=[cos(q1).*cos(q2+2.*(q3+q5)),cos(q2+2.*(q3+q5)).*sin(q1),sin(q2+ ...
        2.*(q3+q5))]';
    y6=[(-1).*sin(q1),cos(q1),0]';
    q7=atan2(dot(-y6,e(:,5)),dot(-x6,e(:,5)));
    q8=qs6/2;
    
    q_last=q;
    q=[q1;q2;q3;q5;q7;q8];
    if norm(q-q_last)<step_tol
        break
    end
    
    l(1)=d1; l(2)=a2+a3/(2*cos(q3));
    l(3)=a4+a3/(2*cos(q3))+a5/(2*cos(q5));
    l(4)=d7+a5/(2*cos(q5))+a8/(2*cos(q8));
    l(5)=d11+a8/(2*cos(q8));
    x2=[(-1).*cos(q1).*sin(q2),(-1).*sin(q1).*sin(q2),cos(q2)]';
    x4=[(-1).*cos(q1).*sin(q2+2.*q3),(-1).*sin(q1).*sin(q2+2.*q3),cos(q2+ ...
        2.*q3)]';
    z6=[(-1).*cos(q1).*sin(q2+2.*(q3+q5)),(-1).*sin(q1).*sin(q2+2.*(q3+ ...
        q5)),cos(q2+2.*(q3+q5))]';
    z10=[sin(q1).*sin(q7).*sin(2.*q8)+(-1).*cos(q1).*(cos(2.*q8).*sin(q2+ ...
        2.*(q3+q5))+cos(q2+2.*(q3+q5)).*cos(q7).*sin(2.*q8)),(-1).*cos(2.* ...
        q8).*sin(q1).*sin(q2+2.*(q3+q5))+(-1).*(cos(q2+2.*(q3+q5)).*cos( ...
        q7).*sin(q1)+cos(q1).*sin(q7)).*sin(2.*q8),cos(q2+2.*(q3+q5)).* ...
        cos(2.*q8)+(-1).*cos(q7).*sin(q2+2.*(q3+q5)).*sin(2.*q8)]';
    e(:,1)=z0; e(:,2)=x2; e(:,3)=x4;
    e(:,4)=z6; e(:,5)=z10;
    for i=2:5
        N(:,i)=N(:,i-1)+e(:,i)*l(i);
    end
end

k1_tol=k;

%% Branch change measure
q=q0(1:6);
q1=q(1); q2=q(2); q3=q(3);
q5=q(4); q7=q(5); q8=q(6);
k=0;
N=zeros(3,5);
N(:,1)=[0;0;d1];
l(1)=d1; l(2)=a2+a3/(2*cos(q3));
l(3)=a4+a3/(2*cos(q3))+a5/(2*cos(q5));
l(4)=d7+a5/(2*cos(q5))+a8/(2*cos(q8));
l(5)=d11+a8/(2*cos(q8));
x2=[(-1).*cos(q1).*sin(q2),(-1).*sin(q1).*sin(q2),cos(q2)]';
x4=[(-1).*cos(q1).*sin(q2+2.*q3),(-1).*sin(q1).*sin(q2+2.*q3),cos(q2+ ...
    2.*q3)]';
z6=[(-1).*cos(q1).*sin(q2+2.*(q3+q5)),(-1).*sin(q1).*sin(q2+2.*(q3+ ...
    q5)),cos(q2+2.*(q3+q5))]';
z10=[sin(q1).*sin(q7).*sin(2.*q8)+(-1).*cos(q1).*(cos(2.*q8).*sin(q2+ ...
    2.*(q3+q5))+cos(q2+2.*(q3+q5)).*cos(q7).*sin(2.*q8)),(-1).*cos(2.* ...
    q8).*sin(q1).*sin(q2+2.*(q3+q5))+(-1).*(cos(q2+2.*(q3+q5)).*cos( ...
    q7).*sin(q1)+cos(q1).*sin(q7)).*sin(2.*q8),cos(q2+2.*(q3+q5)).* ...
    cos(2.*q8)+(-1).*cos(q7).*sin(q2+2.*(q3+q5)).*sin(2.*q8)]';
e(:,1)=z0; e(:,2)=x2; e(:,3)=x4;
e(:,4)=z6; e(:,5)=z10;
for i=2:5
    N(:,i)=N(:,i-1)+e(:,i)*l(i);
end

while toc(t_fabrik)<max_time
    %% If converge
    dp=N(:,5)-pt;
    do=acos(dot(e(:,5),et));
    epo=[dp;do];
    t=toc(t_fabrik);
    if norm(epo)<eps_po && t<max_time
        % Solve q11
        x11=Td(1:3,1);
        x10=[(1/4).*((-2).*cos(q1+(-1).*q7).*cos(q7).*cos(2.*q8)+(-1).*((-2)+ ...
            cos(2.*q8)).*sin(q1).*sin(2.*q7)+2.*cos(q1).*((1+2.*cos(q2+2.*(q3+ ...
            q5))).*cos(q7).^2.*cos(2.*q8)+2.*cos(q2+2.*(q3+q5)).*sin(q7).^2+( ...
            -2).*cos(q7).*sin(q2+2.*(q3+q5)).*sin(2.*q8))),(-1).*cos(q1).*sin( ...
            2.*q7).*sin(q8).^2+sin(q1).*(cos(q2+2.*(q3+q5)).*(cos(q7).^2.*cos( ...
            2.*q8)+sin(q7).^2)+(-1).*cos(q7).*sin(q2+2.*(q3+q5)).*sin(2.*q8)), ...
            sin(q2+2.*(q3+q5)).*(cos(q7).^2.*cos(2.*q8)+sin(q7).^2)+cos(q2+2.* ...
            (q3+q5)).*cos(q7).*sin(2.*q8)]';
        q11=atan2(dot(cross(x10,x11),z10),dot(x10,x11));
        q=[q1;q2;q3;q5;q7;q8;q11];
        exitflag=1;
        k=k+k1_tol;
        return
    end
    
    k=k+1;
    
    if k>k_max
        k=k-1;
        break
    end
    
    %% Forward reaching phase
    N1(:,5)=pt; e(:,5)=et;
    N1(:,4)=N1(:,5)-l(5)*e(:,5);
    % Normal vector of motion plane (np) in branch 2 is opposite to that in branch 1.
    np=-cross(N(:,1),N1(:,4));
    np=np/norm(np);
    Nh3=N(:,3)-dot(N(:,3),np)*np;
    e(:,4)=N1(:,4)-Nh3;
    e(:,4)=e(:,4)/norm(e(:,4));
    qs6=acos(dot(e(:,4),e(:,5)));
    if qs6>qslim(6,2)
        qs6=qslim(6,2);
        vn4=cross(e(:,5),e(:,4));
        e(:,4)=rot(vn4,qs6,e(:,5));
    end
    N1(:,3)=N1(:,4)-l(4)*e(:,4);
    Nh2=N(:,2)-dot(N(:,2),np)*np;
    e(:,3)=N1(:,3)-Nh2;
    e(:,3)=e(:,3)/norm(e(:,3));
    qs4=atan2(dot(cross(e(:,3),e(:,4)),np),dot(e(:,3),e(:,4)));
    if qs4>qslim(4,2)
        qs4=qslim(4,2);
        vn3=cross(e(:,4),e(:,3));
        e(:,3)=rot(vn3,qs4,e(:,4));
    elseif qs4<qslim(4,1)
        qs4=qslim(4,1);
        vn3=cross(e(:,4),e(:,3));
        e(:,3)=rot(vn3,qs4,e(:,4));
    end
    N1(:,2)=N1(:,3)-l(3)*e(:,3);
    
    %% Backward reaching phase
    N2(:,1)=N(:,1);
    e(:,2)=N1(:,2)-N2(:,1);
    e(:,2)=e(:,2)/norm(e(:,2));
    qs2=atan2(dot(cross(e(:,1),e(:,2)),np),dot(e(:,1),e(:,2)));
    if qs2>qslim(2,2)
        qs2=qslim(2,2);
        vn1=cross(e(:,1),e(:,2));
        e(:,2)=rot(vn1,qs2,e(:,1));
    elseif qs2<qslim(2,1)
        qs2=qslim(2,1);
        vn1=cross(e(:,1),e(:,2));
        e(:,2)=rot(vn1,qs2,e(:,1));
    end
    N2(:,2)=N2(:,1)+l(2)*e(:,2);
    e(:,3)=N1(:,3)-N2(:,2);
    e(:,3)=e(:,3)/norm(e(:,3));
    qs3=atan2(dot(cross(e(:,2),e(:,3)),np),dot(e(:,2),e(:,3)));
    if qs3>qslim(3,2)
        qs3=qslim(3,2);
        vn2=cross(e(:,2),e(:,3));
        e(:,3)=rot(vn2,qs3,e(:,2));
    elseif qs3<qslim(3,1)
        qs3=qslim(3,1);
        vn2=cross(e(:,2),e(:,3));
        e(:,3)=rot(vn2,qs3,e(:,2));
    end
    N2(:,3)=N2(:,2)+l(3)*e(:,3);
    e(:,4)=N1(:,4)-N2(:,3);
    e(:,4)=e(:,4)/norm(e(:,4));
    qs4=atan2(dot(cross(e(:,3),e(:,4)),np),dot(e(:,3),e(:,4)));
    if qs4>qslim(4,2)
        qs4=qslim(4,2);
        vn3=cross(e(:,3),e(:,4));
        e(:,4)=rot(vn3,qs4,e(:,3));
    elseif qs4<qslim(4,1)
        qs4=qslim(4,1);
        vn3=cross(e(:,3),e(:,4));
        e(:,4)=rot(vn3,qs4,e(:,3));
    end
    N2(:,4)=N2(:,3)+l(4)*e(:,4);
    e(:,5)=N1(:,5)-N2(:,4);
    e(:,5)=e(:,5)/norm(e(:,5));
    qs6=acos(dot(e(:,4),e(:,5)));
    if qs6>qslim(6,2)
        qs6=qslim(6,2);
        vn4=cross(e(:,4),e(:,5));
        e(:,5)=rot(vn4,qs6,e(:,4));
    end
    N2(:,5)=N2(:,4)+l(5)*e(:,5);
    
    %% State update phase
    q1=atan2(dot(x0,np),dot(-y0,np));
    q2=qs2; q3=qs3/2; q5=qs4/2;
    x6=[cos(q1).*cos(q2+2.*(q3+q5)),cos(q2+2.*(q3+q5)).*sin(q1),sin(q2+ ...
        2.*(q3+q5))]';
    y6=[(-1).*sin(q1),cos(q1),0]';
    q7=atan2(dot(-y6,e(:,5)),dot(-x6,e(:,5)));
    q8=qs6/2;
    
    q_last=q;
    q=[q1;q2;q3;q5;q7;q8];
    if norm(q-q_last)<step_tol
        break
    end
    
    l(1)=d1; l(2)=a2+a3/(2*cos(q3));
    l(3)=a4+a3/(2*cos(q3))+a5/(2*cos(q5));
    l(4)=d7+a5/(2*cos(q5))+a8/(2*cos(q8));
    l(5)=d11+a8/(2*cos(q8));
    x2=[(-1).*cos(q1).*sin(q2),(-1).*sin(q1).*sin(q2),cos(q2)]';
    x4=[(-1).*cos(q1).*sin(q2+2.*q3),(-1).*sin(q1).*sin(q2+2.*q3),cos(q2+ ...
        2.*q3)]';
    z6=[(-1).*cos(q1).*sin(q2+2.*(q3+q5)),(-1).*sin(q1).*sin(q2+2.*(q3+ ...
        q5)),cos(q2+2.*(q3+q5))]';
    z10=[sin(q1).*sin(q7).*sin(2.*q8)+(-1).*cos(q1).*(cos(2.*q8).*sin(q2+ ...
        2.*(q3+q5))+cos(q2+2.*(q3+q5)).*cos(q7).*sin(2.*q8)),(-1).*cos(2.* ...
        q8).*sin(q1).*sin(q2+2.*(q3+q5))+(-1).*(cos(q2+2.*(q3+q5)).*cos( ...
        q7).*sin(q1)+cos(q1).*sin(q7)).*sin(2.*q8),cos(q2+2.*(q3+q5)).* ...
        cos(2.*q8)+(-1).*cos(q7).*sin(q2+2.*(q3+q5)).*sin(2.*q8)]';
    e(:,1)=z0; e(:,2)=x2; e(:,3)=x4;
    e(:,4)=z6; e(:,5)=z10;
    for i=2:5
        N(:,i)=N(:,i-1)+e(:,i)*l(i);
    end
end

%% If TPGI solver fails then the initial guess is returned.
q=q0;