% Apply proposed method to spatial 6-DOF CDM-PRJ
function [q,exitflag,t,k]=TPGI_S6D(Td,q0,para)
t_fabrik=tic;
exitflag=0;
max_time=50e-3;
k_max=1000;
eps_po=1e-6;
overrun_threshold=2;

qlim=para.qlim;
pt=Td(1:3,4);
et=Td(1:3,3);

d=para.d; a=para.a;
d1=d(1); a2=a(2); a3=a(3);
d5=d(5); a6=a(6); d9=d(9);
qslim=qlim;
qslim(3,:)=2*qslim(3,:);
qslim(5,:)=2*qslim(5,:);
q=q0;
q1=q(1); q2=q(2); q3=q(3);
q5=q(4); q6=q(5);
overlimit_flag=zeros(5,k_max/2);
k=0;
t=999;
N=zeros(3,4);
N(:,1)=[0;0;d1];
x0=[1;0;0]; y0=[0;1;0]; z0=[0;0;1];
l(1)=d1; l(2)=a2+a3/(2*cos(q3));
l(3)=a3/(2*cos(q3))+d5+a6/(2*cos(q6));
l(4)=a6/(2*cos(q6))+d9;
x2=[(-1).*sin(q1).*sin(q2),cos(q1).*sin(q2),cos(q2)]';
z4=[(-1).*sin(q1).*sin(q2+2.*q3),cos(q1).*sin(q2+2.*q3),cos(q2+2.*q3) ...
    ]';
z7=[(-1).*cos(2.*q6).*sin(q1).*sin(q2+2.*q3)+(-1).*(cos(q2+2.*q3).* ...
    cos(q5).*sin(q1)+cos(q1).*sin(q5)).*sin(2.*q6),(-1).*sin(q1).*sin( ...
    q5).*sin(2.*q6)+cos(q1).*(cos(2.*q6).*sin(q2+2.*q3)+cos(q2+2.*q3) ...
    .*cos(q5).*sin(2.*q6)),cos(q2+2.*q3).*cos(2.*q6)+(-1).*cos(q5).* ...
    sin(q2+2.*q3).*sin(2.*q6)]';
e(:,1)=z0; e(:,2)=x2; e(:,3)=z4; e(:,4)=z7;
for i=2:4
    N(:,i)=N(:,i-1)+e(:,i)*l(i);
end

while toc(t_fabrik)<max_time/2
    %% If converge
    dp=N(:,end)-pt;
    do=acos(dot(e(:,end),et));
    epo=[dp;do];
    t=toc(t_fabrik);
    if norm(epo)<eps_po && t<max_time
        % Solve q9
        x8=[(-1).*cos(q1).*sin(2.*q5).*sin(q6).^2+sin(q1).*(cos(q2+2.*q3).*( ...
            cos(q5).^2.*cos(2.*q6)+sin(q5).^2)+(-1).*cos(q5).*sin(q2+2.*q3).* ...
            sin(2.*q6)),(-1).*sin(q1).*sin(2.*q5).*sin(q6).^2+cos(q1).*((-1).* ...
            cos(q2+2.*q3).*(cos(q5).^2.*cos(2.*q6)+sin(q5).^2)+cos(q5).*sin( ...
            q2+2.*q3).*sin(2.*q6)),sin(q2+2.*q3).*(cos(q5).^2.*cos(2.*q6)+sin( ...
            q5).^2)+cos(q2+2.*q3).*cos(q5).*sin(2.*q6)]';
        x9=Td(1:3,1);
        z8=[(-1).*cos(2.*q6).*sin(q1).*sin(q2+2.*q3)+(-1).*(cos(q2+2.*q3).* ...
            cos(q5).*sin(q1)+cos(q1).*sin(q5)).*sin(2.*q6),(-1).*sin(q1).*sin( ...
            q5).*sin(2.*q6)+cos(q1).*(cos(2.*q6).*sin(q2+2.*q3)+cos(q2+2.*q3) ...
            .*cos(q5).*sin(2.*q6)),cos(q2+2.*q3).*cos(2.*q6)+(-1).*cos(q5).* ...
            sin(q2+2.*q3).*sin(2.*q6)]';
        q9=atan2(dot(cross(x8,x9),z8),dot(x8,x9));
        q=[q1;q2;q3;q5;q6;q9];
        exitflag=1;
        return
    end
    
    k=k+1;
    
    if k>k_max/2
        k=k-1;
        break
    end
    
    %% Random disturbance measure
    if k>overrun_threshold
        temp=overlimit_flag(:,k-1);
        if overrun_threshold>=2
            for i=2:overrun_threshold
                temp=temp & overlimit_flag(:,k-i);
            end
        end
        deadlock=any(temp);
        if deadlock
            q=qlim(:,1)+(qlim(:,2)-qlim(:,1)).*rand(6,1);
            q1=q(1); q2=q(2); q3=q(3); q5=q(4); q6=q(5);
            x2=[(-1).*sin(q1).*sin(q2),cos(q1).*sin(q2),cos(q2)]';
            z4=[(-1).*sin(q1).*sin(q2+2.*q3),cos(q1).*sin(q2+2.*q3),cos(q2+2.*q3) ...
                ]';
            z7=[(-1).*cos(2.*q6).*sin(q1).*sin(q2+2.*q3)+(-1).*(cos(q2+2.*q3).* ...
                cos(q5).*sin(q1)+cos(q1).*sin(q5)).*sin(2.*q6),(-1).*sin(q1).*sin( ...
                q5).*sin(2.*q6)+cos(q1).*(cos(2.*q6).*sin(q2+2.*q3)+cos(q2+2.*q3) ...
                .*cos(q5).*sin(2.*q6)),cos(q2+2.*q3).*cos(2.*q6)+(-1).*cos(q5).* ...
                sin(q2+2.*q3).*sin(2.*q6)]';
            e(:,2)=x2; e(:,3)=z4; e(:,4)=z7;
            for i=2:4
                N(:,i)=N(:,i-1)+e(:,i)*l(i);
            end
        end
    end
    
    %% Forward reaching phase
    N1(:,4)=pt; e(:,4)=et;
    N1(:,3)=N1(:,4)-l(4)*e(:,4);
    np=cross(N(:,1),N1(:,3));
    np=np/norm(np);
    Nh2=N(:,2)-dot(N(:,2),np)*np;
    e(:,3)=N1(:,3)-Nh2;
    e(:,3)=e(:,3)/norm(e(:,3));
    qs5=acos(dot(e(:,3),e(:,4)));
    if qs5>qslim(5,2)
        overlimit_flag(5,k)=1;
        qs5=qslim(5,2);
        vn3=cross(e(:,4),e(:,3));
        e(:,3)=rot(vn3,qs5,e(:,4));
    end
    N1(:,2)=N1(:,3)-l(3)*e(:,3);
    
    %% Backward reaching phase
    N2(:,1)=N(:,1);
    e(:,2)=N1(:,2)-N2(:,1);
    e(:,2)=e(:,2)/norm(e(:,2));
    qs2=atan2(dot(cross(e(:,1),e(:,2)),np),dot(e(:,1),e(:,2)));
    if qs2>qslim(2,2)
        overlimit_flag(2,k)=1;
        qs2=qslim(2,2);
        vn1=cross(e(:,1),e(:,2));
        e(:,2)=rot(vn1,qs2,e(:,1));
    elseif qs2<qslim(2,1)
        overlimit_flag(2,k)=1;
        qs2=qslim(2,1);
        vn1=cross(e(:,1),e(:,2));
        e(:,2)=rot(vn1,qs2,e(:,1));
    end
    N2(:,2)=N2(:,1)+l(2)*e(:,2);
    e(:,3)=N1(:,3)-N2(:,2);
    e(:,3)=e(:,3)/norm(e(:,3));
    qs3=atan2(dot(cross(e(:,2),e(:,3)),np),dot(e(:,2),e(:,3)));
    if qs3>qslim(3,2)
        overlimit_flag(3,k)=1;
        qs3=qslim(3,2);
        vn2=cross(e(:,2),e(:,3));
        e(:,3)=rot(vn2,qs3,e(:,2));
    elseif qs3<qslim(3,1)
        overlimit_flag(3,k)=1;
        qs3=qslim(3,1);
        vn2=cross(e(:,2),e(:,3));
        e(:,3)=rot(vn2,qs3,e(:,2));
    end
    N2(:,3)=N2(:,2)+l(3)*e(:,3);
    e(:,4)=N1(:,4)-N2(:,3);
    e(:,4)=e(:,4)/norm(e(:,4));
    qs5=acos(dot(e(:,3),e(:,4)));
    if qs5>qslim(5,2)
        overlimit_flag(5,k)=1;
        qs5=qslim(5,2);
        vn3=cross(e(:,3),e(:,4));
        e(:,4)=rot(vn3,qs5,e(:,3));
    end
    N2(:,4)=N2(:,3)+l(4)*e(:,4);
    
    %% State update phase
    q1=atan2(dot(-y0,np),dot(-x0,np));
    q2=qs2; q3=qs3/2;
    x4=[cos(q2+2.*q3).*sin(q1),(-1).*cos(q1).*cos(q2+2.*q3),sin(q2+2.*q3) ...
        ]';
    y4=[cos(q1),sin(q1),0]';
    q5=atan2(dot(-y4,e(:,4)),dot(-x4,e(:,4)));
    q6=qs5/2;
    
    l(1)=d1; l(2)=a2+a3/(2*cos(q3));
    l(3)=a3/(2*cos(q3))+d5+a6/(2*cos(q6));
    l(4)=a6/(2*cos(q6))+d9;
    x2=[(-1).*sin(q1).*sin(q2),cos(q1).*sin(q2),cos(q2)]';
    z4=[(-1).*sin(q1).*sin(q2+2.*q3),cos(q1).*sin(q2+2.*q3),cos(q2+2.*q3) ...
        ]';
    z7=[(-1).*cos(2.*q6).*sin(q1).*sin(q2+2.*q3)+(-1).*(cos(q2+2.*q3).* ...
        cos(q5).*sin(q1)+cos(q1).*sin(q5)).*sin(2.*q6),(-1).*sin(q1).*sin( ...
        q5).*sin(2.*q6)+cos(q1).*(cos(2.*q6).*sin(q2+2.*q3)+cos(q2+2.*q3) ...
        .*cos(q5).*sin(2.*q6)),cos(q2+2.*q3).*cos(2.*q6)+(-1).*cos(q5).* ...
        sin(q2+2.*q3).*sin(2.*q6)]';
    e(:,1)=z0; e(:,2)=x2; e(:,3)=z4; e(:,4)=z7;
    for i=2:4
        N(:,i)=N(:,i-1)+e(:,i)*l(i);
    end
end

k1_tol=k;

%% Branch change measure
q=q0;
q1=q(1); q2=q(2); q3=q(3);
q5=q(4); q6=q(5);
overlimit_flag=zeros(5,k_max/2);
k=0;
N=zeros(3,4);
N(:,1)=[0;0;d1];
x0=[1;0;0]; y0=[0;1;0]; z0=[0;0;1];
l(1)=d1; l(2)=a2+a3/(2*cos(q3));
l(3)=a3/(2*cos(q3))+d5+a6/(2*cos(q6));
l(4)=a6/(2*cos(q6))+d9;
x2=[(-1).*sin(q1).*sin(q2),cos(q1).*sin(q2),cos(q2)]';
z4=[(-1).*sin(q1).*sin(q2+2.*q3),cos(q1).*sin(q2+2.*q3),cos(q2+2.*q3) ...
    ]';
z7=[(-1).*cos(2.*q6).*sin(q1).*sin(q2+2.*q3)+(-1).*(cos(q2+2.*q3).* ...
    cos(q5).*sin(q1)+cos(q1).*sin(q5)).*sin(2.*q6),(-1).*sin(q1).*sin( ...
    q5).*sin(2.*q6)+cos(q1).*(cos(2.*q6).*sin(q2+2.*q3)+cos(q2+2.*q3) ...
    .*cos(q5).*sin(2.*q6)),cos(q2+2.*q3).*cos(2.*q6)+(-1).*cos(q5).* ...
    sin(q2+2.*q3).*sin(2.*q6)]';
e(:,1)=z0; e(:,2)=x2; e(:,3)=z4; e(:,4)=z7;
for i=2:4
    N(:,i)=N(:,i-1)+e(:,i)*l(i);
end

while toc(t_fabrik)<max_time
    %% If converge
    dp=N(:,end)-pt;
    do=acos(dot(e(:,end),et));
    epo=[dp;do];
    t=toc(t_fabrik);
    if norm(epo)<eps_po && t<max_time
        % Solve q9
        x8=[(-1).*cos(q1).*sin(2.*q5).*sin(q6).^2+sin(q1).*(cos(q2+2.*q3).*( ...
            cos(q5).^2.*cos(2.*q6)+sin(q5).^2)+(-1).*cos(q5).*sin(q2+2.*q3).* ...
            sin(2.*q6)),(-1).*sin(q1).*sin(2.*q5).*sin(q6).^2+cos(q1).*((-1).* ...
            cos(q2+2.*q3).*(cos(q5).^2.*cos(2.*q6)+sin(q5).^2)+cos(q5).*sin( ...
            q2+2.*q3).*sin(2.*q6)),sin(q2+2.*q3).*(cos(q5).^2.*cos(2.*q6)+sin( ...
            q5).^2)+cos(q2+2.*q3).*cos(q5).*sin(2.*q6)]';
        x9=Td(1:3,1);
        z8=[(-1).*cos(2.*q6).*sin(q1).*sin(q2+2.*q3)+(-1).*(cos(q2+2.*q3).* ...
            cos(q5).*sin(q1)+cos(q1).*sin(q5)).*sin(2.*q6),(-1).*sin(q1).*sin( ...
            q5).*sin(2.*q6)+cos(q1).*(cos(2.*q6).*sin(q2+2.*q3)+cos(q2+2.*q3) ...
            .*cos(q5).*sin(2.*q6)),cos(q2+2.*q3).*cos(2.*q6)+(-1).*cos(q5).* ...
            sin(q2+2.*q3).*sin(2.*q6)]';
        q9=atan2(dot(cross(x8,x9),z8),dot(x8,x9));
        q=[q1;q2;q3;q5;q6;q9];
        exitflag=1;
        k=k+k1_tol;
        return
    end
    
    k=k+1;
    
    if k>k_max/2
        k=k+k1_tol;
        k=k-1;
        break
    end
    
    %% Random disturbance measure
    if k>overrun_threshold
        temp=overlimit_flag(:,k-1);
        if overrun_threshold>=2
            for i=2:overrun_threshold
                temp=temp & overlimit_flag(:,k-i);
            end
        end
        deadlock=any(temp);
        if deadlock
            q=qlim(:,1)+(qlim(:,2)-qlim(:,1)).*rand(6,1);
            q1=q(1); q2=q(2); q3=q(3); q5=q(4); q6=q(5);
            x2=[(-1).*sin(q1).*sin(q2),cos(q1).*sin(q2),cos(q2)]';
            z4=[(-1).*sin(q1).*sin(q2+2.*q3),cos(q1).*sin(q2+2.*q3),cos(q2+2.*q3) ...
                ]';
            z7=[(-1).*cos(2.*q6).*sin(q1).*sin(q2+2.*q3)+(-1).*(cos(q2+2.*q3).* ...
                cos(q5).*sin(q1)+cos(q1).*sin(q5)).*sin(2.*q6),(-1).*sin(q1).*sin( ...
                q5).*sin(2.*q6)+cos(q1).*(cos(2.*q6).*sin(q2+2.*q3)+cos(q2+2.*q3) ...
                .*cos(q5).*sin(2.*q6)),cos(q2+2.*q3).*cos(2.*q6)+(-1).*cos(q5).* ...
                sin(q2+2.*q3).*sin(2.*q6)]';
            e(:,2)=x2; e(:,3)=z4; e(:,4)=z7;
            for i=2:4
                N(:,i)=N(:,i-1)+e(:,i)*l(i);
            end
        end
    end
    
    %% Forward reaching phase
    N1(:,4)=pt; e(:,4)=et;
    N1(:,3)=N1(:,4)-l(4)*e(:,4);
    % Reverse np
    np=-cross(N(:,1),N1(:,3));
    np=np/norm(np);
    Nh2=N(:,2)-dot(N(:,2),np)*np;
    e(:,3)=N1(:,3)-Nh2;
    e(:,3)=e(:,3)/norm(e(:,3));
    qs5=acos(dot(e(:,3),e(:,4)));
    if qs5>qslim(5,2)
        overlimit_flag(5,k)=1;
        qs5=qslim(5,2);
        vn3=cross(e(:,4),e(:,3));
        e(:,3)=rot(vn3,qs5,e(:,4));
    end
    N1(:,2)=N1(:,3)-l(3)*e(:,3);
    
    %% Backward reaching phase
    N2(:,1)=N(:,1);
    e(:,2)=N1(:,2)-N2(:,1);
    e(:,2)=e(:,2)/norm(e(:,2));
    qs2=atan2(dot(cross(e(:,1),e(:,2)),np),dot(e(:,1),e(:,2)));
    if qs2>qslim(2,2)
        overlimit_flag(2,k)=1;
        qs2=qslim(2,2);
        vn1=cross(e(:,1),e(:,2));
        e(:,2)=rot(vn1,qs2,e(:,1));
    elseif qs2<qslim(2,1)
        overlimit_flag(2,k)=1;
        qs2=qslim(2,1);
        vn1=cross(e(:,1),e(:,2));
        e(:,2)=rot(vn1,qs2,e(:,1));
    end
    N2(:,2)=N2(:,1)+l(2)*e(:,2);
    e(:,3)=N1(:,3)-N2(:,2);
    e(:,3)=e(:,3)/norm(e(:,3));
    qs3=atan2(dot(cross(e(:,2),e(:,3)),np),dot(e(:,2),e(:,3)));
    if qs3>qslim(3,2)
        overlimit_flag(3,k)=1;
        qs3=qslim(3,2);
        vn2=cross(e(:,2),e(:,3));
        e(:,3)=rot(vn2,qs3,e(:,2));
    elseif qs3<qslim(3,1)
        overlimit_flag(3,k)=1;
        qs3=qslim(3,1);
        vn2=cross(e(:,2),e(:,3));
        e(:,3)=rot(vn2,qs3,e(:,2));
    end
    N2(:,3)=N2(:,2)+l(3)*e(:,3);
    e(:,4)=N1(:,4)-N2(:,3);
    e(:,4)=e(:,4)/norm(e(:,4));
    qs5=acos(dot(e(:,3),e(:,4)));
    if qs5>qslim(5,2)
        overlimit_flag(5,k)=1;
        qs5=qslim(5,2);
        vn3=cross(e(:,3),e(:,4));
        e(:,4)=rot(vn3,qs5,e(:,3));
    end
    N2(:,4)=N2(:,3)+l(4)*e(:,4);
    
    %% State update phase
    q1=atan2(dot(-y0,np),dot(-x0,np));
    q2=qs2; q3=qs3/2;
    x4=[cos(q2+2.*q3).*sin(q1),(-1).*cos(q1).*cos(q2+2.*q3),sin(q2+2.*q3) ...
        ]';
    y4=[cos(q1),sin(q1),0]';
    q5=atan2(dot(-y4,e(:,4)),dot(-x4,e(:,4)));
    q6=qs5/2;
    
    l(1)=d1; l(2)=a2+a3/(2*cos(q3));
    l(3)=a3/(2*cos(q3))+d5+a6/(2*cos(q6));
    l(4)=a6/(2*cos(q6))+d9;
    x2=[(-1).*sin(q1).*sin(q2),cos(q1).*sin(q2),cos(q2)]';
    z4=[(-1).*sin(q1).*sin(q2+2.*q3),cos(q1).*sin(q2+2.*q3),cos(q2+2.*q3) ...
        ]';
    z7=[(-1).*cos(2.*q6).*sin(q1).*sin(q2+2.*q3)+(-1).*(cos(q2+2.*q3).* ...
        cos(q5).*sin(q1)+cos(q1).*sin(q5)).*sin(2.*q6),(-1).*sin(q1).*sin( ...
        q5).*sin(2.*q6)+cos(q1).*(cos(2.*q6).*sin(q2+2.*q3)+cos(q2+2.*q3) ...
        .*cos(q5).*sin(2.*q6)),cos(q2+2.*q3).*cos(2.*q6)+(-1).*cos(q5).* ...
        sin(q2+2.*q3).*sin(2.*q6)]';
    e(:,1)=z0; e(:,2)=x2; e(:,3)=z4; e(:,4)=z7;
    for i=2:4
        N(:,i)=N(:,i-1)+e(:,i)*l(i);
    end
end

%% If fail then return the initial guess
q=q0;