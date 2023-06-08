% Apply proposed method to planar 4-DOF CDM-PRJ
function [q,exitflag,t,k]=TPGI_P4D(Td,q0,para)
t_fabrik=tic;
exitflag=0;
max_time=50e-3;
k_max=1000;
eps_po=1e-6;
overrun_threshold=2;

qlim=para.qlim;
pt=Td(1:3,4);
et=Td(1:3,1);   % The end direction in the plane system is xe!

a=para.a;
a1=a(1); a2=a(2); a3=a(3); a4=a(4); a5=a(5); a6=a(6); a7=a(7);
qslim=qlim;
qslim(2,:)=2*qslim(2,:);
qslim(3,:)=2*qslim(3,:);
qslim(4,:)=2*qslim(4,:);
q=q0;
q1=q(1); q2=q(2); q4=q(3); q6=q(4);
overlimit_flag=zeros(4,k_max);
k=0;
t=999;
N=zeros(3,4);
x0=[1;0;0]; z0=[0;0;1];
l(1)=a1+a2/(2*cos(q2));
l(2)=a2/(2*cos(q2))+a3+a4/(2*cos(q4));
l(3)=a4/(2*cos(q4))+a5+a6/(2*cos(q6));
l(4)=a6/(2*cos(q6))+a7;
e(:,1)=[cos(q1);sin(q1);0];
e(:,2)=[cos(q1+2*q2);sin(q1+2*q2);0];
e(:,3)=[cos(q1+2*q2+2*q4);sin(q1+2*q2+2*q4);0];
e(:,4)=[cos(q1+2*q2+2*q4+2*q6);sin(q1+2*q2+2*q4+2*q6);0];
N(:,1)=e(:,1)*l(1);
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
        q=[q1;q2;q4;q6];
        exitflag=1;
        return
    end
    
    k=k+1;
    
    if k>k_max
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
            q=qlim(:,1)+(qlim(:,2)-qlim(:,1)).*rand(4,1);
            q1=q(1); q2=q(2); q4=q(3); q6=q(4);
            e(:,1)=[cos(q1);sin(q1);0];
            e(:,2)=[cos(q1+2*q2);sin(q1+2*q2);0];
            e(:,3)=[cos(q1+2*q2+2*q4);sin(q1+2*q2+2*q4);0];
            e(:,4)=[cos(q1+2*q2+2*q4+2*q6);sin(q1+2*q2+2*q4+2*q6);0];
            N(:,1)=e(:,1)*l(1);
            for i=2:4
                N(:,i)=N(:,i-1)+e(:,i)*l(i);
            end
        end
    end
    
    %% Forward reaching phase
    N1(:,4)=pt; e(:,4)=et;
    N1(:,3)=N1(:,4)-l(4)*e(:,4);
    e(:,3)=N1(:,3)-N(:,2);
    e(:,3)=e(:,3)/norm(e(:,3));
    qs4=atan2(dot(cross(e(:,3),e(:,4)),z0),dot(e(:,3),e(:,4)));
    if qs4>qslim(4,2)
        overlimit_flag(4,k)=1;
        qs4=qslim(4,2);
        e(:,3)=rot(-z0,qs4,e(:,4));
    elseif qs4<qslim(4,1)
        overlimit_flag(4,k)=1;
        qs4=qslim(4,1);
        e(:,3)=rot(-z0,qs4,e(:,4));
    end
    N1(:,2)=N1(:,3)-l(3)*e(:,3);
    e(:,2)=N1(:,2)-N(:,1);
    e(:,2)=e(:,2)/norm(e(:,2));
    qs3=atan2(dot(cross(e(:,2),e(:,3)),z0),dot(e(:,2),e(:,3)));
    if qs3>qslim(3,2)
        overlimit_flag(3,k)=1;
        qs3=qslim(3,2);
        e(:,2)=rot(-z0,qs3,e(:,3));
    elseif qs3<qslim(3,1)
        overlimit_flag(3,k)=1;
        qs3=qslim(3,1);
        e(:,2)=rot(-z0,qs3,e(:,3));
    end
    N1(:,1)=N1(:,2)-l(2)*e(:,2);
    
    %% Backward reaching phase
    e(:,1)=N1(:,1);
    e(:,1)=e(:,1)/norm(e(:,1));
    qs1=atan2(dot(cross(x0,e(:,1)),z0),dot(x0,e(:,1)));
    if qs1>qslim(1,2)
        overlimit_flag(1,k)=1;
        qs1=qslim(1,2);
        e(:,1)=rot(z0,qs1,x0);
    elseif qs1<qslim(1,1)
        overlimit_flag(1,k)=1;
        qs1=qslim(1,1);
        e(:,1)=rot(z0,qs1,x0);
    end
    N2(:,1)=l(1)*e(:,1);
    e(:,2)=N1(:,2)-N2(:,1);
    e(:,2)=e(:,2)/norm(e(:,2));
    qs2=atan2(dot(cross(e(:,1),e(:,2)),z0),dot(e(:,1),e(:,2)));
    if qs2>qslim(2,2)
        overlimit_flag(2,k)=1;
        qs2=qslim(2,2);
        e(:,2)=rot(z0,qs2,e(:,1));
    elseif qs2<qslim(2,1)
        overlimit_flag(2,k)=1;
        qs2=qslim(2,1);
        e(:,2)=rot(z0,qs2,e(:,1));
    end
    N2(:,2)=N2(:,1)+l(2)*e(:,2);
    e(:,3)=N1(:,3)-N2(:,2);
    e(:,3)=e(:,3)/norm(e(:,3));
    qs3=atan2(dot(cross(e(:,2),e(:,3)),z0),dot(e(:,2),e(:,3)));
    if qs3>qslim(3,2)
        overlimit_flag(3,k)=1;
        qs3=qslim(3,2);
        e(:,3)=rot(z0,qs3,e(:,2));
    elseif qs3<qslim(3,1)
        overlimit_flag(3,k)=1;
        qs3=qslim(3,1);
        e(:,3)=rot(z0,qs3,e(:,2));
    end
    N2(:,3)=N2(:,2)+l(3)*e(:,3);
    e(:,4)=N1(:,4)-N2(:,3);
    e(:,4)=e(:,4)/norm(e(:,4));
    qs4=atan2(dot(cross(e(:,3),e(:,4)),z0),dot(e(:,3),e(:,4)));
    if qs4>qslim(4,2)
        overlimit_flag(4,k)=1;
        qs3=qslim(4,2);
        e(:,4)=rot(z0,qs4,e(:,3));
    elseif qs4<qslim(4,1)
        overlimit_flag(4,k)=1;
        qs4=qslim(4,1);
        e(:,4)=rot(z0,qs4,e(:,3));
    end
    N2(:,4)=N2(:,3)+l(4)*e(:,4);
    
    %% State update phase
    q1=qs1;
    q2=qs2/2;
    q4=qs3/2;
    q6=qs4/2;
    l(1)=a1+a2/(2*cos(q2));
    l(2)=a2/(2*cos(q2))+a3+a4/(2*cos(q4));
    l(3)=a4/(2*cos(q4))+a5+a6/(2*cos(q6));
    l(4)=a6/(2*cos(q6))+a7;
    e(:,1)=[cos(q1);sin(q1);0];
    e(:,2)=[cos(q1+2*q2);sin(q1+2*q2);0];
    e(:,3)=[cos(q1+2*q2+2*q4);sin(q1+2*q2+2*q4);0];
    e(:,4)=[cos(q1+2*q2+2*q4+2*q6);sin(q1+2*q2+2*q4+2*q6);0];
    N(:,1)=e(:,1)*l(1);
    for i=2:4
        N(:,i)=N(:,i-1)+e(:,i)*l(i);
    end
end

%% If fail then return the initial guess
q=q0;