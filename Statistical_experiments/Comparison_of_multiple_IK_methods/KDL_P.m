function [q,exitflag,t,j]=KDL_P(Td,qc,para)
% Plain version of KDL without random disturbance measure.
% Input
% Td: desired end pose;
% qc: initial guess;
% para: parameters of manipulator.
% Output
% q: IK solution;
% exitflag: exit flag whose value is 0 (failure) or 1 (success);
% t: runtime;
% j: number of iterations.

t_kdlp=tic;
exitflag=0;
max_time=50e-3;     % Max runtime
j_max=1000;     % Max number of iterations
eps_po=1e-6;    % Error threshold
alpha=1;    % Learning rate
step_tol=1e-7;  % Threshold of joint angle variation

a=para.a;
d=para.d;
a2=a(2);a3=a(3);a4=a(4);a5=a(5);a8=a(8);
d1=d(1);d7=d(7);d11=d(11);
qlim=para.qlim;

q=qc;
overlimit_flag=zeros(7,j_max);
j=1;
while toc(t_kdlp)<max_time
    if j>j_max
        break
    end
    
    q1=qc(1);q2=qc(2);q3=qc(3);q5=qc(4);q7=qc(5);q8=qc(6);q11=qc(7);
    % Forward kinematics
    Tc=[(-1).*sin(q11+(-1).*q7).*(cos(q7).*sin(q1)+cos(q1).*cos(q2+2.*( ...
        q3+q5)).*sin(q7))+cos(q11+(-1).*q7).*(cos(2.*q8).*(cos(q1).*cos( ...
        q2+2.*(q3+q5)).*cos(q7)+(-1).*sin(q1).*sin(q7))+(-1).*cos(q1).* ...
        sin(q2+2.*(q3+q5)).*sin(2.*q8)),(-1).*cos(q11+(-1).*q7).*(cos(q7) ...
        .*sin(q1)+cos(q1).*cos(q2+2.*(q3+q5)).*sin(q7))+(-1).*sin(q11+(-1) ...
        .*q7).*(cos(2.*q8).*(cos(q1).*cos(q2+2.*(q3+q5)).*cos(q7)+(-1).* ...
        sin(q1).*sin(q7))+(-1).*cos(q1).*sin(q2+2.*(q3+q5)).*sin(2.*q8)), ...
        sin(q1).*sin(q7).*sin(2.*q8)+(-1).*cos(q1).*(cos(2.*q8).*sin(q2+ ...
        2.*(q3+q5))+cos(q2+2.*(q3+q5)).*cos(q7).*sin(2.*q8)),(-1).*cos(q1) ...
        .*(a2.*sin(q2)+a3.*sin(q2+q3)+a4.*sin(q2+2.*q3)+a5.*sin(q2+2.*q3+ ...
        q5)+(d7+a8.*cos(q8)+d11.*cos(2.*q8)).*sin(q2+2.*(q3+q5)))+(-1).* ...
        cos(q1).*cos(q2+2.*(q3+q5)).*cos(q7).*(a8+2.*d11.*cos(q8)).*sin( ...
        q8)+(a8+2.*d11.*cos(q8)).*sin(q1).*sin(q7).*sin(q8);(-1).*sin(q11+ ...
        (-1).*q7).*((-1).*cos(q1).*cos(q7)+cos(q2+2.*(q3+q5)).*sin(q1).* ...
        sin(q7))+cos(q11+(-1).*q7).*(cos(2.*q8).*(cos(q2+2.*(q3+q5)).*cos( ...
        q7).*sin(q1)+cos(q1).*sin(q7))+(-1).*sin(q1).*sin(q2+2.*(q3+q5)).* ...
        sin(2.*q8)),cos(q11+(-1).*q7).*(cos(q1).*cos(q7)+(-1).*cos(q2+2.*( ...
        q3+q5)).*sin(q1).*sin(q7))+(-1).*sin(q11+(-1).*q7).*(cos(2.*q8).*( ...
        cos(q2+2.*(q3+q5)).*cos(q7).*sin(q1)+cos(q1).*sin(q7))+(-1).*sin( ...
        q1).*sin(q2+2.*(q3+q5)).*sin(2.*q8)),(-1).*cos(2.*q8).*sin(q1).* ...
        sin(q2+2.*(q3+q5))+(-1).*(cos(q2+2.*(q3+q5)).*cos(q7).*sin(q1)+ ...
        cos(q1).*sin(q7)).*sin(2.*q8),(-1).*sin(q1).*(a2.*sin(q2)+a3.*sin( ...
        q2+q3)+a4.*sin(q2+2.*q3)+a5.*sin(q2+2.*q3+q5)+(d7+a8.*cos(q8)+ ...
        d11.*cos(2.*q8)).*sin(q2+2.*(q3+q5)))+(-1).*cos(q2+2.*(q3+q5)).* ...
        cos(q7).*(a8+2.*d11.*cos(q8)).*sin(q1).*sin(q8)+(-1).*cos(q1).*( ...
        a8+2.*d11.*cos(q8)).*sin(q7).*sin(q8);(-1).*sin(q2+2.*(q3+q5)).* ...
        sin(q11+(-1).*q7).*sin(q7)+cos(q11+(-1).*q7).*(cos(q7).*cos(2.*q8) ...
        .*sin(q2+2.*(q3+q5))+cos(q2+2.*(q3+q5)).*sin(2.*q8)),(-1).*cos( ...
        q11+(-1).*q7).*sin(q2+2.*(q3+q5)).*sin(q7)+(-1).*sin(q11+(-1).*q7) ...
        .*(cos(q7).*cos(2.*q8).*sin(q2+2.*(q3+q5))+cos(q2+2.*(q3+q5)).* ...
        sin(2.*q8)),cos(q2+2.*(q3+q5)).*cos(2.*q8)+(-1).*cos(q7).*sin(q2+ ...
        2.*(q3+q5)).*sin(2.*q8),d1+a2.*cos(q2)+a3.*cos(q2+q3)+a4.*cos(q2+ ...
        2.*q3)+a5.*cos(q2+2.*q3+q5)+cos(q2+2.*(q3+q5)).*(d7+a8.*cos(q8)+ ...
        d11.*cos(2.*q8))+(-1).*cos(q7).*(a8+2.*d11.*cos(q8)).*sin(q2+2.*( ...
        q3+q5)).*sin(q8);0,0,0,1];
    
    % End pose error
    ep=Td(1:3,4)-Tc(1:3,4);
    eo=0.5*(cross(Tc(1:3,1),Td(1:3,1))+cross(Tc(1:3,2),Td(1:3,2))+cross(Tc(1:3,3),Td(1:3,3)));
    epo=[ep;eo];
    
    t=toc(t_kdlp);
    if norm(epo)<eps_po && t<max_time
        q=qc;
        exitflag=1;
        return
    end
    
    % Jacobian matrix
    Js=[sin(q1).*(a2.*sin(q2)+a3.*sin(q2+q3)+a4.*sin(q2+2.*q3)+a5.*sin(q2+2.*q3+q5)+(d7+a8.*cos(q8)+d11.*cos(2.*q8)).*sin(q2+2.*(q3+q5)))+cos(q2+2.*(q3+q5)).*cos(q7).*(a8+2.*d11.*cos(q8)).*sin(q1).*sin(q8)+cos(q1).*(a8+2.*d11.*cos(q8)).*sin(q7).*sin(q8),(-1).*cos(q1).*(a2.*cos(q2)+a3.*cos(q2+q3)+a4.*cos(q2+2.*q3)+a5.*cos(q2+2.*q3+q5)+cos(q2+2.*(q3+q5)).*(d7+a8.*cos(q8)+d11.*cos(2.*q8))+(-1).*cos(q7).*(a8+2.*d11.*cos(q8)).*sin(q2+2.*(q3+q5)).*sin(q8)),(-1).*cos(q1).*(a3.*cos(q2+q3)+2.*(a4.*cos(q2+2.*q3)+a5.*cos(q2+2.*q3+q5)+cos(q2+2.*(q3+q5)).*(d7+a8.*cos(q8)+d11.*cos(2.*q8))+(-1).*cos(q7).*(a8+2.*d11.*cos(q8)).*sin(q2+2.*(q3+q5)).*sin(q8))),(-1).*cos(q1).*(a5.*cos(q2+2.*q3+q5)+2.*cos(q2+2.*(q3+q5)).*(d7+a8.*cos(q8)+d11.*cos(2.*q8))+(-2).*cos(q7).*(a8+2.*d11.*cos(q8)).*sin(q2+2.*(q3+q5)).*sin(q8)),(a8+2.*d11.*cos(q8)).*(cos(q7).*sin(q1)+cos(q1).*cos(q2+2.*(q3+q5)).*sin(q7)).*sin(q8),(-1).*cos(q1).*cos(q2+2.*(q3+q5)).*cos(q7).*(a8.*cos(q8)+2.*d11.*cos(2.*q8))+(a8.*cos(q8)+2.*d11.*cos(2.*q8)).*sin(q1).*sin(q7)+cos(q1).*(a8+4.*d11.*cos(q8)).*sin(q2+2.*(q3+q5)).*sin(q8),0;(-1).*cos(q1).*(a2.*sin(q2)+a3.*sin(q2+q3)+a4.*sin(q2+2.*q3)+a5.*sin(q2+2.*q3+q5)+(d7+a8.*cos(q8)+d11.*cos(2.*q8)).*sin(q2+2.*(q3+q5)))+(-1).*cos(q1).*cos(q2+2.*(q3+q5)).*cos(q7).*(a8+2.*d11.*cos(q8)).*sin(q8)+(a8+2.*d11.*cos(q8)).*sin(q1).*sin(q7).*sin(q8),(-1).*sin(q1).*(a2.*cos(q2)+a3.*cos(q2+q3)+a4.*cos(q2+2.*q3)+a5.*cos(q2+2.*q3+q5)+cos(q2+2.*(q3+q5)).*(d7+a8.*cos(q8)+d11.*cos(2.*q8))+(-1).*cos(q7).*(a8+2.*d11.*cos(q8)).*sin(q2+2.*(q3+q5)).*sin(q8)),(-1).*sin(q1).*(a3.*cos(q2+q3)+2.*(a4.*cos(q2+2.*q3)+a5.*cos(q2+2.*q3+q5)+cos(q2+2.*(q3+q5)).*(d7+a8.*cos(q8)+d11.*cos(2.*q8))+(-1).*cos(q7).*(a8+2.*d11.*cos(q8)).*sin(q2+2.*(q3+q5)).*sin(q8))),(-1).*sin(q1).*(a5.*cos(q2+2.*q3+q5)+2.*cos(q2+2.*(q3+q5)).*(d7+a8.*cos(q8)+d11.*cos(2.*q8))+(-2).*cos(q7).*(a8+2.*d11.*cos(q8)).*sin(q2+2.*(q3+q5)).*sin(q8)),(a8+2.*d11.*cos(q8)).*((-1).*cos(q1).*cos(q7)+cos(q2+2.*(q3+q5)).*sin(q1).*sin(q7)).*sin(q8),(-1).*cos(q2+2.*(q3+q5)).*cos(q7).*(a8.*cos(q8)+2.*d11.*cos(2.*q8)).*sin(q1)+(-1).*cos(q1).*(a8.*cos(q8)+2.*d11.*cos(2.*q8)).*sin(q7)+(a8+4.*d11.*cos(q8)).*sin(q1).*sin(q2+2.*(q3+q5)).*sin(q8),0;0,(-1).*a2.*sin(q2)+(-1).*a3.*sin(q2+q3)+(-1).*a4.*sin(q2+2.*q3)+(-1).*a5.*sin(q2+2.*q3+q5)+(-1).*(d7+a8.*cos(q8)+d11.*cos(2.*q8)).*sin(q2+2.*(q3+q5))+(-1).*cos(q2+2.*(q3+q5)).*cos(q7).*(a8+2.*d11.*cos(q8)).*sin(q8),(-1).*a3.*sin(q2+q3)+(-2).*(a4.*sin(q2+2.*q3)+a5.*sin(q2+2.*q3+q5)+(d7+a8.*cos(q8)+d11.*cos(2.*q8)).*sin(q2+2.*(q3+q5))+cos(q2+2.*(q3+q5)).*cos(q7).*(a8+2.*d11.*cos(q8)).*sin(q8)),(-1).*a5.*sin(q2+2.*q3+q5)+(-2).*(d7+a8.*cos(q8)+d11.*cos(2.*q8)).*sin(q2+2.*(q3+q5))+(-2).*cos(q2+2.*(q3+q5)).*cos(q7).*(a8+2.*d11.*cos(q8)).*sin(q8),(a8+2.*d11.*cos(q8)).*sin(q2+2.*(q3+q5)).*sin(q7).*sin(q8),(-1).*cos(q7).*(a8.*cos(q8)+2.*d11.*cos(2.*q8)).*sin(q2+2.*(q3+q5))+(-1).*cos(q2+2.*(q3+q5)).*(a8+4.*d11.*cos(q8)).*sin(q8),0;0,sin(q1),2.*sin(q1),2.*sin(q1),(-2).*cos(q1).*sin(q2+2.*(q3+q5)).*sin(q8).^2+(cos(q1).*cos(q2+2.*(q3+q5)).*cos(q7)+(-1).*sin(q1).*sin(q7)).*sin(2.*q8),2.*(cos(q7).*sin(q1)+cos(q1).*cos(q2+2.*(q3+q5)).*sin(q7)),sin(q1).*sin(q7).*sin(2.*q8)+(-1).*cos(q1).*(cos(2.*q8).*sin(q2+2.*(q3+q5))+cos(q2+2.*(q3+q5)).*cos(q7).*sin(2.*q8));0,(-1).*cos(q1),(-2).*cos(q1),(-2).*cos(q1),(-2).*sin(q1).*sin(q2+2.*(q3+q5)).*sin(q8).^2+(cos(q2+2.*(q3+q5)).*cos(q7).*sin(q1)+cos(q1).*sin(q7)).*sin(2.*q8),(-2).*cos(q1).*cos(q7)+2.*cos(q2+2.*(q3+q5)).*sin(q1).*sin(q7),(-1).*cos(2.*q8).*sin(q1).*sin(q2+2.*(q3+q5))+(-1).*(cos(q2+2.*(q3+q5)).*cos(q7).*sin(q1)+cos(q1).*sin(q7)).*sin(2.*q8);1,0,0,0,2.*cos(q2+2.*(q3+q5)).*sin(q8).^2+cos(q7).*sin(q2+2.*(q3+q5)).*sin(2.*q8),2.*sin(q2+2.*(q3+q5)).*sin(q7),cos(q2+2.*(q3+q5)).*cos(2.*q8)+(-1).*cos(q7).*sin(q2+2.*(q3+q5)).*sin(2.*q8)];
    Jmp=pinv(Js);
    dqs=Jmp*epo;
    qc=qc+alpha*dqs;
    for i=1:7
        if qc(i)<qlim(i,1)
            qc(i)=qlim(i,1);
            overlimit_flag(i,j)=1;
        elseif qc(i)>qlim(i,2)
            qc(i)=qlim(i,2);
            overlimit_flag(i,j)=1;
        end
    end
    
    if norm(qc-q)<step_tol
        break
    end
    q=qc;
    
    j=j+1;
end