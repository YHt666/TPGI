% Generality test for success rate, runtime and iterations.
close all; clear; clc

% Options
if_display_progress=1;
if_display_result=1;
if_plot=1;

% Choose ID of manipulator to be tested
% ID         Manipulator
% 1          planar 3-DOF
% 2          planar 4-DOF
% 3          spatial 6-DOF
% 4          spatial 7-DOF
ID_manipulator=4;

switch ID_manipulator
    case 1
        para=fun_para_3R;
    case 2
        para=fun_para_4R;
    case 3
        para=fun_para_6R;
    case 4
        para=fun_para_7R;
end

DH=para.DH;
qlim=para.qlim;
Uk=para.Uk;
robot=para.robot;

max_time=50e-3;
N=1e4;  % Number of cases
n=size(qlim,1);
q_rand=qlim(:,1)+(qlim(:,2)-qlim(:,1)).*rand(n,N);
q0_rand=qlim(:,1)+(qlim(:,2)-qlim(:,1)).*rand(n,N);
Td_set=robot.fkine((Uk*q_rand)').T;

q=zeros(n,N);
t=zeros(1,N);
k=zeros(1,N);
exitflag=zeros(1,N);
for i=1:N
    if if_display_progress
        if mod(i,100)==0
            disp(['Progress: ',num2str(i/N*100),'%'])
        end
    end
    Td=Td_set(:,:,i);
    q0=q0_rand(:,i);
    switch ID_manipulator
        case 1
            [q(:,i),exitflag(i),t(i),k(i)]=TPGI_P3D(Td,q0,para);
        case 2
            [q(:,i),exitflag(i),t(i),k(i)]=TPGI_P4D(Td,q0,para);
        case 3
            [q(:,i),exitflag(i),t(i),k(i)]=TPGI_S6D(Td,q0,para);
        case 4
            [q(:,i),exitflag(i),t(i),k(i)]=TPGI_S7D(Td,q0,para);
    end
end

solve_id=t<max_time & exitflag;
solve_rate=sum(solve_id)/N*100;
solve_time=t(solve_id)*1e3;
mean_solve_time=mean(solve_time,2);
solve_steps=k(solve_id);
mean_solve_steps=mean(solve_steps,2);
Q_solve_time=prctile(solve_time,[25 50 75]);
median_solve_time=Q_solve_time(2);
IQR_solve_time=Q_solve_time(3)-Q_solve_time(1);
Q_solve_steps=prctile(solve_steps,[25 50 75]);
median_solve_steps=Q_solve_steps(2);
IQR_solve_steps=Q_solve_steps(3)-Q_solve_steps(1);

if if_display_result
    disp('Success rate (%)')
    disp(solve_rate)
    disp('Runtime (ms), M (IQR)')
    disp([median_solve_time,IQR_solve_time])
    disp('Number of iterations, M (IQR)')
    disp([median_solve_steps,IQR_solve_steps])
end

if if_plot
    figure
    title('Runtime')
    histogram(solve_time,'Normalization','probability')
    figure
    title('Runtime')
    boxchart(solve_time)
    set(gca,'YScale','log')
    figure
    title('Iterations')
    histogram(solve_steps,'Normalization','probability')
    figure
    title('Iterations')
    boxchart(solve_steps)
    set(gca,'YScale','log')
end
