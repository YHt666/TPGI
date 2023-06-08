close all; clear; clc

% Options
if_display_progress=1;
if_save_data=0;
if_display_result=1;
if_plot=1;

para=fun_para;
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
method_sta(6).q=0;
method_name={'KDL-P','KDL','SQP-P','SQP','TPGI-P','TPGI'};

for m_id=6  % Select methods to be tested
% for m_id=1:6  % Test all methods
    disp(['Start program of ',method_name{m_id},'.'])
    q=zeros(n,N);
    t=zeros(1,N);
    k=zeros(1,N);
    exitflag=zeros(1,N);
    for i=1:N
        if if_display_progress
            if mod(i,100)==0
                disp(['Progress of ',method_name{m_id},': ',num2str(i/N*100),'%'])
            end
        end
        
        Td=Td_set(:,:,i);
        q0=q0_rand(:,i);
        
        % Methods
        switch m_id
            case 1
                [q(:,i),exitflag(i),t(i),k(i)]=KDL_P(Td,q0,para);
            case 2
                [q(:,i),exitflag(i),t(i),k(i)]=KDL(Td,q0,para);
            case 3
                [q(:,i),exitflag(i),t(i),k(i)]=SQP_P(Td,q0,para);
            case 4
                [q(:,i),exitflag(i),t(i),k(i)]=SQP(Td,q0,para);
            case 5
                [q(:,i),exitflag(i),t(i),k(i)]=TPGI_P(Td,q0,para);
            case 6
                [q(:,i),exitflag(i),t(i),k(i)]=TPGI(Td,q0,para);
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
    
    % Save data
    if if_save_data
        method_sta(m_id).name=method_name{m_id};
        method_sta(m_id).q=q;
        method_sta(m_id).t=t;
        method_sta(m_id).k=k;
        method_sta(m_id).solve_id=solve_id;
        method_sta(m_id).solve_rate=solve_rate;
        method_sta(m_id).solve_time=solve_time;
        method_sta(m_id).mean_solve_time=mean_solve_time;
        method_sta(m_id).Q_solve_time=Q_solve_time;
        method_sta(m_id).median_solve_time=median_solve_time;
        method_sta(m_id).IQR_solve_time=IQR_solve_time;
        method_sta(m_id).solve_steps=solve_steps;
        method_sta(m_id).mean_solve_steps=mean_solve_steps;
        method_sta(m_id).Q_solve_steps=Q_solve_steps;
        method_sta(m_id).median_solve_steps=median_solve_steps;
        method_sta(m_id).IQR_solve_steps=IQR_solve_steps;
        save method_sta.mat method_sta
        disp(['Statistics of ',method_name{m_id},' is saved.'])
    end
    
    if if_display_result
        disp(['Results of ',method_name{m_id},':'])
        disp('Success rate (%)')
        disp(solve_rate)
        disp('Runtime (ms), M (IQR)')
        disp([median_solve_time,IQR_solve_time])
        disp('Number of iterations, M (IQR)')
        disp([median_solve_steps,IQR_solve_steps])
    end
    
    if if_plot
        figure
        title(['Runtime of ',method_name{m_id}])
        histogram(solve_time,'Normalization','probability')
        figure
        title(['Runtime of ',method_name{m_id}])
        boxchart(solve_time)
        set(gca,'YScale','log')
        figure
        title(['Iterations of ',method_name{m_id}])
        histogram(solve_steps,'Normalization','probability')
        figure
        title(['Iterations of ',method_name{m_id}])
        boxchart(solve_steps)
        set(gca,'YScale','log')
    end
end
