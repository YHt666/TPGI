close all; clear; clc

max_time=50e-3;

para=fun_para;
robot=para.robot;
Uk=para.Uk;
Ut=para.Ut;
qlim=para.qlim;

%% Runtime distribution
load cell_set.mat
load box_struct.mat
box_center=box_struct.center;
lc=box_struct.lc;
nc=box_struct.nc;
cs=cell_set;
Nc=size(cs(1).q,2);
t=zeros(1,Nc);
exitflag=zeros(1,Nc);
for i=1:length(cs)
    if mod(i,100)==0
        disp([num2str(i/length(cs)*100),'%'])
    end
    for j=1:Nc
        qd=cs(i).q(:,j);
        Td=robot.fkine(Uk*qd).T;
        q0=qlim(:,1)+(qlim(:,2)-qlim(:,1)).*rand(7,1);
        [~,exitflag(j),t(j),~]=TPGI(Td,q0,para);
    end
    % Runtime of poses that cannot be solved is set to the upper limit
    t(~exitflag)=max_time;
    % Runtime for each cell is the average time taken to solve its poses.
    Do=mean(t,2)*1e3;
    cs(i).Do=Do;
end

%% Save data
save cs.mat cs

%% Plot and print
close all
Do_list=sort([cs(:).Do],'descend');
p=zeros(3,length(cs));
for i=1:length(cs)
    p(:,i)=id2word(cs(i).id,box_center,lc,nc);
end
% Axonometric view
figure
scatter3(p(1,:),p(2,:),p(3,:),10,[cs(:).Do],'filled')
axis equal
caxis([0 max_time*1e3])
colorbar
xlabel('X (m)'); ylabel('Y (m)'); zlabel('Z (m)')
% X-Z sectional view
select=p(2,:)>0 & p(2,:)<lc;
cs_xz=cs(select);
p_xz=p(:,select);
Do_xz=[cs_xz(:).Do];
figure
scatter3(p_xz(1,:),p_xz(2,:),p_xz(3,:),20,Do_xz,'filled')
axis equal
caxis([0 max_time*1e3])
colorbar
xlabel('X (m)'); ylabel('Y (m)'); zlabel('Z (m)')
view(0,0)
% X-Y sectional view
select=p(3,:)>0 & p(3,:)<lc;
cs_xy=cs(select);
p_xy=p(:,select);
Do_xy=[cs_xy(:).Do];
figure
scatter3(p_xy(1,:),p_xy(2,:),p_xy(3,:),10,Do_xy,'filled')
axis equal
caxis([0 max_time*1e3])
colorbar
xlabel('X (m)'); ylabel('Y (m)'); zlabel('Z (m)')
view(2)

workspace_coverage=(1-sum(Do_list>30)/length(Do_list))*100;
disp('Workspace coverage (%)')
disp(workspace_coverage)

%% function
% Mapping of world coordinates to struct index coordinates
function id=word2id(p,box_center,lc,nc)
g=word2cell(p,box_center,lc);
id=cell2id(g,nc);
end

% Mapping of index coordinates to world coordinates
function p=id2word(id,box_center,lc,nc)
g=id2cell(id,nc);
p=cell2word(g,box_center,lc);
end

% Mapping of cellular coordinate system to world coordinate system
function p=cell2word(g,box_center,lc)
p=box_center+sign(g).*(abs(g)-0.5)*lc;
end

% Mapping of cell coordinates to structure index coordinates
function id=cell2id(g,nc)
id=g+nc/2+(1-sign(g))/2;
end

% Mapping of structure index coordinates to cell coordinates
function g=id2cell(id,nc)
g=id-(nc+1-sign(id-(nc+1)/2))/2;
end

% Mapping of world coordinate system to cellular coordinate system
function g=word2cell(p,box_center,lc)
delta=p-box_center;
g=sign(delta).*ceil(abs(delta)/lc);
end
