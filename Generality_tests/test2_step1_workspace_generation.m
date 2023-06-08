close all; clear; clc

% Choose paremeters of manipulator to be tested
para=fun_para_7R;
robot=para.robot;
Uk=para.Uk;
Ut=para.Ut;
qlim=para.qlim;
row=size(Ut,1);
col=size(Ut,2);

%% Generate the seed workspace by Monte Carlo method
N=1e4;
q_rand_set=zeros(row,N);
for i=1:row
    q_rand_set(i,:)=qlim(i,1)+(qlim(i,2)-qlim(i,1)).*rand(1,N);
end
X_rand_set=zeros(3,N);
for i=1:N
    Te=robot.fkine(Uk*q_rand_set(:,i)).T;
    X_rand_set(:,i)=Te(1:3,4);
end

%% Gaussian growth of the seed worksapce
% Estimate bounding box size
box=[min(X_rand_set,[],2),max(X_rand_set,[],2)];
box_center=mean(box,2);
box_length=box(:,2)-box(:,1);
box_length=box_length*1.05;     % Increase by 5% to compensate for simulation error
box=[box_center-box_length/2,box_center+box_length/2];
% disp(box)
% a=para.a;
% d=para.d;
% sum(a+d)

% Discretization into cubic cells
lc=0.05;    % The cell side length is empirically determined to be about 1/20 of the arm span length
nc=ceil(box_length/lc);
nc=2*ceil(nc/2);    % Adjust the number of cells in each dimension to an even number

Nc=10;
% Build cellular data structure
cell0.id=[0;0;0];   % ID coordinate of the cell
cell0.np=0;     % Number of poses in the cell
cell0.X=[];     % End poses inside the cell
cell0.q=[];     % Corresponding joint angles
cell_dataset=repmat(cell0,nc(1),nc(2),nc(3));
for i=1:N
    id=word2id(X_rand_set(:,i),box_center,lc,nc);
    % Each cell stores a maximum of Nc poses
    if cell_dataset(id(1),id(2),id(3)).np<Nc
        cell_dataset(id(1),id(2),id(3)).np=cell_dataset(id(1),id(2),id(3)).np+1;
        cell_dataset(id(1),id(2),id(3)).X=[cell_dataset(id(1),id(2),id(3)).X,X_rand_set(:,i)];
        cell_dataset(id(1),id(2),id(3)).q=[cell_dataset(id(1),id(2),id(3)).q,q_rand_set(:,i)];
    end
end

nf_max=5;
w=1.01;
PC_list=[];
for i=1:nc(1)
    for j=1:nc(2)
        for k=1:nc(3)
            cell_dataset(i,j,k).id=[i;j;k];
            cell_dataset(i,j,k).p=id2word([i;j;k],box_center,lc,nc);
            if cell_dataset(i,j,k).np>0 && cell_dataset(i,j,k).np<Nc
                PC_list=[PC_list,[i;j;k]];
            end
        end
    end
end

%% Densification and growth
while ~isempty(PC_list)
    id=PC_list(:,1);
    sigma=(qlim(:,2)-qlim(:,1))/6;
    nf=0;
    if mod(size(PC_list,2),100)==0
        fprintf('Number of remaining PCs：%d\n',size(PC_list,2))
        fprintf('ID of the cell being filled：[%d,%d,%d]\n',id(1),id(2),id(3))
    end
    while cell_dataset(id(1),id(2),id(3)).np<Nc
        if nf>nf_max
            nf=0;
            sigma=sigma/w;
        end
        nf=nf+1;
        q0=cell_dataset(id(1),id(2),id(3)).q(:,randperm(size(cell_dataset(id(1),id(2),id(3)).q,2),1));
        qs=normrnd(q0,sigma);
        if all(qs>=qlim(:,1) & qs<=qlim(:,2))
            Ts=robot.fkine(Uk*qs).T;
            Xs=Ts(1:3,4);
            if all(Xs>=box(:,1) & Xs<=box(:,2))
                ids=word2id(Xs,box_center,lc,nc);
                cells=cell_dataset(ids(1),ids(2),ids(3));
                if cells.np<Nc
                    if cells.np==0
                        PC_list=[PC_list,ids];
                    end
                    cell_dataset(ids(1),ids(2),ids(3)).np=cells.np+1;
                    cell_dataset(ids(1),ids(2),ids(3)).X=[cells.X,Xs];
                    cell_dataset(ids(1),ids(2),ids(3)).q=[cells.q,qs];
                    if cell_dataset(ids(1),ids(2),ids(3)).np==Nc
                        PC_list(:,all(PC_list==ids))=[];
                    end
                end
                if all(ids==id)
                    nf=0;
                end
            end
        end
    end
end

cell_set=[];
ws=[];
for i=1:nc(1)
    for j=1:nc(2)
        for k=1:nc(3)
            if cell_dataset(i,j,k).np>0
                cell_set=[cell_set,cell_dataset(i,j,k)];
                ws=[ws,cell_dataset(i,j,k).X];
            end
        end
    end
end

%% save data
box_struct.center=box_center;
box_struct.lc=lc;
box_struct.nc=nc;
save box_struct.mat box_struct
save cell_dataset.mat cell_dataset
save cell_set.mat cell_set

%% Visualization
figure
robot.plot(zeros(1,col))
hold on
scatter3(ws(1,:),ws(2,:),ws(3,:),1,'filled')

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
