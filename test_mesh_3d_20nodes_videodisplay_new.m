%coordinates是27节点，nodes是基于27节点的编号。接下来我们只要把节点编号和形函数对应即可。

% Purpose:
%         To Mesh a square/Rectangular plate to use in FEM Analysis
% Variable Description:
%           L - Length of the Plate along X-axes
%           B - Breadth of the Plate along Y-axes
%           Nx - Number of Elements along X-axes
%           Ny - Number of Elements along Y-axes
%           coordinates - The nodal coordinates of the mesh
%           -----> coordinates = [node X Y]
%           nodes - The nodal connectivity of the elements
%           -----> nodes = [node1 node2......]
%--------------------------------------------------------------------------
clc ; clear;
% Variables which can be changed
% Dimensions of the plate
L = 1 ;             % Length of the Plate along X-axes
B = 1 ;             % Breadth of the Plate along Y-axes
T = 1;
% Number of Elements required
Nx = 3 ;            % Number of Elements along X-axes
Ny = 3 ;            % Number of Elements along Y-axes
Nz = 3;
%---------------------------------------
% From here dont change
nel = Nx*Ny*Nz ;        % Total Number of Elements in the Mesh
nnel = 20 ;           % Number of nodes per Element
% Number of points on the Length and Breadth
npx = 2*Nx+1 ;
npy = 2*Ny+1 ;
npz = 2*Nz+1;
nnode = npx*npy*npz-Nx*Ny*Nz-(Nx+1)*Ny*Nz-(Ny+1)*Nx*Nz-(Nz+1)*Nx*Ny ;      % Total Number of Nodes in the Mesh
nnode_fake = npx*npy*npz;
% Discretizing the Length and Breadth of the plate

% To get the Nodal Connectivity Matrix
% Valid point index in each 20-nodes element, comparing with a 27-nodes element
%3
%|
%2
%|
%1---2----3
oklst=[1 1 1;2 1 1;3 1 1;1 2 1;3 2 1;1 3 1;2 3 1;3 3 1; 1 1 3;2 1 3;3 1 3;1 2 3;3 2 3;1 3 3;2 3 3;3 3 3;1 1 2;3 1 2;3 3 2;1 3 2];
%show the valid points
figure,plot3(oklst(:,1),oklst(:,2),oklst(:,3),'o')

i=0;j=0;k=0;
%index是循环中用于统计20节点而不是27节点的计数器
index=0;
index_fake=0;
fakelst=[];
for z=0:1/(npz-1):T
    k=k+1;
    if k > 3
        k=2;
    end
    j=0;
    for y=0:1/(npy-1):B
        j=j+1;
        if j > 3
            j=2;
        end
        i=0;
        for x=0:1/(npx-1):L
            i=i+1;
            index_fake=index_fake+1;
            
            if i > 3
                i=2;
            end
            for okindex=1:length(oklst)
                if [i,j,k]==oklst(okindex,:)
                    index=index+1;
                    coordinates(index,1)=x;coordinates(index,2)=y;coordinates(index,3)=z;
                    flag = 1;
                end
            end
            if flag == 0
                fakelst=[fakelst index_fake];
            end
            flag = 0;
        end
    end
end
NodeNo = 1:nnode_fake ;
nodes = zeros(nel,nnel);

NodeNo = reshape(NodeNo,npx,npy,npz);
nodes(:,1) = reshape(NodeNo(3:2:npx,1:2:npy-2,3:2:npz),nel,1);
nodes(:,2) = reshape(NodeNo(3:2:npx,3:2:npy,3:2:npz),nel,1);
nodes(:,3) = reshape(NodeNo(1:2:npx-2,3:2:npy,3:2:npz),nel,1);
nodes(:,4) = reshape(NodeNo(1:2:npx-2,1:2:npy-2,3:2:npz),nel,1);
nodes(:,5) = reshape(NodeNo(3:2:npx,1:2:npy-2,1:2:npz-2),nel,1);
nodes(:,6) = reshape(NodeNo(3:2:npx,3:2:npy,1:2:npz-2),nel,1);
nodes(:,7) = reshape(NodeNo(1:2:npx-2,3:2:npy,1:2:npz-2),nel,1);
nodes(:,8) = reshape(NodeNo(1:2:npx-2,1:2:npy-2,1:2:npz-2),nel,1);

nodes(:,9) = reshape(NodeNo(3:2:npx,2:2:npy-1,3:2:npz),nel,1);
nodes(:,10) = reshape(NodeNo(2:2:npx-1,3:2:npy,3:2:npz),nel,1);
nodes(:,11) = reshape(NodeNo(1:2:npx-2,2:2:npy-1,3:2:npz),nel,1);
nodes(:,12) = reshape(NodeNo(2:2:npx-1,1:2:npy-2,3:2:npz),nel,1);

nodes(:,13) = reshape(NodeNo(3:2:npx,2:2:npy-1,1:2:npz-2),nel,1);
nodes(:,14) = reshape(NodeNo(2:2:npx-1,3:2:npy,1:2:npz-2),nel,1);
nodes(:,15) = reshape(NodeNo(1:2:npx-2,2:2:npy-1,1:2:npz-2),nel,1);
nodes(:,16) = reshape(NodeNo(2:2:npx-1,1:2:npy-2,1:2:npz-2),nel,1);

nodes(:,17) = reshape(NodeNo(3:2:npx,1:2:npy-2,2:2:npz-1),nel,1);
nodes(:,18) = reshape(NodeNo(3:2:npx,3:2:npy,2:2:npz-1),nel,1);
nodes(:,19) = reshape(NodeNo(1:2:npx-2,3:2:npy,2:2:npz-1),nel,1);
nodes(:,20) = reshape(NodeNo(1:2:npx-2,1:2:npy-2,2:2:npz-1),nel,1);

% test=coordinates(reshape(nodes,[],1),:)
% plot3(test(:,1),test(:,2),test(:,3),'*')
% for i = 1:length(reshape(nodes,[],1))
%     text(test(i,1),test(i,2),test(i,3),num2str(nodes(i)));
% end
nodes_new=zeros(nel,nnel);
for i = 1:nel
    for j = 1:nnel
        nodes_new(i,j) = nodes(i,j)-length(fakelst(fakelst<nodes(i,j)));
        
    end
end
test=coordinates(reshape(nodes_new,[],1),:);
plot3(test(:,1),test(:,2),test(:,3),'o','Markersize',20);
hold on;plot3([0 1],[0 0],[0 0],'LineWidth',3,'color','black');
hold on;plot3([0 0],[0 1],[0 0],'LineWidth',3,'color','black');
hold on;plot3([1 0],[1 1],[0 0],'LineWidth',3,'color','black');
hold on;plot3([1 1],[0 1],[0 0],'LineWidth',3,'color','black');

hold on;plot3([0 1],[0 0],[1 1],'LineWidth',3,'color','black');
hold on;plot3([0 0],[0 1],[1 1],'LineWidth',3,'color','black');
hold on;plot3([1 0],[1 1],[1 1],'LineWidth',3,'color','black');
hold on;plot3([1 1],[0 1],[1 1],'LineWidth',3,'color','black');

hold on;plot3([0 0],[0 0],[0 1],'LineWidth',3,'color','black');
hold on;plot3([1 1],[0 0],[0 1],'LineWidth',3,'color','black');
hold on;plot3([0 0],[1 1],[0 1],'LineWidth',3,'color','black');
hold on;plot3([1 1],[1 1],[0 1],'LineWidth',3,'color','black');
grid on; grid minor;
xlabel('x');ylabel('y');
for i = 1:length(reshape(nodes_new,[],1))
    text(test(i,1),test(i,2),test(i,3),num2str(nodes_new(i)),'Fontsize',20);
end

touch_x = [];
touch_y = [];
touch_z = [];


touch_all_x = [];
touch_all_y = [];
touch_all_z = [];

for ele_index = 1:Nx * Ny * Nz
    %取得该单元的所有节点
    ele_nodes = nodes_new(ele_index,:);
    %取得该单元右上角节点的坐标
    check_point =  coordinates(ele_nodes(2),:);
    %取得该单元的xyz三个相邻的单元的索引，如果该单元在边界上，那么可能没有相邻单元，需要后期判断是否在边缘来纠正
    
    x_check = check_point(1);
    y_check = check_point(2);
    z_check = check_point(3);
    
    touch_x = ele_index + 1;
    touch_y = ele_index + Ny;
    touch_z = ele_index + Nx * Ny;
    %按单元顺序储存xyz方向的相邻单元
    
    if x_check  < L
        touch_all_x = [touch_all_x touch_x];
    else
        touch_all_x = [touch_all_x -1];
    end
    if y_check  < B
        touch_all_y = [touch_all_y touch_y];
    else
        touch_all_y = [touch_all_y -1];
    end
    if z_check  < T
        touch_all_z = [touch_all_z touch_z];
    else
        touch_all_z = [touch_all_z -1];
    end
end
%检查是否等长，如果不等长会报错
touch_all_xyz = [touch_all_x;touch_all_y;touch_all_z];
%从一号开始，递归创建一个要画的列表，先不急于处理动画
begin = 1;
%创建一个当前迭代的要生长出的单元索引的列表，以及当前迭代刚生长出来的单元索引的列表
print_ind_lst_now = [begin];
print_ind_lst_next = [];

nstep_iter = round(2*sqrt(3)*max([Nx,Ny,Nz]));
len_once_print = 10;
store_matrix = zeros(len_once_print,nstep_iter);


%先处理初始化
begin_ind =1;



store_matrix(1,1) = begin_ind;
end_index_tag = 1;
for it = 1:nstep_iter
    for n = 1:length(print_ind_lst_now)
        
        if print_ind_lst_now(n) ~= -1
            ind_x = touch_all_xyz(1,print_ind_lst_now(n));
            ind_y = touch_all_xyz(2,print_ind_lst_now(n));
            ind_z = touch_all_xyz(3,print_ind_lst_now(n));
            add_item = [ind_x ind_y ind_z];
            
            %删除不用画的假单元
            add_item(add_item == -1) = [];
            print_ind_lst_next = [print_ind_lst_next add_item];
            print_ind_lst_next = unique(print_ind_lst_next);
        end
        
    end
    if isempty(print_ind_lst_next)
        break
    else
        store_matrix(1:length(print_ind_lst_next),it+1) = print_ind_lst_next;
        print_ind_lst_now = print_ind_lst_next;
        print_ind_lst_next = [];
        end_index_tag = end_index_tag +1;
    end
end

store_matrix(:,end_index_tag+1:end) = [];



i=0;j=0;k=0;
index=0;
index_fake=0;

fakelst=[];

indexlst = nodes_new
positionlst = coordinates;

out = VideoWriter('vdo_show.avi');
out.FrameRate = 30;
open(out);

figure('color','none');
set(gca,'Visible','off')
set(gcf, 'units','normalized','position',[0 0 1 1]);

box off
set(gca,'xtick',[])
set(gca,'ytick',[])
set(gca,'ztick',[])

xlim([-.1 1.1]*L);
ylim([-.1 1.1]*B);
zlim([-.1 1.1]*T);
axis equal;
view(45,30);



hold on;plot3([0 1],[0 0],[0 0],'LineWidth',3,'color','white');
hold on;plot3([0 0],[0 1],[0 0],'LineWidth',3,'color','white');
hold on;plot3([1 0],[1 1],[0 0],'LineWidth',3,'color','white');
hold on;plot3([1 1],[0 1],[0 0],'LineWidth',3,'color','white');

hold on;plot3([0 1],[0 0],[1 1],'LineWidth',3,'color','white');
hold on;plot3([0 0],[0 1],[1 1],'LineWidth',3,'color','white');
hold on;plot3([1 0],[1 1],[1 1],'LineWidth',3,'color','white');
hold on;plot3([1 1],[0 1],[1 1],'LineWidth',3,'color','white');

hold on;plot3([0 0],[0 0],[0 1],'LineWidth',3,'color','white');
hold on;plot3([1 1],[0 0],[0 1],'LineWidth',3,'color','white');
hold on;plot3([0 0],[1 1],[0 1],'LineWidth',3,'color','white');
hold on;plot3([1 1],[1 1],[0 1],'LineWidth',3,'color','white');
% grid on; grid minor;
xlabel('x');ylabel('y');





%这里稍微改一下，把这个三重循环写成针对一个单元的，然后再套一个大循环
view_idx = 0;
zoom_index = 0;
zoom_ind = 1;
for step = 1:size(store_matrix,2)
    
    print_this_step = store_matrix(:,step);
    print_this_step(print_this_step == 0) = [];
    %开始画该步骤的每个元素，从这里开始才开始套原来的程序
                                
                                if step > 1 && step < size(store_matrix,2)
                                zoom_ind = zoom_ind +1;
                                end
                              zoom_rate = (size(store_matrix,2)-2);
                                                                xlim([-.1 1.1]*L);
                                ylim([-.1 1.1]*B);
                                zlim([-.1 1.1]*T);
%                                 xlim([-.1 1.1]*zoom_ind*L/zoom_rate);
%                                 ylim([-.1 1.1]*zoom_ind*B/zoom_rate);
%                                 zlim([-.1 1.1]*zoom_ind*T/zoom_rate);
    for ele_print = 1:length(print_this_step)
        i=0;j=0;k=0;
        index=0;
        index_fake=0;
        fakelst=[];
       nodes_in_ele = nodes_new(print_this_step(ele_print),:);
    %零点起始点的节点
    ind_node_begin = nodes_in_ele(8);
    position_begin = coordinates(ind_node_begin,:);
    x_begin = position_begin(1);
    y_begin = position_begin(2);
    z_begin = position_begin(3);     
        x_last = x_begin;
        y_last = y_begin;
        z_last = z_begin;
        all_index = 0;
        
        
        
        for z=0+z_begin:1/(npz-1):T/Nz+z_begin
            k=k+1;
            
            j=0;
            for y=0+y_begin:1/(npy-1):B/Ny+y_begin
                j=j+1;
                
                i=0;
                for x=0+x_begin:1/(npx-1):L/Nx+x_begin
                    i=i+1;
                    index_fake=index_fake+1;
                    %挨个搜索确定
                    
                    for okindex=1:length(oklst)
                        %这是用于27节点修正为20节点，所以需要保留
                        if [i,j,k] == oklst(okindex,:)
                            %need x,x_last,y,y_last,z,z_last
                            index = index + 1;
                            %动画控制
                            animation_speed = 30;

                            for kkkk = 1:animation_speed
                                view_idx = view_idx +1;
                                hhh = plot3([x_last x_last+kkkk*(x-x_last)/animation_speed],[y_last y_last+kkkk*(y-y_last)/animation_speed],[z_last z_last+kkkk*(z-z_last)/animation_speed],'--o','color','yellow','markersize',15,'linewidth',3);
                                zoom_rate = (size(store_matrix,2)-2);
                                
zoom_index = zoom_index +1;

                                %axis equal;
                                view(45+view_idx/3,30);
                                %pause(0.1);
                                hold on;
                                %这里是要做什么，为什么要等于30的时候检测？
                                if kkkk == animation_speed
                                    
                                    plot3(x,y,z,'o','markeredgecolor','white','markerfacecolor','red','markersize',15);
                                    %indexlst的第二列是选取的铆钉点
                                    if index == 20%ismember(nodes_in_ele(index),nodes_in_ele(:,2))
                                       % el_idx = find(indexlst(:,2)==nodes_in_ele(index));
                                        nodes_in_el_idx = nodes_in_ele;
                                        p1 = positionlst(nodes_in_el_idx(1),:);
                                        p2 = positionlst(nodes_in_el_idx(2),:);
                                        p3 = positionlst(nodes_in_el_idx(3),:);
                                        p4 = positionlst(nodes_in_el_idx(4),:);
                                        p5 = positionlst(nodes_in_el_idx(5),:);
                                        p6 = positionlst(nodes_in_el_idx(6),:);
                                        p7 = positionlst(nodes_in_el_idx(7),:);
                                        p8 = positionlst(nodes_in_el_idx(8),:);
                                        hold on;plot3([p1(1) p2(1)],[p1(2) p2(2)],[p1(3) p2(3)],'--','LineWidth',1,'color','white');
                                        hold on;plot3([p2(1) p3(1)],[p2(2) p3(2)],[p2(3) p3(3)],'--','LineWidth',1,'color','white');
                                        hold on;plot3([p3(1) p4(1)],[p3(2) p4(2)],[p3(3) p4(3)],'--','LineWidth',1,'color','white');
                                        hold on;plot3([p4(1) p1(1)],[p4(2) p1(2)],[p4(3) p1(3)],'--','LineWidth',1,'color','white');
                                        
                                        hold on;plot3([p5(1) p6(1)],[p5(2) p6(2)],[p5(3) p6(3)],'--','LineWidth',1,'color','white');
                                        hold on;plot3([p6(1) p7(1)],[p6(2) p7(2)],[p6(3) p7(3)],'--','LineWidth',1,'color','white');
                                        hold on;plot3([p7(1) p8(1)],[p7(2) p8(2)],[p7(3) p8(3)],'--','LineWidth',1,'color','white');
                                        hold on;plot3([p8(1) p5(1)],[p8(2) p5(2)],[p8(3) p5(3)],'--','LineWidth',1,'color','white');
                                        
                                        hold on;plot3([p1(1) p5(1)],[p1(2) p5(2)],[p1(3) p5(3)],'--','LineWidth',1,'color','white');
                                        hold on;plot3([p2(1) p6(1)],[p2(2) p6(2)],[p2(3) p6(3)],'--','LineWidth',1,'color','white');
                                        hold on;plot3([p3(1) p7(1)],[p3(2) p7(2)],[p3(3) p7(3)],'--','LineWidth',1,'color','white');
                                        hold on;plot3([p4(1) p8(1)],[p4(2) p8(2)],[p4(3) p8(3)],'--','LineWidth',1,'color','white');
                                    end
                                    
                                end
                                F = getframe(gcf);
                                writeVideo(out,F);
                                delete(hhh);
                            end
                            
                            
                            x_last = x;
                            y_last = y;
                            z_last = z;
                            
                            flag = 1;
                            
                            all_index = all_index+1;
                        end
                    end
                    if flag == 0
                        fakelst=[fakelst index_fake];
                    end
                    
                    flag = 0;
                    
                end
            end
        end
    end
end
close(out);
