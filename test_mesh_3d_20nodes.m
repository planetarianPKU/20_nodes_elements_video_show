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
clc ; clear all ;
% Variables which can be changed
% Dimensions of the plate
L = 1 ;             % Length of the Plate along X-axes
B = 1 ;             % Breadth of the Plate along Y-axes
T = 1;
% Number of Elements required
Nx = 2 ;            % Number of Elements along X-axes
Ny = 2 ;            % Number of Elements along Y-axes
Nz = 2;
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
 plot3(test(:,1),test(:,2),test(:,3),'o','Markersize',20)
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
 
writematrix(nodes_new,'nodes_new.txt');
writematrix(coordinates,'coordinates.txt');
