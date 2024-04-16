% Read Meshes
 fileID = fopen('MESH.dat','r');

        %%% Just for fun meshes %%%
% fileID = fopen('./meshes/0.01/square.mesh','r');
%  fileID = fopen('./meshes/0.001/square.mesh','r');
% fileID = fopen('./meshes/0.0001/square.mesh','r');
% commasign = char(44); % ASCII for ','

% obtain numbers of all nodes and inner nodes
allV = fscanf(fileID,'%d',[1 1]);
innerV = fscanf(fileID,'%d',[1 1]);

% obtain (x,y) coordinates
formatCoordinate = '%f %f'; 
sizeCoor = [2 allV];
coordinates = fscanf(fileID,formatCoordinate,sizeCoor);
coordinates = coordinates';

% obtain a number of triangles
tri_num = fscanf(fileID,'%d',[1 1]);

% obtain node numbers
formatNodes = '%d %d %d'; 
sizeNodes = [3 tri_num];
ndc = fscanf(fileID,formatNodes,sizeNodes);
ndc = ndc';

fclose(fileID);

% Assign coordinates to x and y
x = coordinates(:,1);
y = coordinates(:,2);


