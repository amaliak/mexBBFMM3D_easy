function gridmesh = CreateRegMesh(grid)
% Function that creates vectors of (x,y,z) coordinate triples if given the unique x,y,z values
% Used in example1.m to generate the grid

nx = length(grid.x); ny = length(grid.y); nz=length(grid.z);

[X,Y,Z]=meshgrid(grid.x,grid.y,grid.z);
gridmesh.x = reshape(X,nx*ny*nz,1);
gridmesh.y = reshape(Y,nx*ny*nz,1);
gridmesh.z = reshape(Z,nx*ny*nz,1);
