function QH = example1(ExecName,grid,Kernel,corlength,H,TestingMode)
% This is an example file for ufing mexBBFMM3D to compute the product of a
% covariance kernel Q with matrix H, where Q is defined on an regular isotropic grid
% The grid must be provided as unique x,y and z locations and all (x,y,z)
% triplets will be created automatically
% -----------
%
% Example usage: 
%                grid.x = -62:4:62; grid.y = -62:4:62; grid.z = -9:3:9;
%                QH = example1('TESTNAME',grid,'GAUSSIAN',50,ones(7168,1),1)
% -----------
% Input:
%     ExecName : the name of the mexfile for the Kernel chosen
%     grid     : structure with vectors grid.x, grid.y, grid.z
%                each vector containing x,y and z coordinates
%                respectively
%     Kernel   : covariance type, e.g. 'GAUSSIAN'
%     corlength: correlation length, isotropic
%                anisotropy in z direction supported, see code
%     H        : matrix by which Kernel is multiplied
%     TestingMode: if set to 1, BBFMM is recompiled and runs in TestingMode in order to
%     determine parameters (nCheb) for desired accuracy. if set to 0, the 
%
% Output:
%       ExecName.mexmaci64: executable for mex file for given configuration
%       QH       : Product of Kernel chosen by matrix H specified in input
%
% Type "help compilemex" for more information on the covariance options
% Type "help runmexBFMM3D" for more information on the input variables in this file 

recompile = 0;
if isempty(TestingMode)
    TestingMode = true;
end

if TestingMode
    recompile = 1;
end

if recompile
    compilemex(ExecName,Kernel,corlength)
    save lastcompiled.mat ExecName Kernel corlength
else
    disp('ExecName, Kernel, Corlength ignored')
    disp('Using executable last compiled')
    disp('------------------------------')
    load lastcompiled.mat
    disp(['Executable Name is ',ExecName])
    disp(['Kernel type is ',Kernel])
    disp(['Correlation length in x,y,z is ',num2str(corlength)])
end

% FMM parameters
% Run code in TestingMode to get relative error of FMM vs direct
% multiplication; increase nCheb and Level to increase accuracy *and* comp.cost
nCheb = 4;          % Number of Chebyshev nodes per dimension
level = 5;          % Level of FMM tree
use_chebyshev = 1;  % 1: chebyshev interpolation; 0: uniform interpolation

%%%%%%%%% THE REST OF THIS M-FILE DOES NOT NEED TO BE MODIFIED %%%%%%%%%%

% Setting up the regular grid

nx = length(grid.x); ny = length(grid.y); nz=length(grid.z);
Ns  = nx*ny*nz;    % Number of sources in simulation cell

% Length of simulation cell (assumed to be a cube)
L = max(max(max(grid.x) - min(grid.x), max(grid.y)-min(grid.y)),max(grid.z)-min(grid.z));

gridmesh = CreateRegMesh(grid);
source = [gridmesh.x, gridmesh.y, gridmesh.z];

% check if H and Q are consistent
if abs((Ns-size(H,1)))>0
    disp('Dimension of H does not agree with Kernel dimenstion')
    keyboard
end



% run mexBBFMM2D

QH = runmexBBFMM3D(source,H,nCheb,L,level,ExecName,use_chebyshev,TestingMode);