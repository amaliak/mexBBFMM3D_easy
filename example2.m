function QH = example2(ExecName,grid,Kernel,corlength,H,TestingMode)
% This is an example file for ufing mexBBFMM3D to compute the product of a
% covariance kernel Q with matrix H, where Q is defined on an IRREGULAR grid
% with the option for anisotropy in the z direction
% Grid must be provided as (x,y,z) triplets
% -----------
% Example usage: 
%                load('./coord_htr.mat')
%                grid.x = x_htr; grid.y = y_htr; grid.z = z_htr;
%                QH = example1('TESTNAME',grid,'GAUSSIAN',50,ones(7168,1))
% -----------
%     ExecName : the name of the mexfile for the Kernel chosen
%     grid     : structure with vectors grid.x, grid.y, grid.z
%                each vector containing all x,y and z coordinates
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

%%%%%%%%%%%%%%%%%%% PART MODIFIED BY THE USER  %%%%%%%%%%%%%%%%%%%%%%%%%%
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
% increase to decrease relative error and increase accuracy and comp.cost
nCheb = 4;          % Number of Chebyshev nodes per dimension
level = 5;          % Level of FMM tree
use_chebyshev = 1;  % 1: chebyshev interpolation; 0: uniform interpolation


%%%%%%%%% THE REST OF THIS M-FILE DOES NOT NEED TO BE MODIFIED %%%%%%%%%%

% Setting up the regular grid
% Check that grid coordinates are the same size
if (length(grid.x)-length(grid.y))>0 || (length(grid.z)-length(grid.y))>0 || (length(grid.x)-length(grid.z))>0
    disp('Problem in grid dimensions, length of x,y,z vectors must agree')
end

Ns  = length(grid.x);    
Nf  = Ns;          
% Length of simulation cell (assumed to be a cube)
L = max(max(max(grid.x) - min(grid.x), max(grid.y)-min(grid.y)),max(grid.z)-min(grid.z));
source = [grid.x, grid.y, grid.z];

% Anisotropy
% Scaling the vertical coordinates with the correlation lengths 
% is equivalent to having anisotropic correlation lengths
lx = corlength; lz=10;
grid.z = grid.z * (lx/lz);

% check if H and Q are consistent
if abs((Ns-size(H,1)))>0
    disp('Dimension of H does not agree with Kernel dimenstion')
    keyboard
end

% Compute matrix-vectors product QH
QH = runmexBBFMM3D(source,H,nCheb,L,level,ExecName,use_chebyshev,TestingMode);

 
 

