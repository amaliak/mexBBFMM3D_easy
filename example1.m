function QH = example1(ExecName,Kernel,corlength,H)

% This is an example file with the input for using mexBBFMM3D to multiply
% the covariance matrix specified in kernelfun.hpp
% for a regular grid, with matrix H. 
% -----------
% Input:
%     ExecName : the name of the mexfile for the Kernel chosen
%     Kernel   : covariance type, e.g. 'GAUSSIAN'
%     corlength: correlation length, isotropic
%                anisotropy in z direction supported, see code
%     H        : matrix by which Kernel is multiplied
%
% Output:
%       ExecName.mexmaci64: executable for mex file for given configuration
%       QH       : Product of Kernel chosen by matrix H specified in input
%
% Type "help compilemex" for more information on the covariance options
% Type "help runmexBFMM3D" for more information on the input variables in this file 
clearvars -except H ExecName Kernel corlength

%%%%%%%%%%%%%%%%%%% PART MODIFIED BY THE USER  %%%%%%%%%%%%%%%%%%%%%%%%%%
TestingMode = true;
recompile = 1;

if recompile
    compilemex(ExecName,Kernel,corlength)
end


% Modify xloc and yloc for your own example
% 3-D locations, stored column-wise i.e. (x | y | z)
x = -62:4:62;
y = -62:4:62;
z = -9:3:9;
%7168

% FMM parameters
% increase to decrease relative error and increase accuracy and comp.cost
nCheb = 4;          % Number of Chebyshev nodes per dimension
level = 5;          % Level of FMM tree
use_chebyshev = 1;  % 1: chebyshev interpolation; 0: uniform interpolation

%%%%%%%%%%%%%%%%%%% END OF PART MODIFIED BY THE USER  %%%%%%%%%%%%%%%%%%%

%%%%%%%%% THE REST OF THIS M-FILE DOES NOT NEED TO BE MODIFIED %%%%%%%%%%

% Info on dimensions

nx = length(x); ny = length(y); nz=length(z);
Ns  = nx*ny*nz;    % Number of sources in simulation cell
Nf  = Ns;          % Number of fields in simulation cell
% Length of simulation cell (assumed to be a cube)
L = max(max(max(x) - min(x), max(y)-min(y)),max(z)-min(z));
source = zeros(nx*ny*nz,3);
for i = 1:ny
    for j = 1:nx
        for k = 1:nz
            source((i-1)*nx*nz + (j-1)*nz + k ,1) = x(i);
            source((i-1)*nx*nz + (j-1)*nz + k ,2) = y(j);
            source((i-1)*nx*nz + (j-1)*nz + k ,3) = z(k);
        end
    end
end

% check if H and Q are consistent
if abs((Ns-size(H,1)))>0
    disp('Dimension of H does not agree with Kernel dimenstion')
    keyboard
end

% run mexBBFMM2D

QH = runmexBBFMM3D(source,H,nCheb,L,level,ExecName,use_chebyshev,TestingMode);