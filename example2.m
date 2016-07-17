function QH = example2(ExecName,Kernel,corlength,H)
% This is an anisotropic, irregular grid example file 
% with the input for using mexBBFMM3D to multiply
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
% %% Compile mex code if new kernel
TestingMode = true;
recompile = 1;

if recompile
    compilemex(ExecName,Kernel,corlength)
end


% Modify xloc and yloc for your own example
% 3-D locations, stored column-wise i.e. (x | y | z)
load('./coord_htr.mat')
x_loc =  x_htr;
y_loc =  y_htr;
z_loc = z_htr;

% anisotropy
% Scaling the vertical coordinates with the correlation lengths 
% is equivalent to having anisotropic correlation lengths
lx = 60; lz=10;
z_loc = z_loc * (lx/lz);

% FMM parameters
% increase to decrease relative error and increase accuracy and comp.cost
nCheb = 4;          % Number of Chebyshev nodes per dimension
level = 5;          % Level of FMM tree
use_chebyshev = 1;  % 1: chebyshev interpolation; 0: uniform interpolation

%%%%%%%%%%%%%%%%%%% END OF PART MODIFIED BY THE USER  %%%%%%%%%%%%%%%%%%%

%%%%%%%%% THE REST OF THIS M-FILE DOES NOT NEED TO BE MODIFIED %%%%%%%%%%

L = max(max(max(x_loc) - min(x_loc), max(y_loc)-min(y_loc)),max(z_loc)-min(z_loc));% Length of simulation cell (assumed to be a cube)
source = [x_loc, y_loc,z_loc];

% Info on dimensions
Ns  = length(source);    % Number of sources in simulation cell
Nf  = Ns;    % Number of fields in simulation cell

% check if H and Q are consistent
if abs((Ns-size(H,1)))>0
    disp('Dimension of H does not agree with Kernel dimenstion')
    keyboard
end

% Compute matrix-vectors product QH
QH = runmexBBFMM3D(source,H,nCheb,L,level,ExecName,use_chebyshev,TestingMode);

 
 

