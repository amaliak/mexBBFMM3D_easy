function [U,S,V] = RandomizedCondSVDFMM(grid,Kernel,Corlength,Corlengthz,N)

tic
% Code for low-rank decomposition of covariance kernel Q
% via randomized algorithm
% No parallelization for small problems
% Performs fast SVD for covariance matrix and finds approximate U, S, and V
% To reduce computational cost, randSVD is combined with BBFMM3D
% Care should be taken to compile  BBFMM3D for the covariance
% matrix we are interested in - the Q matrix is calculated internally in
% the BBFMM3D code which performs fast matrix multiplication of the
% covariance with random vectors, as required by randSVD
%
% Documentation: randomized SVD algorithm
%
% Amalia Kokkinaki
% modification of code written by P. Kitanidis with input from Harry J. Lee's
% function rnd_eig_fft.m and help from Ruoxi Wang on the use of mexBBFMM3D
%
% Example usage:
%                grid.x = -12:6:12; grid.y = -12:6:12; grid.z = -6:3:6;
%                gridmesh = CreateRegMesh(grid);
%                [U,S,V] = RandomizedCondSVDFMM(gridmesh,'GAUSSIAN',100,10,3,9);
% Input:
%        grid     : structure with vectors grid.x, grid.y, grid.z
%                each vector containing all (x,y,z) triplets 
%               The meshed grid can be created from x,y and z vectors by using 
%          function gridmesh = CreateRegMesh(grid);
%        N : rank of reduced rank svd 
%        Kernel Q: covariance type, see compilemex for options
%        corlength: correlation length in x and y isotropic
%        corlengthz: correlation length in z
%        
% Output: 
%       U,V N first eigenbases
%       S: diagonal matrix with N first eigenvectors

%% setting up BBFMM3D

% compile mexBBFMM3D for covariance of interest
% This creates the executable that will evaluate products QH
% needed for randSVD of matrix Q
% For help with compilation type help compilemex
ExecName = 'RandSVDBBF3D';
compilemex(ExecName,Kernel,Corlength)

m  = length(grid.x);

% Anisotropy
% Scaling the vertical coordinates with the correlation lengths 
% is equivalent to having anisotropic correlation lengths
lx = Corlength; lz=Corlengthz;
grid.z_an = grid.z * (lx/lz);

% Length of simulation cell (assumed to be a cube)
L = max(max(max(grid.x) - min(grid.x), max(grid.y)-min(grid.y)),max(grid.z_an)-min(grid.z_an));
source = [grid.x, grid.y, grid.z_an];

% FMM parameters
Ns  = m;            % Number of sources in simulation cell
nCheb = 4;          % Number of Chebyshev nodes per dimension
level = ceil(log10(Ns));
use_chebyshev = 1;  % 1: chebyshev interpolation; 0: uniform interpolation


% Run BBFMM3D in testing mode to determine nCheb for desired accuracy

disp('Running test for accuracy')
disp('Also creating the Pre-computation files')
TestingMode = 1;
testOmega=randn(Ns,1);
runmexBBFMM3D(source,testOmega(:,1),nCheb,L,level,ExecName,use_chebyshev,TestingMode);

check = 1;%input('Is accuracy good enough?');
while ~check
	disp(['Current number of Cheb nodes is',num2str(nCheb)])
	nCheb = input('New number of Cheb nodes is');
	disp('Running test for accuracy with new nCheb')
    runmexBBFMM3D(source,testOmega(:,1),nCheb,L,level,ExecName,use_chebyshev,TestingMode);
	check = input('Is accuracy good enough?');
end

clear  QH QHexact

%% Randomized SVD 
disp('Starting randSVD')

if N>m
    warning('N set to minimum of dimensions')
    N=m;
end
 

% Y = (A*A')^q*A*Omega; % A(A(AO))) 3x(mlogm N)
% since A = Cov, A=A', Y=(Q^2q)*Q*Omega  
% for q = 1 Y = Q*Q*Q*Omega
q = 1;

rng(5);

P = min(2*N,m);
Omega = randn(m,P);
           
Y=zeros(m,P);    

% Y = Q*Omega;
disp('Computation of Y, Step 1/5 starting'); 
Y = runmexBBFMM3D(source,Omega,nCheb,L,level,ExecName,use_chebyshev,0);

% Y = (Q*Q)^q*(Q*Omega);
for d = 1 : 2*q
   disp(['Computation of Y, Step ',num2str(d+1),'/5 starting'])
   Y= runmexBBFMM3D(source,Y,nCheb,L,level,ExecName,use_chebyshev,0);    
end
% step 2
[R,~,~] =svd(Y,0); 

% step 3  B = R'*Q --> B' = Q'R = QR
disp('Computation of B, Last Step starting');
B = runmexBBFMM3D(source,R,nCheb,L,level,ExecName,use_chebyshev,0);

% step 4
[V,S,Ut] = svd(B,0);

% step 5
U = R*Ut;
U = U(:,1:N); S = S(1:N,1:N); V = V(:,1:N);



toc
%Plotting only works for imported grid
plotflag=true;
zlevel=mean(grid.z);
if plotflag
%Ploting the first three bases
figure; plotU(grid,U(:,1),zlevel)
title('Basis 1')
figure; plotU(grid,U(:,2),zlevel)
title('Basis 2')
figure; plotU(grid,U(:,3),zlevel)
title('Basis 3')

end
end
