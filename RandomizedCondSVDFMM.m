function [U,S,V] = RandomizedCondSVDFMM(grid,Kernel,Corlength,Corlengthz,N,a)
tic
% Parallel code for low-rank decomposition of covariance kernel Q
% via randomized algorithm
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
%                load('./coord_htr2.mat')
%                grid.x = x_htr2; grid.y = y_htr2; grid.z = z_htr2;
%               [U,S,V] = RandomizedCondSVDFMM(grid,'GAUSSIAN',100,10,3,9);
% Input:
%        grid     : structure with vectors grid.x, grid.y, grid.z
%                each vector containing all x,y and z coordinates
%                respectively
%        m : size of covariance matrix number of unknowns
%        N : rank of reduced rank svd 
%        a : oversampling parameter for randSVD
%        q is 1, or 2 (hardcoded below to = 1)
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


% set up parallel pool
% matrix H is split in submatrices depenting on how many processors are
% available, and then the QH matrix is constucted by assumptling the
% submatrices i.e. QH = [QH1 QH2 ... QHn].

% find how many processors are available
disp('Setting up parfor')
% check if there is existing parallel job running
% if so, delete it
p = gcp('nocreate');
if ~isempty(p); delete(gcp); end
myCluster  = parcluster('local');
delete(myCluster.Jobs);

% start parallel job, find out number of workers available
poolobj = parpool('local');
noproc = poolobj.NumWorkers;
%noproc=12;
subH = (N+a)/noproc;
if mod(N+a,noproc)~=0
    disp(['N is ',num2str(N)])
    disp(['a is ',num2str(a)])
    disp(['noproc is ',num2str(noproc)])
    % decrease a to fit noproc
    % no need to regenerate Omega, just use fewer columns
    disp('Adjusting a in randSVD to fit number of processors') 
    a = floor(subH)*noproc - N; 
    if a < 1
        a = ceil(subH)*noproc-N;
    end
    if (N+a) < (N+1)
        display('Error in randomized svd');
        display('Problem with dimensions of Omega')
        keyboard;
    end
    disp(['New a is ',num2str(a),' and N+a is ',num2str(N+a)])		
    %a = input('Adjust ling a to fit noproc');
end

%subH is the number of submatrices depending on the number of processors
subH = (N+a)/noproc;

if mod(subH,1) > 0; keyboard; end;

% end of setting up parallel runs
% start of SVD


rng(5);
Omega = randn(m,N+a);

           
Y=zeros(m,N+a);    

% Y = Q*Omega;
disp('Computation of Y, Step 1/5 starting'); 
parfor i=1:noproc
    start = subH*(i-1)+1;
    fin = i*subH;
    
    Hi = Omega(:,start:fin);
    res(:,:,i) = runmexBBFMM3D(source,Hi,nCheb,L,level,ExecName,use_chebyshev,0);
    disp(['Multiplication of colums from ',num2str(start),' to ',num2str(fin)])
end

for i=1:noproc
    start=subH*(i-1)+1;
    fin=i*subH;
    Y(:,start:fin)=res(:,:,i);
end

% Y = (Q*Q)^q*(Q*Omega);
for d = 1 : 2*q
   disp(['Computation of Y, Step ',num2str(d+1),'/5 starting'])
   parfor i=1:noproc
    	start = subH*(i-1)+1;
    	fin = i*subH;
        Hi = Y(:,start:fin);
        res(:,:,i) = runmexBBFMM3D(source,Hi,nCheb,L,level,ExecName,use_chebyshev,0);
    	disp(['Multiplication of colums from ',num2str(start),' to ',num2str(fin)])    
    end
    
    for i=1:noproc
    	start = subH*(i-1)+1;
     	fin = i*subH;
     	Y(:,start:fin) = res(:,:,i);
    end    
end
% step 2
[R,~,~] =svd(Y,0); 

% step 3  B = R'*Q --> B' = Q'R = QR
disp('Computation of B, Last Step starting');
parfor i=1:noproc
        start = subH*(i-1)+1;
        fin = i*subH;
        Hi = R(:,start:fin);
        res(:,:,i) = runmexBBFMM3D(source,Hi,nCheb,L,level,ExecName,use_chebyshev,0);
	disp(['Multiplication of colums from ',num2str(start),' to ',num2str(fin)])
end

disp('Computation of B, Last Step finished');
for i=1:noproc
    start=subH*(i-1)+1;
    fin=i*subH;
    B(:,start:fin)=res(:,:,i);
end

% step 4
[V,S,Ut] = svd(B,0);

% step 5
U = R*Ut;
U = U(:,1:N); S = S(1:N,1:N); V = V(:,1:N);

% close parallel pool
p = gcp('nocreate');
if ~isempty(p)
    delete(gcp)
end


toc
plotflag=false;
if plotflag
%Ploting the first three bases
figure; plotU(grid,U(:,1),1005)
title('Basis 1')
figure; plotU(grid,U(:,2),1005)
title('Basis 2')
figure; plotU(grid,U(:,3),1005)
title('Basis 3')

end
end
