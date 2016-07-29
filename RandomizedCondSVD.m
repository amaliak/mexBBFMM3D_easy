function [UN,SN,VN] = RandomizedCondSVD(A,N,q,TestingMode,CompareMode)
% Function that performs truncated SVD by a randomized algorithm
% based on Halko et al., 2009
% Example usage:
%       grid.x = -12:6:12; grid.y = -12:6:12; grid.z = -6:3:6;
%       [Q,~]=cov_reg(grid,'GAUSSIAN',6,6,6,[]); imagesc(Q);
%       [UN,SN,VN] = RandomizedCondSVD(Q,10,1,1);
% Input:
%       A:           Covariance matrix. Can be created by cov_reg.m or cov_irg.m
%       N:           Number of components of SVD needed, rank of reduced
%                    rank svd
%       q:           default is 1, or 2 (try 1)
%       Testingmode: if 1, the error compared to full svd will be
%                    calculated. Caution, do not use for very large matrices as it will
%                    take a very long time to perform the full svd
%       CompareMode: if 1, the times for different methods will be compared
%                    to evaluate efficiency: randomized Svd, full svd, matlab's svd(Q,0)
%                    and matlab's svd(Q,'econ')
% Output: 
%       UN: First N left singular vectors, returned as the columns of a
%           matrix
%       SN: First N singular values, returned as a column vector. 
%           The singular values are nonnegative real numbers listed in decreasing order.
%       VN: First N right singular vectors, returned as the columns of a
%           matrix

tic
[m,n] = size(A);
if N>m || N>n
    warning('N set to minimum of n and m')
    N = min(n,m);
end

Omega = randn(n,N+10);
Y = (A*A')^q*A*Omega;
[Q,~,~] =svd(Y,0);

B = Q'*A;

[V,S,Ut] = svd(B',0);
U = Q*Ut;
UN = U(:,1:N); SN = S(1:N,1:N); VN = V(:,1:N);

disp('For randomized svd ->')
toc

if TestingMode
    TestError = norm(UN*SN*VN'-A)/norm(A);
    disp('The error is:')
    disp([num2str(100*TestError),' %'])
end

SN=diag(SN);

% to compare times
if CompareMode

    tic
    svd(A);
    disp('For full svd ->')
    toc

%     tic
%     svd(A,0);
%     disp('For Matlab truncated svd->')
%     toc
% 
%     tic
%     svd(A,'econ');
%     disp('For Matlab truncated svd->')
%     toc
    
end
