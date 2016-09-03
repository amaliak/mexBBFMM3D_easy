function [U,S,V,Size,TimeElapsed,TestError] = rsvd(A,K,TestingMode)
%-------------------------------------------------------------------------------------
% random SVD
% Extremely fast computation of the truncated Singular Value Decomposition, using
% randomized algorithms as described in Halko et al. 'finding structure with randomness
%
% usage : 
%
%  input:
%  * A : matrix whose SVD we want
%  * K : number of components to keep
%
%  output:
%  * U,S,V : classical output as the builtin svd matlab function
%-------------------------------------------------------------------------------------
% Antoine Liutkus  (c) Inria 2014
tStart = tic;
datestr(now) 
[M,N] = size(A);
P = min(2*K,N);
X = randn(N,P);
Y = A*X;
W1 = orth(Y);
B = W1'*A;
[W2,S,V] = svd(B,'econ');
U = W1*W2;
K=min(K,size(U,2));
UK = U(:,1:K);
SK = S(1:K,1:K);
VK=V(:,1:K);
datestr(now) 
TimeElapsed = toc(tStart);  % TOC, pair 2  
Size=N;

if TestingMode
    TestError = norm(UK*SK*VK'-A)/norm(A);
    disp('The error is:')
    disp([num2str(100*TestError),' %'])
end

