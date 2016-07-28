function D=cov_irg(varargin)

% function to calculate covariance matrix of irregular grid

% x is the VECTOR of all x's 
% y is the VECTOR of all y's
% e.g. for a 2x2 grid  x=[1 1 2 2] y=[1 2 1 2]

% if grid is regular this can be done by: 
% [x,y]=meshgrid(x,y);
% [xj,xl]=meshgrid(x(:),x(:));
% [yj,yl]=meshgrid(y(:),y(:));
% D=cov_irg(xl(:,1),yl(:,1));

if nargin == 2
    % 2D problem
    x=varargin{1};
    y=varargin{2};
    z=zeros(size(y));
elseif nargin == 3
    % 3D problem
    x=varargin{1};
    y=varargin{2};
    z=varargin{3};
end

m=length(x);
if (m-length(y)) || (m-length(z))
    disp('x, y and z need to be the same size')
end

D=zeros(m,m);
tic

for k=1:m
    for l=1:m
        D(k,l)=sqrt((x(k)-x(l)).^2 + (y(k)-y(l)).^2 + (z(k)-z(l)).^2);
    end
end

toc
end
