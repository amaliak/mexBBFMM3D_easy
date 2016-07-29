function [Q,QH]=cov_reg(grid,kerneltype,lx,ly,lz,H)
% Function that creates covariance matrix and optionally 
% performs direct multiplication with matrix H
% Example use for creating a covariance matrix:
%           grid.x = -24:4:24; grid.y = -24:4:24; grid.z = -9:3:9;
%           [Q,~]=cov_reg(grid,'GAUSSIAN',6,6,6,[]); imagesc(Q);
% Example use for creating a covariance matrix and multiplying with H:
%           grid.x = -24:4:24; grid.y = -24:4:24; grid.z = -9:3:9;
%           [Q,QH]=cov_reg(grid,'GAUSSIAN',6,6,6,H); imagesc(Q);
%
% Input: 
%      grid: structure with x,y and z coordinates
%      kerneltype: choose from 'GAUSSIAN', 'LINEAR' , 'ONEOVERR2'
%      lx,ly,lz: correlation lengths in x,y and z
%      H: matrix for multiplication Q * H
% Output: 
%      Q: Covariance matrix Q
%      QH: Q * H

[x,y,z] = meshgrid(grid.x,grid.y,grid.z);
[xj,xl]=meshgrid(x(:),x(:));
[yj,yl]=meshgrid(y(:),y(:));
[zj,zl]=meshgrid(z(:),z(:));

rx = (xj-xl).^2.*(1.0/lx)*(1.0/lx);
ry = (yj-yl).^2.*(1.0/ly)*(1.0/ly);
rz = (zj-zl).^2.*(1.0/lz)*(1.0/lz);
r = sqrt( rx + ry + rz );

switch kerneltype
    case 'GAUSSIAN'
    Q = exp(-r.*r);               
    case 'LINEAR'
    Q = r;
    case 'ONEOVERR2'
    Q = -1./(r.^2);
    otherwise 
    disp('Create own kernel')
end

if ~isempty(H)
    % check for dimension agreement
    m = size(Q);
    n = size(H,1);
    if (m-n)
        disp('Dimension mismatch please correct H')
        disp(['m= ',num2str(m)])
        disp(['n= ',num2str(n)])
        QH = NaN;
    else
        QH = Q*H;
    end
else
    QH = NaN;
end
    

end