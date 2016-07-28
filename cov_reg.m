function [Q,QH]=cov_reg(grid,kerneltype,lx,ly,lz,H)
% grid.x = -24:4:24; grid.y = -24:4:24; grid.z = -9:3:9;
% function to calculate covariance for regular grid using vectorized form
% to check how slower my code is

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


QH = Q*H;

end