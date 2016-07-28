function [Xr,Yr,Zr,Cr]=vec2mat(xt,yt,zt,ct)

% this function takes the 3D irregular grid vector produced by TOUGH2 and 
% creates a 3D matrix of a regular grid 2x2x2 for plotting

%% new code for irregular grid
 % objective is to have a regular grid for plotting
 xr = [min(xt):2:max(xt)];
 yr = [min(yt):2:max(yt)];
 zr = [min(zt):2:max(zt)];
 [Xr,Yr,Zr] = meshgrid(xr,yr,zr);
 
 Cr = zeros(size(Xr));
 
 % source data
 Cr = griddata(xt,yt,zt,ct,Xr,Yr,Zr) ;