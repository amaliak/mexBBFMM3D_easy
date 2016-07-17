function [Xr,Yr,Zr,Cr]=vec2mat(xt,yt,zt,ct)

% this function takes the 3D irregular grid vector produced by TOUGH2 and 
% creates a 3D matrix of a regular grid 2x2x2 for plotting

%% new code for irregular grid
 % objective is to have a regular grid for plotting
 xr = [257:2:382];
 yr = [257:2:382];
 zr = [1019:-2:1001];
 [Xr,Yr,Zr] = meshgrid(xr,yr,zr);
 
 Cr = zeros(size(Xr));
 
 % source data
 Cr = griddata(xt,yt,zt,ct,Xr,Yr,Zr) ;