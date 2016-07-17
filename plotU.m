function plotU(U,level)
% example usage

% first call 
% eSol = qreadSAVE(); to create structure with solution 
% and then
% plotX(eSol,1010,'pressure')

% variable can be pressure, T (for temperature) or X (for concentration)

% G.z, G.elem
load MESHmap_tSols_pump_small.mat

% Get permeability field from MESH file or load it from somewhere
%G = qreadMESH('./',[]); %G.pmx has the permeability
%[X,Y,Z,K]=vec2mat(xt,yt,zt,G.pmx(1:end));

% interpolate to regular grid
[X,Y,Z,E]=vec2mat(xt,yt,zt,U(:,1)');

figure;
h=slice(X,Y,Z,E,[],[],level);

set(h,'LineStyle','none');
colormap('hot')
colorbar('Location','SouthOutside')
title('- Basis')
xlabel(['Level',num2str(level)])
set(h,'LineStyle','none');
hold on
%plotwells()
hold off
axis tight
daspect([1 1 1])
view(gca,[-54 24]);

end
