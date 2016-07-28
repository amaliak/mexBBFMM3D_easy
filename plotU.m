function plotU(grid,U,zlevel)
% function for three dimensional plotting

% interpolate to regular grid
[X,Y,Z,E]=vec2mat(grid.x,grid.y,grid.z,U(:,1));

figure;
h=slice(X,Y,Z,E,[],[],zlevel);

set(h,'LineStyle','none');
colormap('hot')
colorbar('Location','SouthOutside')
title('- Basis')
xlabel(['Level',num2str(zlevel)])
set(h,'LineStyle','none');
hold on
%plotwells()
hold off
axis tight
daspect([1 1 1])
view(gca,[-54 24]);

end
