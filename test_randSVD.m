function [Size,TimeElapsed,TestError]=test_randSVD(dx,method)

grid.x = -24:dx:24; grid.y = -24:dx:24; grid.z = -9:3:9;
[Q,~]=cov_reg(grid,'GAUSSIAN',6,6,6,[]);

compression = 100; 
if method == 1 
[UN,SN,VN,Size,TimeElapsed,TestError] = rsvd(Q,compression,1);
else
[UN,SN,VN,Size,TimeElapsed,TestError] = RandomizedCondSVD(Q,compression,1,0);
end