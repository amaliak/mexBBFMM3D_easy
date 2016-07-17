function QH = runmexBBFMM3D(coord,H,nCheb,L,level,execname,use_chebyshev,TestingMode)
% Step 2 of running mexBBFMM3D
% Requires setting up coordinates - see example1.m
% -----------
%Input:  
% coord :These are the coordinates of each  point of the 3D domain provided as vectors. 
%        The grid does not have to be regular or structured.
%        3-D locations, stored column-wise i.e. (x | y | z)
% H:     The matrix we want to multiply Q with (multiply Q with H)
% nCheb: Number of Chebyshev nodes (the higher nCheb, the higher the accuracy,
%        but higher the computational cost), integers between 4 to 6 are suggested.
%        Try starting with 4 and increase to achieve desirable accuracy.
% level: The higher the level the higher the accuracy but higher
%        computational cost. Try starting at 4  and increase to achieve desirable accuracy.
% L:     The length of the cube defining the domain. This can be calculated
%        from the coordinates, see example provided
% use_chebyshev: FMM parametet that should be set to 1. For more advanced
%           users, see documentation of BBFMM3D
% execname: Use name chosen during compilation in compilemex.m
% TestingMode: If TestingMode=1, the code multiplies Q and H with both fast (BBFMM2D) 
%              and direct method. It outputs the runnting time for each method as well 
%              as the relative error. This mode is useful when one wants to determine 
%              how many Chebyshev nodes (nCheb) to use for a desired accuracy. 
%              It will take a long time because the direct multiplication is also performed. 
%              In Application Mode (TestingMode=0), the code uses BBFMM2D only and does not 
%              compare with the direct approach. This mode should be used for large problems.
% -----------
% Output: product QH with BBFMM2D (application mode) 
%         product QH with BBFMM2D and QHexact with direct multiplication (testing mode)

if TestingMode
  %% testing mode 
  % Compute matrix-matrix product QH of dimension 10000x100
  eval(['[QH,QHexact] = ',execname,'(coord,coord,H,nCheb,level, L, use_chebyshev);']);
else
  %% application mode
  eval(['QH = ',execname,'(coord,coord,H,nCheb,level, L, use_chebyshev);']);
end
