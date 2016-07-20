function compilemex(execname,kerneltype,corlength)
% This m-file compiles the mex code for new kernel
% Step 1 of setting up mexBBFMM3D
% This step should be repeated for every new kernel
% Example usage: 
%        compilemex('case1','GAUSSIAN',100)
%
% Input: 
%        execname:  Name of mex-file generated after compilation
%               Example: if execname = 'case1_test'; the file case1_test.mexmaci64
%               will be generated on a mac computer
%               The compilation produces the file ./BBFMM3D/include/kernelfun.hpp
%        kerneltype: covariance function
%               Options: GAUSSIAN, LINEAR,ONEOVERR2
%               more kernels can be added in the if statment below
%               for kernels supported see https://github.com/sivaramambikasaran/BBFMM2D
%        corlength: correlation length - isotropic. Anisotropy can be
%        included by scaling, see example2.m
%        homogen:  K(ax, ay) = a^m K(x,y),=> homogen = m
%        symmetry = 1;    symmetric: 1; non-symmetric: 0; % anti-symmetric: 0
%
% Output: Running compilemex.m will compile the source code and will
%         generate a MEX-file with name as given in outputfile (e.g. case1.mexmaci64). 
%         The extension (.mexmaci64) will depend on your platform
clearvars -except execname kerneltype corlength
delete('./*.o');
delete('./BBFMM3D/output/*.bin'); % make sure to delete output file if kernel is changed
syms r;                           % distance of two points (radius)

switch kerneltype
    case 'GAUSSIAN'
    kernel = exp(-r.*r);               
    homogen = 0;
    symmetry = 1;
    case 'LINEAR'
    kernel = r;
    homogen = 1;
    symmetry = 1;
    case 'ONEOVERR2'
    kernel = -1./(r.^2);
    homogen = -2;
    symmetry = 1;
    otherwise 
    disp('Create own kernel')
end

                                
make(r,kernel,corlength,homogen,symmetry,execname);




