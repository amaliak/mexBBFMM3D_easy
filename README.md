##Quick Start Guide for mexBBFMM3D, a MATLAB interface for <a href="http://www.sciencedirect.com/science/ article/pii/S0021999109004665"> Black Box Fast Multipole Method (BBFMM3D)  

----------------
###Function
BBFMM3D: Perform fast (linear) multiplication of a kernel matix Q with a vector or matrix: P = QH
randSVD: Performs approximate singular value decomposition 
randSVD with BBFMM3D: Performs approximate singular value decomposition for large covariance kernels Q. BBFMM3D is used for fast multiplication of the Q matrix with random vectors within the randSVD code. 

###Overview
The Fast Multipole Method (FMM) is an algorithm that performs fast multiplication of an N x N dense matrix Q(x,y) where N is the number of unknown values at points (x,y) in a 2D domain, with a matrix H of size N x m (N>>m). The direct multiplication approach has complexity O(N^2), while the BBFMM3D has linear complexity O(N). The table below shows computation times on a single core CPU when using the BBFMM3D and when using direct multiplication.

|   N      |  Time for BBFMM3D (sec)  |  Time for direct multiplication (sec) |   
| -------: |:------------------------:|:-------------------------------------:| 
| 12,500   |                       4.3|  22.8                                 |
| 50,000   |                      5.27|  367.9                                |
| 100,000  |                       9.6|  1505.3                               |


The BBFMM3D algorithm uses Chebyshev interpolation to approximate the multiplication. The approximation error can be controlled by the number of Chebyshev nodes used, which can be adjusted to achieve the desired accuracy. mexBBFMM3D provides a Matlab interface for the BBFMM3D package, which is written in C++. Note that the C++ code is faster than the Matlab interface. The corresponding code for 2D cases can be found [here](https://github.com/judithyueli/mexBBFMM2D).


###Disclaimer

This is a quick-start guide with instructions on how to set up and use mexBBFMM3D in MATLAB, with two example m-file that can be used to perform matrix-vector and matrix-matrix multiplication for regular and irregular grids. A reasonably good knowledge of MATLAB is assumed and minimal understanding of the FMM theory is needed.  

For a more involved description of the code and the method please see [here](https://github.com/ruoxi-wang/BBFMM3D), and for a full description of the algorithm see [Fong and Darve 2006] in section [__Reference__](#ref).

In this guide, we will use the example of the multiplication of a Gaussian covariance matrix Q (termed as Gaussian kernel) with a matrix H. The method can also be applied for other smooth kernels (see section [__Appendix__](#ref_app)).

Please cite the following paper if you use this code:

Fong, William, and Eric Darve. "The black-box fast multipole methodshod." Journal of Computational Physics 228, no. 23 (2009): 8712-8725. You can see details <a href="http://www.sciencedirect.com/science/article/pii/S0021999109004665">here</a>.

###2. DIRECTORIES AND FILES


	./example1.m		:	Example of how to use mexBBFMM3D for isotropic, regular grid
        ./example2.m		:	Example of how to use mexBBFMM3D for anisotropic, irregular large scale grid 
	./make.m		:	Makefile 
	./include/		:	Relevant header files  
	./mexFMM3D.cpp	        :	mex function  
	./eigen/		:	Eigen Library  
	./BBFMM3D/		: 	BBFMM3D library
	./README.md		:	This file  
	./Troubleshooting.md    :       Instructions for troubleshooting compilation problems
	
###Quick start guide

####Step 1:  Download the code and supporting software

Supporting software includes:
- The fft library (http://www.fftw.org/fftw-2.1.5.tar.gz)
- Matlab SDK files. If you get errors related to SDK files, see Troubleshooting.md 

Attention Mac users: If you have upgraded to Xcode 7, see [Answer](http://www.mathworks.com/matlabcentral/answers/246507-why-can-t-mex-find-a-supported-compiler-in-matlab-r2015b-after-i-upgraded-to-xcode-7-0) by MathWorks Support Team on 28 Dec 2015.

####Step 1a:  Check if you have MEX and MATLAB Symbolic Math Toolbox set up in Matlab

This package relies on MATLAB MEX functions and MATLAB Symbolic Math Toolbox. In order to use MEX functions, you should setup mex.

Setup MEX by typing the following MATLAB command

```
      mex -setup  
```

Once mex is set up successfully, to ensure that  MEX is installed, try the following commands:
```
copyfile([matlabroot,'/extern/examples/mex/arraySize.c'],'./','f')     
mex -v arraySize.c     
arraySize(5000) 
```

The last two lines of the screen output should read  

```
Dimensions: 5000x5000
Size of  array in kilobytes: 24414
```

If you have trouble with setting up mex, see __Trouble Shooting.md__

####Step 1b:  Download the code from https://github.com/amaliak/randSVD-with-BBFMM3D

A folder called `randSVD-with-BBFMM3D-master` will be downloaded. Copy the folder `randSVD-with-BBFMM3D-master` to your specific working directory. In MATLAB, set the current folder to be the directory where the downloaded folder was copied. You will see the folders
`BBFMM3D\` and `Eigen\`, as well as an m-file called (`compilemex.m`) and two example m-files (`example1.m` and `example2.m`)  that we will use in this quick start quide. These m-files will be used to compile, set-up, test and use the `BBFMM3D`. The user only needs to change the input to the example m-files. No modifications will be needed to other m-files or the contents of the folder `BBFMM3D` which includes the c++ source code. 

Note: For Step 2 you have to operate in the main directory of mexBBFMM3D, which contains `make.m`, For Step 3 or 4, you can call the generated MEX-files (e.g., `ExecName.mexmaci64`) by moving them to your own working directory, or add the main directory of mexBBFMM3D to the path.

####Step 2: Compile the MEX file to make sure compilation works

__A.__ Open file `compilemex.m` and read instructions.

__B.__ Choose your kernel function type, e.g. `GAUSSIAN` and the correlation length and give your BBFMM3D case a name, e.g. `Test1`

Compile by giving the command: 
```
compilemex(ExecName,Kernel,corlength)
```

This will compile the source code and generate a MEX-file with the name you provided (e.g. `ExecName.mexmaci64`). The extension (`.mexmaci64`) will depend on your platform.

NOTE!! Recompile the MEX-file (step 2) when the kernel function is changed.

If compilation is successful, you should see the message `mex compiling is successful!`


####Step 3a: Run example1.m for a regular grid

__A.__ Open file `example1.m` and read instructions. The example recompiles the code and then computes QH, unless TestingMode=0.

__B.__ Choose input variables for `QH = example1(ExecName,grid,Kernel,corlength,H,TestingMode)`. Avoid using a very large grid while in Testing Mode, as the code performs the direct multiplication for comparison and it may take a very long time to be completed. The grid must be provided as unique x,y and z locations and (x,y,z) triplets will be created automatically.


Input:
    - `ExecName` : the name of the mexfile for the Kernel chosen
    
    - `grid`     : structure with vectors grid.x, grid.y, grid.z each vector containing x,y and z coordinates respectively
    
    - `Kernel`   : covariance type, e.g. 'GAUSSIAN'
    
    - `corlength`: correlation length, isotropic anisotropy in z direction supported, see code
    
    - `H`       : matrix by which Kernel is multiplied
    
    - `TestingMode`: if set to 1, BBFMM is recompiled and runs in TestingMode in order to determine parameters (nCheb) for desired accuracy. if set to 0, the direct multiplication is not performed and the last compiled executable is used for BBFMM3D.

 Output:
    - `ExecName.mexmaci64`: executable for mex file for given configuration
    
    - `QH`       : Product of Kernel chosen by matrix H specified in input

Example usage: 
```
	      grid.x = -62:4:62; grid.y = -62:4:62; grid.z = -9:3:9;
              QH = example1('TESTNAME',grid,'GAUSSIAN',50,ones(7168,1),1)
```
When run in TestingMode (`TestingMode = 1`), the output will give a relative error that compares the accuracy of BBFMM3D with the direct multiplication of Q*H. Example printout:

```
Starting FMM computation...

Pre-computation time: 1.9583
FMM computing time: 4.6788
FMM total time: 6.6370

Starting direct computation...
Direct calculation starts from: 0 to 0.
Exact computing time: 8.2686
Relative Error: 9.724730e-05
```

If the Relative Error is deemed low enough, the code can be run with `TestingMode=0`, in which case the direct multiplication is not performed for comparison. Note that the code will use the last compiled executable, even if otherwise specified in the input variables of the example. 

Example printout when `TestingMode = 0`: 
```
ExecName, Kernel, Corlength ignored
Using executable last compiled
------------------------------
Executable Name is TESTNAME
Kernel type is GAUSSIAN
Correlation length in x,y,z is 30

Starting FMM computation...

Pre-computation time: 0.0604
FMM computing time: 3.9679
FMM total time: 4.0283
```
####Step 3b: Run example2.m for an irregular grid with anisotropy

__A.__ Open file `example2.m` and read instructions. The example recompiles the code and then computes QH, unless TestingMode=0.

__B.__ Choose input variables for `QH = example2(ExecName,grid,Kernel,corlength,H,TestingMode)`. Avoid using a very large grid while in Testing Mode, as the code performs the direct multiplication for comparison and it may take a very long time to be completed. Q is defined on an IRREGULAR grid with the option for anisotropy in the z direction. The grid must be provided as (x,y,z) triplets. 

The input is as in `example1.m`, with the exception that the grid shoudl now contain all x,y and z triplets, so that each is an `Nx1` vector.

Input: 
     - `ExecName` : the name of the mexfile for the Kernel chosen
     
     - `grid`     : structure with vectors grid.x, grid.y, grid.z each vector containing all x,y and z triplets 
     
     - `Kernel`   : covariance type, e.g. 'GAUSSIAN'
     
     - `corlength`: correlation length in x and y, isotropic 
     
     - `corlengthz`: correlation length in z 

     - `H`        : matrix by which Kernel is multiplied
     
     - `TestingMode`: if set to 1, BBFMM is recompiled and runs in TestingMode in order to determine parameters (nCheb) for desired accuracy. 

 Output:
       - `ExecName.mexmaci64`: executable for mex file for given configuration
       
       - `QH`       : Product of Kernel chosen by matrix H specified in input

Example usage: 

```
	        load('./coord_htr.mat')
                grid.x = x_htr; grid.y = y_htr; grid.z = z_htr;
                QH = example1('TESTNAME',grid,'GAUSSIAN',50,10,ones(23910,1),1)
```

When run in TestingMode (TestingMode = 1), the output will give a relative error that compares the accuracy of BBFMM3D with the direct multiplication of Q*H. The screen printout is the same as in example1 above. 

### APPENDIX<a name="ref_app"></a>

__Kernel Options__

`r` : distance between two points

`L` : length scale parameter

`\sigma^2`: variance parameter

#### Example of kernel type:
+ Gaussian kernel 

      ![guasskernel](http://latex.codecogs.com/gif.latex?%5Cdpi%7B150%7D%20Q%28r%29%20%3D%20%5Csigma%5E2%20%5Cexp%28-%5Cdfrac%7Br%5E2%7D%7BL%5E2%7D%29)

+ Exponential kernel

      ![expkernel](http://latex.codecogs.com/gif.latex?%5Cdpi%7B150%7D%20Q%28r%29%20%3D%20%5Cexp%28-%5Cdfrac%7Br%7D%7BL%7D%29)

+ Logrithm kernel

      ![logkernel](http://latex.codecogs.com/gif.latex?%5Cdpi%7B150%7D%20Q%28r%29%20%3D%20A%20%5Clog%28r%29%2C%20A%3E0)

+ Linear kernel

      ![linearkernel](http://latex.codecogs.com/gif.latex?%5Cdpi%7B150%7D%20Q%28r%29%20%3D%20%5Ctheta%20r%2C%20%5Ctheta%20%3E0)

+ Power kernel

      ![powerkernel](http://latex.codecogs.com/gif.latex?%5Cdpi%7B150%7D%20Q%28r%29%20%3D%20%5Ctheta%20r%5Es%2C%20%5Ctheta%20%3E0%2C%200%20%3Cs%20%3C2)
        

#### This package uses:

1. [Eigen](http://eigen.tuxfamily.org/index.php?title=Main_Page)

2. [BBFMM2D](https://github.com/sivaramambikasaran/BBFMM2D)

#### Reference:<a name="ref"></a>
1. Sivaram Ambikasaran, Judith Yue Li, Peter K. Kitanidis, Eric Darve, Large-scale stochastic linear inversion using hierarchical matrices, Computational Geosciences, December 2013, Volume 17, Issue 6, pp 913-927 [link](http://link.springer.com/article/10.1007%2Fs10596-013-9364-0)
2. Judith Yue Li, Sivaram Ambikasaran, Eric F. Darve, Peter K. Kitanidis, A Kalman filter powered by H2-matrices for quasi-continuous data assimilation problems [link](https://www.dropbox.com/s/xxjdvixq7py4bhp/HiKF.pdf)
3. Saibaba, A., S. Ambikasaran, J. Li, P. Kitanidis, and E. Darve (2012), Application of hierarchical matrices to linear inverse problems in geostatistics, OGST Revue d’IFP Energies Nouvelles, 67(5), 857–875, doi:http://dx.doi.org/10.2516/ogst/2012064. [link](http://ogst.ifpenergiesnouvelles.fr/articles/ogst/abs/2012/05/ogst120061/ogst120061.html)
4. Fong, W., and E. Darve (2009), The black-box fast multipole method, Journal of Computational Physics, 228(23), 8712–8725.[link](http://www.logos.t.u-tokyo.ac.jp/~tau/Darve_bbfmm_2009.pdf)

<script type="text/javascript"
   src="http://cdn.mathjax.org/mathjax/latest/MathJax.js?config=TeX-AMS-MML_HTMLorMML"></script>







