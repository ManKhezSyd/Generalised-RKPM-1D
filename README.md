# Generalised-RKPM-1D
The Generalised RKPM is an improved and more versatile version of the original RKPM, which was developed by:

Liu, Wing Kam, Sukky Jun, and Yi Fei Zhang. “Reproducing kernel particle methods.” International journal for numerical methods in fluids 20.8‐9 (1995): 1081-1106.

The generalised RKPM allows to include the derivatives of the field function in the kernel function. Thus, the DOFs also include derivatives of the variables. This makes Generalised RKPM a suitable tool for solving problems in which boundary conditions also include derivatives. 
 A few studies that are conducted using GenRKPM are as follows:
[*] Khezri, M., Gharib, M., Bradford, M. A., & Uy, B. (2015). A state space augmented generalised RKPM for three-dimensional analysis of thick and laminated composite plates. Computers & Structures, 158, 225-239.
[*] Khezri, M., M. Gharib, and K. J. R. Rasmussen. “A unified approach to meshless analysis of thin to moderately thick plates based on a shear-locking-free Mindlin theory formulation.” Thin-Walled Structures 124 (2018): 161-179.
[*] Gharib, M., Khezri, M., Foster, S. J., & Castel, A. (2017). Application of the meshless generalised RKPM to the transient advection-diffusion-reaction equation. Computers & Structures, 193, 172-186.
[*] Gharib, M., Khezri, M., & Foster, S. J. (2017). Meshless and analytical solutions to the time-dependent advection-diffusion-reaction equation with variable coefficients and boundary conditions. Applied Mathematical Modelling, 49, 220-242.
[*] Shodja, H. M., Khezri, M., Hashemian, A., & Behzadan, A. (2010). RKPM with augmented corrected collocation method for treatment of material discontinuities. Computer Modeling in Engineering & Sciences(CMES), 62(2), 171-204.

- The codes uploaded here are for one-dimensional shape functions of generalised RKPM. The codes are written using intel Fortran and are independent modules that should be included in a project. 
Don’t hesitate to contact me if you require further information….
I will gradually upload 2D and 3D shape functions too. 
