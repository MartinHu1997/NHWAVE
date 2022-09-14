# NHWAVE
NHWAVE Model is Developed by Kirby, Ma, et. al.

And, I updated pieces of codes for user-friendliness.

Please refer to the document (./src/code_update_record (by HU).txt) for updating details.
Please refer to the those .txt files (in ./Examples of INPUT Files) if you would make input file.

Lihan Hu (Hohai Univ.),
Email: hulihan@hhu.edu.cn

================================================================

■ NOTICE 1: Debugs for Vegetation Model (Dec 7, 2021)

· A more reasonable parameter definition of 'Cfk' and 'Cfe' is applied, which can avoid the conflicts between global and local variables.

· To fix a bug caused by confusion of variable names. We replace the variable 'Vegbv' with 'StemD' at two places (in the subroutine kepsilon_3D of nhwave.F)

NOW, the Vegetable Model can be normally compiled!

Thanks to Y. Chen (Dalian University of Technology) who found these bugs first.
//
■ NOTICE 2: A new module, PRESSUREACC, has been added (Sept, 2022)

· This module can actually accelerate the computing process of poisson equation solving.

· The PDI method is real equipment. And the details of this method can be found in Shi et. al (2015).
REFERENCE: Shi, J., Shi, F., Kirby, J. T., Ma, G., Wu, G., Tong, C., & Zheng, J. (2015). Pressure Decimation and Interpolation (PDI) method for a baroclinic non-hydrostatic model. Ocean Modelling, 96, 265–279. https://doi.org/10.1016/j.ocemod.2015.09.010

· The pressure-subgrid method is real equipment now but not tested, which is more flexible than the PDI approach above.
//