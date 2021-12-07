# NHWAVE
NHWAVE Model is Developed by Kirby, Ma, et. al.

And, I updated pieces of codes for user-friendliness.

Please refer to the document (./src/code_update_record (by HU).txt) for updating details.
Please refer to the those .txt files (in ./Examples of INPUT Files) if you would make input file.

Lihan Hu (Hohai Univ.),
Email: hulihan@hhu.edu.cn

================================================================
■ NOTICE: Debugs for Vegetation Model (Dec 7, 2021)
· More reasonable parameter definitions of 'Cfk' and 'Cfe' are applied, which can avoid the conflicts between global and local variables.
· To fix a bug caused by confusion of variable names. We replace the variable 'Vegbv' with 'StemD' at two places (in the subroutine kepsilon_3D of nhwave.F)

NOW, the Vegetable Model can be normally compiled!

Thanks to Y. Chen (Dalian University of Technology) who found these bugs first.
//