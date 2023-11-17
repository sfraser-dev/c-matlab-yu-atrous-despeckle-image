# C Mex MATLAB despeckling algorithm

The hole_wavelet2d.c mex file must be compiled and
available for "despeckleExample.m" and "yu_at4SPEK.m" to work.
To compile it, type the following at the MATLAB prompt:
>> mex hole_wavelet2d.c  


Run "despeckleExample.m" in MATLAB. This m-file creates a
speckle corrupted Lena image ("speckledImage") and passes
it to the novel denoising algorithm ("yu_at4SPEK.m" - so 
called as it is based upon Yu's algoritm, it utilises
the A Trous wavelet transform (4 levels) and it is set up
for despeckling images). The noisy and cleaned images are
then displayed.
