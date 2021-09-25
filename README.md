To use the faster C-compiled code, you need to compile the mexrrqr2.c using MATLAB MEX compiler.

Start by running the main_ICASSP13.m, after correcting the paths inside the m-file. This main file calls the rrqr1 function and rrqr1 calls the mexrrqr2 function which is also optionally implemented in MATLAB (inside rrqr2.m). However, the MATLAB implementation in rrqr2.m is much slower than the mex version. It is provided in case you cannot compile mexrrqr2.c and the goal is just to test the code on a small dataset.

There is also a rrqr_vectorized function which implements rrqr1+rrqr2 in vectorized MATLAB but that is still slower than rrqr1+mexrrqr2. Again, this is a good option if you cannot get mexrrqr2 to work.

Attached code also includes some benchmark (SKR code) that was used to produce the paper results.