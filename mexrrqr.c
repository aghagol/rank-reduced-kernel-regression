#include "mex.h"
#include <math.h>
#include "matrix.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    double *pkp     = mxGetPr(prhs[0]);
    double *g       = mxGetPr(prhs[1]);
    double *alpha   = mxGetPr(prhs[2]);
    double *thresh  = mxGetPr(prhs[3]);
    double *f       = mxGetPr(prhs[4]);
    double *xk      = mxGetPr(prhs[5]);
    double *yk      = mxGetPr(prhs[6]);
    double *u1      = mxGetPr(prhs[7]);
    double *u2      = mxGetPr(prhs[8]);
    double *sigma2  = mxGetPr(prhs[9]);
    double *e1      = mxGetPr(prhs[10]);
    double *e2      = mxGetPr(prhs[11]);
    int m = (int)mxGetM(prhs[1]);
    int n = (int)mxGetN(prhs[1]);
    double *z, *q;
    int pk = 9;
    int np = 81;
    double rb[81], gb[81], bb[81], d1[81], d2[81], Wk[81], Xg1[81], Xg2[81];
    double X11, X12, X22, R11, R12, R22, determinan;
    int i, j, ii, jj, counter;
    
    mwSize dims[3];
    dims[0] = (mwSize)m;
    dims[1] = (mwSize)n;
    dims[2] = 3;
    plhs[0] = mxCreateNumericArray(3, dims, mxDOUBLE_CLASS, mxREAL);
    plhs[1] = mxCreateNumericArray(3, dims, mxDOUBLE_CLASS, mxREAL);
    z = mxGetPr(plhs[0]);
    q = mxGetPr(plhs[1]);
    mexPrintf("%d %d \n",m,n);
    
    for(j=(pk-1)/2; j<(n-(pk-1)/2); j++)
    {
        for(i=(pk-1)/2; i<(m-(pk-1)/2); i++)
        {
            if (*(alpha + i + j*m) > *thresh)
            {
                mexPrintf("%d %d \n",i,j);
                counter = 0;
                for(jj=(1-pk)/2; jj<=(pk-1)/2; jj++)
                {
                    for(ii=(1-pk)/2; ii<=(pk-1)/2; ii++)
                    {
//                         mexPrintf("%d \n",(ii+(pk-1)/2) + (jj+(pk-1)/2)*pk);
                        rb[counter] = *(f + i+ii + (j+jj)*m + 0*m*n);
                        gb[counter] = *(f + i+ii + (j+jj)*m + 1*m*n);
                        bb[counter] = *(f + i+ii + (j+jj)*m + 2*m*n);
                        counter++;
                    }
                }
//                 for(ii=0; ii<np; ii++)
//                 {
//                     d1[ii] = pow(*(xk+ii) * *(u1 + i + j*n) + *(yk+ii) * *(u2 + i + j*n), 2);
//                     d2[ii] = pow(*(yk+ii) * *(u1 + i + j*n) - *(xk+ii) * *(u2 + i + j*n), 2);
//                     Wk[ii] = exp(-.5 / *sigma2 / *(alpha+i+j*n) * (*(e1+i+j*n) * d1[ii] + *(e2+i+j*n) * d2[ii]));
//                 }
//                 X11 = 0;
//                 X12 = 0;
//                 X22 = 0;
//                 for(ii=0; ii<np; ii++)
//                 {
//                     X11 += Wk[ii];
//                     X12 += Wk[ii] * d1[ii];
//                     X22 += Wk[ii] * d1[ii] * d1[ii];
//                 }
//                 mexPrintf("X11=%f, X12=%f, X22=%f\n", X11, X12, X22);
//                 determinan = X11*X22-X12*X12;
//                 R11 = X22 / (determinan + 0.0000001);
//                 R22 = X11 / (determinan + 0.0000001);
//                 R12 = -1 * X12 / (determinan + 0.0000001);
//                 for(ii=0; ii<np; ii++)
//                 {
//                     Xg1[ii] = R11 * Wk[ii] + R12 * Wk[ii] * d1[ii];
//                     Xg2[ii] = R12 * Wk[ii] + R22 * Wk[ii] * d1[ii];
//                 }
//                 for(ii=0; ii<np; ii++)
//                 {
//                     *(z + i + j*n + 0*m*n) = Xg1[ii] * rb[ii];
//                     *(z + i + j*n + 1*m*n) = Xg1[ii] * gb[ii];
//                     *(z + i + j*n + 2*m*n) = Xg1[ii] * bb[ii];
//                     *(q + i + j*n + 0*m*n) = Xg2[ii] * rb[ii];
//                     *(q + i + j*n + 1*m*n) = Xg2[ii] * gb[ii];
//                     *(q + i + j*n + 2*m*n) = Xg2[ii] * bb[ii];
//                 }
            }
            else
            {
                *(z + i + j*n + 0*m*n) = *(f + i + j*n + 0*m*n);
                *(z + i + j*n + 1*m*n) = *(f + i + j*n + 1*m*n);
                *(z + i + j*n + 2*m*n) = *(f + i + j*n + 2*m*n);
                *(q + i + j*n + 0*m*n) = 0;
                *(q + i + j*n + 1*m*n) = 0;
                *(q + i + j*n + 2*m*n) = 0;

            }
        }
    }
}