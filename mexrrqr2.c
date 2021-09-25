#include "mex.h"
#include <math.h>
#include "matrix.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    double *pkd     = mxGetPr(prhs[0]);
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
    int pk = (int)pkd[0];
    int np = pk * pk;
    double rb[121], gb[121], bb[121];
    double d1[121], d2[121], Wk[121], Xg1[121], Xg2[121];
    double X11, X12, X22, R11, R12, R22, det;
    int i, j, ii, jj, counter;
    
    mwSize dims[3];
    dims[0] = (mwSize)m;
    dims[1] = (mwSize)n;
    dims[2] = 3;
    plhs[0] = mxCreateNumericArray(3, dims, mxDOUBLE_CLASS, mxREAL);
    plhs[1] = mxCreateNumericArray(3, dims, mxDOUBLE_CLASS, mxREAL);
    z = mxGetPr(plhs[0]);
    q = mxGetPr(plhs[1]);

    for(j=(pk-1)/2; j<(n-(pk-1)/2); j++)
    {
        for(i=(pk-1)/2; i<(m-(pk-1)/2); i++)
        {
            if (alpha[i + j*m] > thresh[0])
            {
                counter = 0;
                for(jj=(1-pk)/2; jj<=(pk-1)/2; jj++)
                {
                    for(ii=(1-pk)/2; ii<=(pk-1)/2; ii++)
                    {
                        rb[counter] = f[i+ii + (j+jj)*m + 0*m*n];
                        gb[counter] = f[i+ii + (j+jj)*m + 1*m*n];
                        bb[counter] = f[i+ii + (j+jj)*m + 2*m*n];
                        counter++;
                    }
                }
                X11 = 0;
                X12 = 0;
                X22 = 0;
                for(ii=0; ii<np; ii++)
                {
                    d1[ii] = pow(xk[ii] * u1[i + j*m] + yk[ii] * u2[i + j*m], 2);
                    d2[ii] = pow(yk[ii] * u1[i + j*m] - xk[ii] * u2[i + j*m], 2);
                    Wk[ii] = exp(-.5 / sigma2[0] / alpha[i+j*m] * (e1[i+j*m] * d1[ii] + e2[i+j*m] * d2[ii]));
                    X11 += Wk[ii];
                    X12 += Wk[ii] * d1[ii];
                    X22 += Wk[ii] * d1[ii] * d1[ii];
                }
                det = X11*X22-X12*X12;
                R11 = X22 / det;
                R22 = X11 / det;
                R12 = -1 * X12 / det;
                for(ii=0; ii<np; ii++)
                {
                    Xg1[ii] = R11 * Wk[ii] + R12 * Wk[ii] * d1[ii];
                    Xg2[ii] = R12 * Wk[ii] + R22 * Wk[ii] * d1[ii];
                    z[i + j*m + 0*m*n] += Xg1[ii] * rb[ii];
                    z[i + j*m + 1*m*n] += Xg1[ii] * gb[ii];
                    z[i + j*m + 2*m*n] += Xg1[ii] * bb[ii];
                    q[i + j*m + 0*m*n] += Xg2[ii] * rb[ii];
                    q[i + j*m + 1*m*n] += Xg2[ii] * gb[ii];
                    q[i + j*m + 2*m*n] += Xg2[ii] * bb[ii];
                }
            }
            else
            {
                z[i + j*m + 0*m*n] = f[i + j*m + 0*m*n];
                z[i + j*m + 1*m*n] = f[i + j*m + 1*m*n];
                z[i + j*m + 2*m*n] = f[i + j*m + 2*m*n];
                q[i + j*m + 0*m*n] = 0;
                q[i + j*m + 1*m*n] = 0;
                q[i + j*m + 2*m*n] = 0;
            }
        }
    }
}