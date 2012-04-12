/* writeKernelToFile.c ---
 *
 * Description: Writes the kernel to LIBSVM compatible structure. For
 * a non-symmetric matrix, have to input K' instead of K.
 *
 * Status:
 * Author: Vinay Jethava
 * Created: Thu Nov 10 02:23:36 2011 (+0100)
 * Last-Updated: Mon Nov 14 17:22:27 2011 (+0100)
 *           By: Vinay Jethava
 *     Update #: 81
 */

/* Change Log:
 * 14-Nov-2011    Vinay Jethava  
 *    Last-Updated: Mon Nov 14 17:12:48 2011 (+0100) #69 (Vinay Jethava)
 *    Added transpose operation for converting K-> K' (for writing to file)
 * 14-Nov-2011    Vinay Jethava  
 *    Last-Updated: Mon Nov 14 17:12:48 2011 (+0100) #69 (Vinay Jethava)
 *    Added labels (as third argument) - optional 
 * 10-Nov-2011    Vinay Jethava
 *    Last-Updated: Thu Nov 10 03:46:50 2011 (+0100) #49 (Vinay Jethava)
 *    Changed name to svmWriteKernel
 * 10-Nov-2011    Vinay Jethava
 *    Last-Updated: Thu Nov 10 03:43:21 2011 (+0100) #45 (Vinay Jethava)
 *    Initial version - does not handle non-symmetric matrices
 *
 */

/* Code: */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include "mex.h"
#include <iostream.h>
#include <math.h>
#if MX_API_VER < 0x07030000
typedef int mwIndex;
#endif

#define CMD_LEN 2048
#define Malloc(type,n) (type *)malloc((n)*sizeof(type))

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {
  char *fileName;
  plhs=NULL;
  if ( (!mxIsDouble(prhs[0])) || (!mxIsChar(prhs[1])) ) {
    mexErrMsgTxt("Usage: svmWriteKernel(K, fileName, [labels])\n");
  }


  // mexWarnMsgTxt("svmWriteKernel() uses K' not K for symmetric matrices.");
  int l =  (int) mxGetN(prhs[1]);
  fileName= Malloc(char, l+1);
  mxGetString(prhs[1], fileName, l+1);
  FILE* fid = fopen(fileName, "w");
  int N = (int) mxGetN(prhs[0]);
  mxArray *ktranspose[1];
  mexCallMATLAB(1, ktranspose, 1, (mxArray**) &prhs[0], "transpose"); 
  double *data = (double*) mxGetPr(*ktranspose);
  int count = 0;
  double *label = (double*) mxGetPr(prhs[2]);
  for(int i=0; i < N; i++) {
    // print the initial id
    if(label!=NULL) {
      int a = (int) floor(*(label+i));
      fprintf(fid, "%d", a);
    } else {
      fprintf(fid, "1");
    }

    fprintf(fid, " 0:%d", (i+1));
    for(int j=0; j < N; j++) {
      fprintf(fid, " %d:%g", (j+1),  *(data+count));
      ++count;
    }
    fprintf(fid, "\n");
  }
  fclose(fid);
}


