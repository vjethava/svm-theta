/* svmReadModelFile.c --- 
 * 
 * Description: Reads libsvm model file into matlab - no error checking!!!
 * 
 * Status: NEEDS ERROR CHECKING
 * Author: Vinay Jethava
 * Created: Thu Nov 10 02:18:45 2011 (+0100)
 * Last-Updated: Wed Jan 18 19:24:27 2012 (+0100)
 *           By: Vinay Jethava
 *     Update #: 18
 */ 

/* Change Log:
 * 10-Nov-2011    Vinay Jethava  
 *    Last-Updated: Thu Nov 10 02:37:23 2011 (+0100) #14 (Vinay Jethava)
 *    Changed name to svmReadModel.c    
 * 10-Nov-2011    Vinay Jethava  
 *    Last-Updated: Thu Nov 10 02:18:55 2011 (+0100) #1 (Vinay Jethava)
 *    Working code for reading libsvm model into matlab 
 * 
 */
 
/* Code: */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include "svm.h"
#include "mex.h"
#include "svm_model_matlab.h"

#if MX_API_VER < 0x07030000
typedef int mwIndex;
#endif


#define CMD_LEN 2048
#define Malloc(type,n) (type *)malloc((n)*sizeof(type))

void print_null(const char *s) {}
void print_string_matlab(const char *s) {mexPrintf(s);}


static void fake_answer(mxArray *plhs[]){
  plhs[0] = mxCreateDoubleMatrix(0, 0, mxREAL);
}

 
void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {
	char *fileName = NULL; 
	if ((nrhs != 1)  || (!mxIsChar(prhs[0]))) {
	  mexErrMsgTxt("Please give model file name."); 
	  fake_answer(plhs);
	  return;
	}
	int l =  (int) mxGetN(prhs[0]); 
	fileName= Malloc(char, l+1); 
	mxGetString(prhs[0], fileName, l+1);
	// mexPrintf("Read model from file: %s\n" , fileName); 
	struct svm_model *model;
	model = svm_load_model(fileName);
	// mexPrintf("model #totalSV: %d\n", model->l);
	model_to_matlab_structure(plhs, 1, model); 
}
