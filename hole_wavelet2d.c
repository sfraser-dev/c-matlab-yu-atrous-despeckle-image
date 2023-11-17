/* a trous wavelet, 2 dimensional, length 3 and zero padded */
/* Based upon the kronecker function */

/* MATLAB USAGE: matrixOut=hole_wavelet2d(matrixIn,j_level) */

#include <math.h>
#include "mex.h"

/* prototype */
void hole_wavelet2d(double**,double**,int,int,int);
double fill_and_process(double**,int,int,int,int,int);

/* mex 'main' function */
void mexFunction(	int nlhs, 
			mxArray *plhs[], 
			int nrhs, 
			const mxArray *prhs[])	{
	int rows, cols;		/* number of rows and columns */
	int j, k, l, m, n;	/* for loop variables */
	int trv;		/* traverse mex arrays */
	int size_element;	/* size of array element (for malloc) */
	double *in_handle, *out_handle;	/* copying mex 'fortran' arrays */ 
	double *in_array, *out_array;	/* dynamic arrays (1-D) */
	double **in, **out;		/* dynamic arrays (quasi 2-D) */
	char *err_msg;			/* an error message string */ 
	int row_a2, col_a2;		/* checking argument two (j_level) */
	int j_level;			/* wavelet scale level */
	
	/* checking number of inputs */
	if(nrhs !=2)
		mexErrMsgTxt("Must have two input arguments");
	if(nlhs !=1)
		mexErrMsgTxt("Must have one output argument");
		
	/* preventing sparse, complex and string matrices */
	if( mxIsComplex(prhs[0])||!(mxIsClass(prhs[0],"double"))
	    || mxIsClass(prhs[0],"sparse") || mxIsChar(prhs[0]) ){
		err_msg="input must be real, double, full and non-string";
		mexErrMsgTxt(err_msg);
	}
				
				
	/* getting the number of rows and columns from input matrix */
	rows=mxGetM(prhs[0]);
	cols=mxGetN(prhs[0]);
	
	/* checking that argument two is a scalar */
	row_a2=mxGetM(prhs[1]);		/* should be 1 */
	col_a2=mxGetN(prhs[1]);		/* should be 1 */
	if((row_a2*col_a2)!=1)
		mexErrMsgTxt("Argument two must be scaler");
	
	/* setting the wavelet scale */
	j_level=(int)mxGetScalar(prhs[1]);
		
	/* creating an output array, giving it a handle */
	plhs[0]=mxCreateDoubleMatrix(rows,cols,mxREAL);
	out_handle=mxGetPr(plhs[0]);
	
	/* finding size of each element of input matrix */
	size_element=mxGetElementSize(prhs[0]);
	
	/* creating dynamic 2d arrays */
	in_array = (double*) mxMalloc (rows*cols*size_element);
	in = (double**) mxMalloc (rows*sizeof(double*));
	out_array = (double*) mxMalloc (rows*cols*size_element);
	out = (double**) mxMalloc (rows*sizeof(double*));
	for(j=0; j<rows; j++){
		in[j]=&(in_array[j*cols]);
		out[j]=&(out_array[j*cols]);
	}
	
	/* creating an input array, giving it a handle */
	in_handle=mxGetPr(prhs[0]);
	
	/*************************************************************************
	** copying 1-d input array into 2-d matrix, note that since MATLAB	**		
	** was originally written in fortran, it treats arrays in a manner 	**
	** similiar to fortran. Thus, in MATLAB, 2-d matrices are stored 	**
	** as 1-d matrices in the following way:				**
	**									**
	** (2-d)                   (1-d)					**
	** 1 2 3								**
	** 4 5 6	= 	1 4 2 5 3 6					**
	*************************************************************************/
	trv=0;	
	for(k=0; k<cols; k++){
		for(l=0; l<rows; l++)
			in[l][k]=in_handle[trv++];
	}
	
	/* PROCESS THE MATRIX HERE */
	hole_wavelet2d(in,out,rows,cols,j_level);
					
	/* copying the dynamic matrix to the output handle */
	trv=0;
	for(m=0; m<cols; m++){
		for(n=0; n<rows; n++)
			out_handle[trv++]=out[n][m];
	}
	
	/* freeing dynamic memory */
	mxFree(in[0]); in_array=0;
	mxFree(in); in=0;
	mxFree(out[0]); out_array=0;
	mxFree(out); out=0;
	
	mexUnlock();		/* allows for re-compiling */
}

void hole_wavelet2d(double** m_in,double** m_out,int no_rows,int no_cols,int J){
        int r,c;

        for(r=0; r<no_rows; r++)
          for(c=0; c<no_cols; c++)
            m_out[r][c]=fill_and_process(m_in,no_rows,no_cols,J,r,c);
}

double fill_and_process(double **img,int no_rows,int no_cols,int J,int r,int c){
        int x;                  /* used in filling kernel values */
        double val[3][3];       /* contains kernel values */
        float total=0.0;                /* used in wavelet */
        float mask_values[3][3];        /* used in wavelet */

        x = (int)pow(2,J);

        /* filling the "kernel with holes" elements */
        if(r-x < 0){
                val[0][2]=0;
                val[0][1]=0;
                val[0][0]=0;
        }
        else{
                if(c+x >= no_cols) val[0][2]=0;
                else val[0][2]=img[r-x][c+x];
                val[0][1]=img[r-x][c];
                if(c-x < 0) val[0][0]=0;
                else val[0][0]=img[r-x][c-x];
        }
        if(c-x < 0){
                val[1][0]=0;
                val[2][0]=0;
        }
        else{
                val[1][0]=img[r][c-x];
                if(r+x >= no_rows) val[2][0]=0;
                else val[2][0]=img[r+x][c-x];
        }
        if(r+x >= no_rows){
                val[2][1]=0;
                val[2][2]=0;
        }
        else{
                val[2][1]=img[r+x][c];
                if(c+x >= no_cols) val[2][2]=0;
                else val[2][2]=img[r+x][c+x];
        }
        if(c+x >= no_cols){
                val[1][2]=0;
        }
        else{
                val[1][2]=img[r][c+x];
        }
	
        val[1][1]=img[r][c];    /* centre value */


        /* processing the kernel */
        total=0.0;
        mask_values[0][0]=1.0/16;
        mask_values[0][1]=2.0/16;
        mask_values[0][2]=1.0/16;
        mask_values[1][0]=2.0/16;
        mask_values[1][1]=4.0/16;
        mask_values[1][2]=2.0/16;
        mask_values[2][0]=1.0/16;
        mask_values[2][1]=2.0/16;
        mask_values[2][2]=1.0/16;
        for(r=0; r<3; r++)
          for(c=0; c<3; c++)
            total+=(((float)val[r][c])*mask_values[r][c]);
        return (double)total;
}

