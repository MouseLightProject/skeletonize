
#include <math.h>
#include <stdio.h>
#include "mex.h"
       

void mexFunction(int nlhs, mxArray *plhs[], int	nrhs, const mxArray *prhs[]) 
{
	int       rows, cols, slices;
	double    *voxel;
        double    *DisMatrix;
        double    *ridge;
        double    *temp=(double *)malloc(27*sizeof(double));
        int       th;
        int       i,j,k;
        double    left,right,ahead,rear,up,down,center;

	/* reading the parameters */
        voxel = (double *) mxGetPr (prhs [0]);
        DisMatrix = (double *) mxGetPr (prhs [1]);
	rows = (int)mxGetScalar (prhs [2]);
	cols = (int)mxGetScalar (prhs [3]);
	slices = (int)mxGetScalar (prhs [4]);
        th = (int)mxGetScalar (prhs [5]);

	/* require memory for return */
	plhs[0] = mxCreateDoubleMatrix(rows*cols*slices, 1, mxREAL);
	ridge = (double*) mxGetPr (plhs [0]);	
	
        for(i=0; i<rows; i++)
		for(j=0; j<cols; j++)
			for(k=0; k<slices; k++)
                        {
                                *(ridge+i+j*rows+k*cols*rows)=0;
                        }


        for(i=1; i<(rows-1); i++)
		for(j=1; j<(cols-1); j++)
			for(k=1; k<(slices-1); k++)
			{
                            left=*(DisMatrix+(i-1)+j*rows+k*cols*rows);
                            right=*(DisMatrix+(i+1)+j*rows+k*cols*rows);   
                            ahead=*(DisMatrix+i+(j-1)*rows+k*cols*rows);
                            rear=*(DisMatrix+i+(j+1)*rows+k*cols*rows);
                            up=*(DisMatrix+i+j*rows+(k-1)*cols*rows);
                            down=*(DisMatrix+i+j*rows+(k+1)*cols*rows);
                            center=*(DisMatrix+i+j*rows+k*cols*rows);

                            if(center>left && center>right && center>ahead && center>rear && center>up && center>down && center>=th) 
                                 *(ridge+i+j*rows+k*cols*rows)=255;
                                
			}


}