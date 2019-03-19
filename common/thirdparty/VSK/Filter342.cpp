#include <math.h>
#include <stdio.h>
#include "mex.h"

double Max;
int total;

void compare(double *a, int m)
{
    int j;
    Max=a[0];    

    for(j=1;j<=m;j++)
    {
       if(Max<a[j])
       { 
          Max=a[j];
       }
       
    }
}


void count(double *a, int m)
{
    int j;
    total=0;
    for(j=1;j<=m;j++)
    {
       if(*(a+j)==Max)
       { 
          total++;
       }
       
    }
}

       

void mexFunction(int nlhs, mxArray *plhs[], int	nrhs, const mxArray *prhs[]) 
{
	int       rows, cols, slices;
	double    *voxel;
        double    *DisMatrix;
        double    *ridge;
        double    *temp=(double *)malloc(343*sizeof(double));
        int       th, contain;
        int       i,j,k;
        int       a,b,c;

	/* reading the parameters */
        voxel = (double *) mxGetPr (prhs [0]);
        DisMatrix = (double *) mxGetPr (prhs [1]);
	rows = (int)mxGetScalar (prhs [2]);
	cols = (int)mxGetScalar (prhs [3]);
	slices = (int)mxGetScalar (prhs [4]);
        th = (int)mxGetScalar (prhs [5]);
        contain = (int)mxGetScalar (prhs [6]);

	/* require memory for return */
	plhs[0] = mxCreateDoubleMatrix(rows*cols*slices, 1, mxREAL);
	ridge = (double*) mxGetPr (plhs [0]);	
	
        for(i=0; i<rows; i++)
		for(j=0; j<cols; j++)
			for(k=0; k<slices; k++)
                        {
                                *(ridge+i+j*rows+k*cols*rows)=0;
                        }


        for(i=0; i<(rows-6); i++)
	         for(j=0; j<(cols-6); j++)
			for(k=0; k<(slices-6); k++)
			{
                            for(a=0;a<7;a++)
                                 for(b=0;b<7;b++)
                                     for(c=0;c<7;c++)
                                        *(temp+a+b*7+c*7*7)=*(DisMatrix+(i+a)+(b+j)*rows+(c+k)*cols*rows);
         
                                         
                            compare(temp, 343);
                            count(temp, 343);

                            if(Max>=th && Max>0 && total<=contain && total>0)
                            {
                                 if( *(temp+171)==Max )
                                     *(ridge+(i+1)+(j+1)*rows+(k+1)*cols*rows)=255;
                            }
                                 
       
                                 
			}                                	 

}