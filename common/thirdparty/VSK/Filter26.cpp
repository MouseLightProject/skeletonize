
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
        //double    *temp=(double *)malloc(27*sizeof(double));
        double    *temp;
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


        for(i=0; i<(rows-2); i++)
	         for(j=0; j<(cols-2); j++)
			for(k=0; k<(slices-2); k++)
			{
                            temp=(double *)malloc(27*sizeof(double));
                            for(a=0;a<3;a++)
                                 for(b=0;b<3;b++)
                                     for(c=0;c<3;c++)
                                        {
                                           *(temp+a+b*3+c*3*3)=*(DisMatrix+(i+a)+(b+j)*rows+(c+k)*cols*rows);
                                        } 
         
                                         
                            compare(temp, 27);
                            count(temp, 27);

                            if(Max>=th && Max>0 && total<=contain && total>0)
                            {
                                 if( *(temp+13)==Max )
                                     *(ridge+(i+1)+(j+1)*rows+(k+1)*cols*rows)=255;
                            }
                                 
       
                                 
			}                                	 

}