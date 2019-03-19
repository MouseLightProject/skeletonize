#include <math.h>
#include <stdio.h>
#include "mex.h"

       

void mexFunction(int nlhs, mxArray *plhs[], int	nrhs, const mxArray *prhs[]) 
{
	int       rows, cols, slices;
        double    *DisMatrix;
        double    *boundary_pos;
        double    *object_pos;
	int       boundary_num;
	int       object_num;
        int       i,j,k;
        int       x1,y1,z1;
        int       x2,y2,z2;
	double    Min;
        double    temp;

	/* reading the parameters */
        boundary_pos = (double *) mxGetPr (prhs [0]);
        object_pos = (double *) mxGetPr (prhs [1]);
        object_num = (int)mxGetScalar (prhs [2]);
        boundary_num = (int)mxGetScalar (prhs [3]);
	rows = (int)mxGetScalar (prhs [4]);
	cols = (int)mxGetScalar (prhs [5]);
	slices = (int)mxGetScalar (prhs [6]);

	/* require memory for return */
	plhs[0] = mxCreateDoubleMatrix(rows*cols*slices, 1, mxREAL);
	DisMatrix = (double *) mxGetPr(plhs [0]);	
	
      for(i=0; i<rows; i++)
	    for(j=0; j<cols; j++)
		 for(k=0; k<slices; k++)
                    {
                          *(DisMatrix+i+j*rows+k*cols*rows)=0;
                    }


        for(i=0; i<object_num; i++)
        {
            x1=*(object_pos+i);
            y1=*(object_pos+i+object_num);
            z1=*(object_pos+i+2*object_num);          
            
            Min=1000;
	    for(j=0; j<boundary_num; j++)
	    {
         	x2=*(boundary_pos+j);
                y2=*(boundary_pos+j+boundary_num);
                z2=*(boundary_pos+j+2*boundary_num);
                //*(distance_temp+j) = (x2-x1)*(x2-x1) + (y2-y1)*(y2-y1) + (z2-z1)*(z2-z1);
                temp = (x2-x1)*(x2-x1) + (y2-y1)*(y2-y1) + (z2-z1)*(z2-z1);
                if(temp<Min)
                   Min=temp;
            }
            
            *(DisMatrix+(x1-1)+(y1-1)*rows+(z1-1)*cols*rows)=Min;	

        }



}