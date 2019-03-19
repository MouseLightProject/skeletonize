
#include <math.h>
#include <stdio.h>
#include "mex.h"

/*
double Max;

void compare(double *a, int m)
{
    int j;
    Max=a[0];
    index=0;

    for(j=1;j<=m;j++)
    {
       if(Max<a[j])
       { 
          Max=a[j];
          index=j;
       }
       
    }
}
*/

void mexFunction(int nlhs, mxArray *plhs[], int	nrhs, const mxArray *prhs[]) 
{
	int       rows, cols, slices;
        double    *voxel;
        double    *ridge_pos;
        double    *Dist;
	int       ridge_num;
        int       i,j,k;
        double    vx,vy,vz;
        double    d;
	double    start_x, start_y, start_z, end_x, end_y, end_z;


	/* reading the parameters */
        voxel = (double *) mxGetPr (prhs [0]);
        ridge_pos = (double *) mxGetPr (prhs [1]);
        ridge_num = (int)mxGetScalar (prhs [2]);
	rows = (int)mxGetScalar (prhs [3]);
	cols = (int)mxGetScalar (prhs [4]);
	slices = (int)mxGetScalar (prhs [5]);

	/* require memory for return */
	plhs[0] = mxCreateDoubleMatrix(ridge_num*ridge_num, 1, mxREAL);
	Dist = (double *) mxGetPr(plhs [0]);	

        
        for(i=0; i<ridge_num; i++)
            for(j=0; j<ridge_num; j++)
            {
               start_x=*(ridge_pos+i);
               start_y=*(ridge_pos+i+ridge_num);
               start_z=*(ridge_pos+i+2*ridge_num);
               end_x=*(ridge_pos+j);  
               end_y=*(ridge_pos+j+ridge_num);
               end_z=*(ridge_pos+j+2*ridge_num);
       
               d = sqrt((end_x-start_x)*(end_x-start_x) + (end_y-start_y)*(end_y-start_y) + (end_z-start_z)*(end_z-start_z));;
               *(Dist+i+j*ridge_num) = d;
               vx=(end_x-start_x)/d;
               vy=(end_y-start_y)/d;
               vz=(end_z-start_z)/d;
               

               while(abs(start_x-end_x)>=abs(vx))
               {
                    start_x=start_x+vx;
                    start_y=start_y+vy;
                    start_z=start_z+vz;
            
                    if( (int)(start_x+0.5)>=rows || (int)(start_y+0.5)>=cols || (int)(start_z+0.5)>=slices )
                    break;

                    if( (int)(start_x+0.5)<0 || (int)(start_y+0.5)<0 || (int)(start_z+0.5)<0 )
                    break;

                    if( *(voxel+(int)(start_x-0.5) + (int)(start_y-0.5)*rows + (int)(start_z-0.5)*cols*rows) == 0 )
                    {
                       *(Dist+i+j*ridge_num) = 100000;
                       break;
                    }
               }
            
             }
       
}
