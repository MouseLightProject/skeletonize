
#include <math.h>
#include <stdio.h>
#include "mex.h"

double *CoorX;
double *CoorY;
double *CoorZ;

double compare(double *a, int m)
{
    int j;
    double Max=a[0];

    for(j=1;j<m;j++)
    {
       if(Max<a[j])
       { 
          Max=a[j];
       }       
    }
    
    return Max;
}



int extract(double *voxel,double x1,double y1,double z1,double x2,double y2,double z2,int up, int rows, int cols, int slices)
{
	int i, step, number;
	CoorX = (double *)malloc(up*sizeof(double));
	CoorY = (double *)malloc(up*sizeof(double));
	CoorZ = (double *)malloc(up*sizeof(double));

    double deltaX = abs(x1 - x2);
    double deltaY = abs(y1 - y2);
    double deltaZ = abs(z1 - z2);
    double *delta = (double *)malloc(3*sizeof(double));

    *(delta)=deltaX;
    *(delta+1)=deltaY;
    *(delta+2)=deltaZ;
   
    double del=compare(delta, 3); 
    
    double sumup = sqrt((x2-x1)*(x2-x1) + (y2-y1)*(y2-y1) + (z2-z1)*(z2-z1));
    double *vector = (double *)malloc(3*sizeof(double));

    *(vector)=(x2-x1)/sumup;
    *(vector+1)=(y2-y1)/sumup;
    *(vector+2)=(z2-z1)/sumup;

    double vec=compare(vector, 3);

    *(CoorX)=x1;
    *(CoorY)=y1;
    *(CoorZ)=z1;

    step=2;
    while(del>=vec)
    {
         x1=x1+*(vector);
         y1=y1+*(vector+1);
         z1=z1+*(vector+2);
                    
         if(step>=up-3)
         break;
    
		 if( (int)(x1+0.5)>=rows || (int)(y1+0.5)>=cols || (int)(z1+0.5)>=slices )
         break;

         if( (int)(x1+0.5)<=0 || (int)(y1+0.5)<=0 || (int)(z1+0.5)<=0 )
         break;
		 
		 if( *(voxel + (int)(x1+0.5) + ((int)(y1+0.5))*rows + ((int)(z1+0.5))*cols*rows )==0 )
         break;
		
		 
         *(CoorX+step-1)=x1;
         *(CoorY+step-1)=y1;
         *(CoorZ+step-1)=z1;
    
         deltaX = abs(x1 - x2);
         deltaY = abs(y1 - y2);
         deltaZ = abs(z1 - z2);
         
         *(delta)=deltaX;
         *(delta+1)=deltaY;
         *(delta+2)=deltaZ;
         del=compare(delta, 3);
		 step++;
     }

     *(CoorX+step-1)=x2;
     *(CoorY+step-1)=y2;
     *(CoorZ+step-1)=z2;
	 number=step;

	 return number;
}




void mexFunction(int nlhs, mxArray *plhs[], int	nrhs, const mxArray *prhs[]) 
{
	    int       rows, cols, slices;
		double    *voxel;
        double    *skeleton_X;
        double    *skeleton_Y;
        double    *skeleton_Z;
        double    *T;
        double    *ridge_pos;
	    int       connect_num;
        int       i,j,k;
        int       count;
	    double    start_x, start_y, start_z, end_x, end_y, end_z;
	    int num;
	    int M,N;


	    /* reading the parameters */
		voxel = (double *) mxGetPr (prhs [0]);
        T = (double *) mxGetPr (prhs [1]);
        ridge_pos = (double *) mxGetPr (prhs [2]);
        connect_num = (int)mxGetScalar (prhs [3]);
	    rows = (int)mxGetScalar (prhs [4]);
	    cols = (int)mxGetScalar (prhs [5]);
	    slices = (int)mxGetScalar (prhs [6]);
        count = (int)mxGetScalar (prhs [7]);

		M=mxGetM(prhs[1]); 
        N=mxGetN(prhs[1]); 


	   /* require memory for return */
	    plhs[0] = mxCreateDoubleMatrix(connect_num*count, 1, mxREAL);
	    skeleton_X = (double *) mxGetPr(plhs [0]);
        plhs[1] = mxCreateDoubleMatrix(connect_num*count, 1, mxREAL);
	    skeleton_Y = (double *) mxGetPr(plhs [1]);
        plhs[2] = mxCreateDoubleMatrix(connect_num*count, 1, mxREAL);
	    skeleton_Z = (double *) mxGetPr(plhs [2]);	

        int index=0;
        for(i=0; i<M; i++)
            for(j=i; j<N; j++)
            {
               if( *(T+i+j*M)==1 )
               {
                  start_x=*(ridge_pos+i);
                  start_y=*(ridge_pos+i+M);
                  start_z=*(ridge_pos+i+2*M);
                  end_x=*(ridge_pos+j);  
                  end_y=*(ridge_pos+j+M);
                  end_z=*(ridge_pos+j+2*M);
       
                  num=extract(voxel,start_x,start_y,start_z,end_x,end_y,end_z,count,rows, cols, slices);
               
                  for(k=0; k<num; k++)
                  {
					  *(skeleton_X+k+index*count) = *(CoorX+k);
					  *(skeleton_Y+k+index*count) = *(CoorY+k);
					  *(skeleton_Z+k+index*count) = *(CoorZ+k);
                  }
				  
                  for(k=num; k<count; k++)
                  {
                      *(skeleton_X+k+index*count) = *(CoorX+num-1);
                      *(skeleton_Y+k+index*count) = *(CoorY+num-1);
                      *(skeleton_Z+k+index*count) = *(CoorZ+num-1);
                  }				  

                  index++;
                }            
             }
       
}



















