
#include <math.h>
#include <stdio.h>
#include "mex.h"
       

void mexFunction(int nlhs, mxArray *plhs[], int	nrhs, const mxArray *prhs[]) 
{
	int       rows, cols, slices;
        double       *voxel;
        double       *boundary;
        int       i,j,k;

	/* reading the parameters */
        voxel = (double *) mxGetPr (prhs [0]);
	rows = (int)mxGetScalar (prhs [1]);
	cols = (int)mxGetScalar (prhs [2]);
	slices = (int)mxGetScalar (prhs [3]);

	/* require memory for return */
	plhs[0] = mxCreateDoubleMatrix(rows*cols*slices, 1, mxREAL);
	boundary = (double *) mxGetPr (plhs [0]);	
	
        for(i=0; i<rows; i++)
		for(j=0; j<cols; j++)
			for(k=0; k<slices; k++)
                        {
                                *(boundary+i+j*rows+k*cols*rows)=*(voxel+i+j*rows+k*cols*rows);
                                //*(boundary+i+j*rows+k*cols*rows)=0;
                        }


    for(i=1; i<(rows-1); i++)
		for(j=1; j<(cols-1); j++)
			for(k=1; k<(slices-1); k++)
			{

                                if( *(voxel+(i-1)+j*rows+k*cols*rows)==0 && *(voxel+i+j*rows+k*cols*rows)==1 )
                                {
                                        *(boundary+i+j*rows+k*cols*rows)=255;
                                }
                                
                                if( *(voxel+(i+1)+j*rows+k*cols*rows)==0 && *(voxel+i+j*rows+k*cols*rows)==1 )
                                {
                                        *(boundary+i+j*rows+k*cols*rows)=255;
                                }

                                if( *(voxel+i+(j-1)*rows+k*cols*rows)==0 && *(voxel+i+j*rows+k*cols*rows)==1 )
                                {
                                        *(boundary+i+j*rows+k*cols*rows)=255;
                                }
                              
                                if( *(voxel+i+(j+1)*rows+k*cols*rows)==0 && *(voxel+i+j*rows+k*cols*rows)==1 )
                                {
                                        *(boundary+i+j*rows+k*cols*rows)=255;
                                }

                                if( *(voxel+i+j*rows+(k-1)*cols*rows)==0 && *(voxel+i+j*rows+k*cols*rows)==1 )
                                {
                                        *(boundary+i+j*rows+k*cols*rows)=255;
                                }

                                if( *(voxel+i+j*rows+(k+1)*cols*rows)==0 && *(voxel+i+j*rows+k*cols*rows)==1 )
                                {
                                        *(boundary+i+j*rows+k*cols*rows)=255;
                                }

			}

		
		for(j=0; j<cols; j++)
			for(k=0; k<slices; k++)
			{
				if( *(voxel+0+j*rows+k*cols*rows)==1 )
					*(boundary+0+j*rows+k*cols*rows)=255;

				if( *(voxel+(rows-1)+j*rows+k*cols*rows)==1 )
                    *(boundary+(rows-1)+j*rows+k*cols*rows)=255;
			}

	    for(i=0; i<rows; i++)
			for(k=0; k<slices; k++)
			{
				if( *(voxel+i+0*rows+k*cols*rows)==1 )
					*(boundary+i+0*rows+k*cols*rows)=255;

				if( *(voxel+i+(cols-1)*rows+k*cols*rows)==1 )
                    *(boundary+i+(cols-1)*rows+k*cols*rows)=255;
			}

       for(i=0; i<rows; i++)
			for(j=0; j<cols; j++)
			{
				if( *(voxel+i+j*rows+0*cols*rows)==1 )
					*(boundary+i+j*rows+0*cols*rows)=255;

				if( *(voxel+i+j*rows+(slices-1)*cols*rows)==1 )
                    *(boundary+i+j*rows+(slices-1)*cols*rows)=255;
			}

}

