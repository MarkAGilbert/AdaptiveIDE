#ifndef KERNELHEADERDEF
#define KERNELHEADERDEF

#include<vector>
#include<iostream>
#include "Parameters.hpp"
#include<cmath>
#define PI 3.14159265358979323846

using namespace std;

// This file generates the HR-local, LR-local and long-distance dispersal kernels. Local_Kernel generates the long distance kernel and Global_Kernel generates the long-distance dispersal kernel.
// The type of kernel is specified by P.kernel_number (1 - tophat, 2 - cone, 3 - hemisphere, 4 - gaussian, 5 - laplace, 6 - test)

void Local_Kernel(Parameters P, vector<fftw_complex*> kernel, vector<double*> lckernel)
{
	fftw_complex* ker;
	fftw_complex* ker2;
	fftw_plan p_fwd;
	double* pyramid;
	
	for(int s=0; s<P.stsq; s++){	

		ker = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*pow(4.0*P.region_width,2));
		ker2 = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*pow(4.0*P.region_width,2));
		p_fwd = fftw_plan_dft_1d(pow(4.0*P.region_width,2),ker2,kernel[s],1,FFTW_ESTIMATE);

		// Assign Values to the Kernel
		double deltasquared = pow(1.0*P.length/(P.region_width*P.landscape_width),2);
		for(int n1=0; n1<4*P.region_width; n1++){
			for(int n2=0; n2<4*P.region_width; n2++){
				ker[n1+4*P.region_width*n2][1] = 0;
				// Assigning Values to real part of each element
				double distancefromorigin = sqrt(pow(1.0*(P.length/(P.landscape_width*P.region_width))*(n1-2*P.region_width),2) + pow(1.0*(P.length/(P.landscape_width*P.region_width))*(n2-2*P.region_width),2));
				// Tophat
				if(P.kernel_no == 1){
					if(distancefromorigin < P.dispersal_distance[s]){
						ker[n1+4*P.region_width*n2][0] = deltasquared/(PI*pow(1.0*P.dispersal_distance[s],2));
					}
					else ker[n1+4*P.region_width*n2][0] = 0;
				}
				// Cone
				if(P.kernel_no == 2){
					if(distancefromorigin < P.dispersal_distance[s]){
						ker[n1+4*P.region_width*n2][0] = 3*deltasquared*(1-distancefromorigin/P.dispersal_distance[s])/(PI*pow(1.0*P.dispersal_distance[s],2));
					}
					else ker[n1+4*P.region_width*n2][0] = 0;
				}
				// Hemisphere
				if(P.kernel_no == 3){
					if(distancefromorigin < P.dispersal_distance[s]){
						ker[n1+4*P.region_width*n2][0] = deltasquared*sqrt(pow(1.0*P.dispersal_distance[s],2)-pow(1.0*distancefromorigin,2))/(2*PI*pow(1.0*P.dispersal_distance[s],3)/3);
					}
					else ker[n1+4*P.region_width*n2][0] = 0;
				}
				// Gaussian
				if(P.kernel_no == 4){
					ker[n1+4*P.region_width*n2][0] = deltasquared*exp(-pow(1.0*distancefromorigin,2)/(2*pow(1.0*P.dispersal_distance[s],2)))/(pow(1.0*P.dispersal_distance[s],2)*2*PI);
				}
				// Laplace
				if(P.kernel_no == 5){
					ker[n1+4*P.region_width*n2][0] = deltasquared*exp(-distancefromorigin/P.dispersal_distance[s])/(2*PI*pow(1.0*P.dispersal_distance[s],2));
				}
				// Test
				if(P.kernel_no == 6){
					ker[n1+4*P.region_width*n2][0] = 0;
					if(n1==2*P.region_width && n2==2*P.region_width+1) ker[n1+4*P.region_width*n2][0]=.5; //.5;
				}
			}
		}
	/*    double summa = 0;
		for(int n1=0; n1<4*P.region_width; n1++){
			for(int n2=0; n2<4*P.region_width; n2++){
				summa += ker[n1+4*P.region_width*n2][0];
			}
		}*/


		// Build the local coarse kernel
		int Pa = pow(2.0*P.region_width,2);
		pyramid = new double[Pa];
		for(int n=0; n<2*P.region_width; n++) for(int m=0; m<2*P.region_width; m++) pyramid[n+2*P.region_width*m] = 0;
		for(int n=1; n<2*P.region_width; n++) for(int m=0; m<2*P.region_width; m++) pyramid[n+2*P.region_width*m] = (P.region_width-abs(n-P.region_width))*(P.region_width-abs(m-P.region_width));
		for(int n=0; n<3; n++){
			for(int m=0; m<3; m++){
				lckernel[s][n+3*m]=0;
				for(int j1=0; j1<2*P.region_width; j1++) for(int j2=0; j2<2*P.region_width; j2++) lckernel[s][n+3*m]+=pyramid[j1+2*P.region_width*j2]*ker[P.region_width*n+j1+4*P.region_width*(P.region_width*m+j2)][0]/pow(1.0*P.region_width,2);
			}
		}

		// Reshape the kernel to have the origin value at [0][0]
		for(int n=0; n<(8*pow(1.0*P.region_width,2)-2*P.region_width); n++){
			int N1 = 8*pow(1.0*P.region_width,2)+2*P.region_width;
			ker2[n][0] = ker[n+N1][0], ker2[n][1] = ker[n+N1][1];
		}
		for(int n=0; n<8*pow(1.0*P.region_width,2)+2*P.region_width; n++){
			int N1 = 8*pow(1.0*P.region_width,2)-2*P.region_width;
			ker2[n+N1][0] = ker[n][0], ker2[n+N1][1] = ker[n][1];
		}
		// FFT it to give kernel
	   fftw_execute(p_fwd);
	   	fftw_free(ker);
		fftw_free(ker2);
		fftw_destroy_plan(p_fwd);
		delete[](pyramid);
	}
}

void Global_Kernel(Parameters P, vector<fftw_complex*> gkernel){

    fftw_complex* gker;
    fftw_complex* gker2;
    fftw_plan p_gfwd;
	for(int s=0; s<P.stsq; s++){
		gker = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*4*P.landscape_length*P.landscape_width);
		gker2 = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*4*P.landscape_length*P.landscape_width);
		p_gfwd = fftw_plan_dft_1d(4*P.landscape_length*P.landscape_width,gker2,gkernel[s],1,FFTW_ESTIMATE);

		// Assign Values to the Kernel
		if(((P.kernel_no == 1) || (P.kernel_no == 2) || (P.kernel_no == 3)) && (P.dispersal_distance[s] > P.length/P.landscape_width)) cout <<"WARNING: Global dispersal for kernels with finite support not yet supported";
		for(int n1=0; n1<2*P.landscape_length; n1++){
			for(int n2=0; n2<2*P.landscape_width; n2++){
				gker[n1+2*P.landscape_length*n2][1] = 0;
				double distancefromorigin =  (P.length/P.landscape_width)*sqrt(pow(1.0*(n1-P.landscape_length),2)+pow(1.0*(n2-P.landscape_width),2));
				double Deltasquared = pow(1.0*P.length/P.landscape_width,2);
				// Compact Support (see above warning)
				if((P.kernel_no == 1) || (P.kernel_no == 2) || (P.kernel_no == 3)){
					gker[n1+2*P.landscape_length*n2][0] = 0;
				}
				// Gaussian
				if(P.kernel_no == 4) gker[n1+2*P.landscape_length*n2][0] = Deltasquared*exp(-pow(1.0*distancefromorigin,2)/(2*pow(1.0*P.dispersal_distance[s],2)))/(pow(1.0*P.dispersal_distance[s],2)*2*PI);

				// Laplace
				if(P.kernel_no == 5) gker[n1+2*P.landscape_length*n2][0] = Deltasquared*exp(-fabs(distancefromorigin/P.dispersal_distance[s]))/(2*PI*pow(1.0*P.dispersal_distance[s],2));

				// Test Case
				if(P.kernel_no == 6){
					gker[n1+2*P.landscape_length*n2][0] = 0;
					if((abs(n1-P.landscape_length) < 1.2)&&(abs(n2-P.landscape_width) <1.5)) gker[n1+2*P.landscape_length*n2][0] = .05;//.5;
				}

				// Get rid of middle elements
				if((abs(n1-P.landscape_length)<1.5)&&(abs(n2-P.landscape_width)<1.5)) gker[n1+2*P.landscape_length*n2][0] = 0;
			}
		}

		// Reshape the kernel to have the origin value at [0][0]
		int N1 = 2*P.landscape_length*P.landscape_width+P.landscape_length;
		int N2 = 2*P.landscape_length*P.landscape_width-P.landscape_length;
		for(int n=0; n<N2; n++){
			gker2[n][0] = gker[n+N1][0], gker2[n][1] = gker[n+N1][1];
		}
		for(int n=0; n<N1; n++){
			gker2[n+N2][0] = gker[n][0], gker2[n+N2][1] = gker[n][1];
		}

	   fftw_execute(p_gfwd);

	   fftw_free(gker);
	   fftw_free(gker2);
	   fftw_destroy_plan(p_gfwd);
	}
}

#endif

