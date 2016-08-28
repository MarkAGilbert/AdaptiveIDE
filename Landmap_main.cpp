#include<iostream>
#include <stdio.h>
#include <stdlib.h>
#include<fstream>
#include <string.h>
#include "fftw3.h"

#include "Simulate_adaptive.hpp"
#include "Parameters.hpp"
//#include "Add_Fine_Detail.hpp"

using namespace std;


int main( ){

    Parameters P;
	
	/*int locx = 1;
	int locy = 1;

	fstream file (P.map_file, ios::in|ios::binary);
    bool read;
	double mmax = 0;
	double mmin = 0;
	for(int s=0; s<P.stsq; s++){
		for(int n=0; n<P.region_width*P.landscape_length; n++){ // Get each row of map data individually
			int Na = n*P.region_width*P.landscape_length; // location of first element of data
			file.seekg(Na*(sizeof read));
			for(int m=0; m<P.region_width*P.landscape_length; m++){
				file.read(reinterpret_cast<char*>( &read ), sizeof read );
				if(double(read) > mmax) mmax = read;
				if(double(read) < mmin) mmin = read;
			}
			cout << mmin << "	" << mmax << endl;
		}
	}*/

	Simulate_adaptive(P);
	//if(P.smooth_output==1) Add_Fine_Detail(P);

}
