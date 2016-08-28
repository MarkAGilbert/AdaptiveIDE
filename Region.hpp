#ifndef REGION1HEADERDEF
#define REGION1HEADERDEF

#include<vector>
#include<iostream>
#include<fstream>
#include "fftw3.h"
#include<math.h>
#include<string.h>
#include "Parameters.hpp"

using namespace std;

class Region
{

public:
  // Constructor
  Region(int no_stages); // Constructor
  Region(const Region& Other); // Copy Constructor

  // Destructor
  ~Region();

  // Elements
  int stages;
  int sqst;
  int width; // Width (in cells) of the square region
  int extended_elements; // Number of cells in the square extended region required for dispersal between the region and its 8 nearest neighbours (equal to the extended_width squared)
  int extended_width; // Width (in cells) of the square extended region required for dispersal between the region and its 8 nearest neighbours
  int elements; // Number of cells in the square region (equal to the width squared)
  bool fine; // Level of discretisation. If fine = 0, the population is represented by a single element. If fine = 1, the population is represented by many elements
  bool used;

  vector<double> final_popn; // The population (which stays fixed after the Region has become 'used')
  vector< double > growth_rate; // The average intrinsic (low population) growth rate for the region
  //vector< double> threshold_weights; // The threshold weights for the region
  vector< double > OldAve; // The average population in the region at the previous time step (used and updated by Landscape::Remesh)

// Functions
  void Grow(); // Growth phase for the region
  void GrowthRate(double* intrinsic_growth_rate, const char * map_file, const char * map_type, int locx, int locy, int fine_width, int map_length); // Set growth_rate by calculating proportion of region that is suitable habitat and multiplying by intrinsic growth rate for suitable habitat and (if fine = 1) generate the local growth array
  void Sum(); // Sums up stage-pairs after dispersal
  void Disperse(vector<fftw_complex*> kernel, Parameters P); // Local dispersal phase for regions with a fine mesh (fine = 1)
  void Disperse(vector<double*> kernel, Parameters P); // Local dispersal phase for regions with a coarse mesh (fine = 0)
  void Share(int n, int m, Region* pNeighbour); // (n,m)th section of extended distribution is added to *pNeighbour's distribution
  void Add(double x, int stage_no); // Adds constant x to each element of the (population) distribution
  double Average(int stage_no); // Returns the average of the (population) distribution
  double Minimum(); // Returns the minimum of the (population) distribution
  void Fine(Parameters P, int N, const char* map_file, const char* map_type, int locx, int locy, int map_length, double* intrinsic_growth_rate); // Changes the region's mesh to fine (from coarse)
  void Coarse(int locx, int locy, Parameters P); // Changes the region's mesh to coarse (from fine)
  
//private:
  //double average;
  vector<fftw_complex*> distribution; // Population distribution on the extended domain
  vector<fftw_complex*> fourier_distribution; // Fourier transform of population distribution
  vector<fftw_plan> p_forward;
  vector<fftw_plan> p_backward; // Paths for fourier transforms between distribution and fourier_distribution using fftw3
  vector<double**> growth; // Array of instrinsic growth rates for the each element in the Region
};

// Constructor
Region::Region(int no_stages)
{	
	stages = no_stages, sqst = stages*stages, fine = 0, width=1, elements=1, extended_elements=9, extended_width=3, used = 0; // Start with coarse mesh (fine=0), a single elements and 9 extended elements (to allow dispersal to neighbours)
	growth_rate.resize(sqst);
	OldAve.resize(sqst);
	for(int s=0; s<sqst; s++) growth_rate[s] = 1, OldAve[s] = 0;
	distribution.resize(sqst);
	fourier_distribution.resize(sqst);
	growth.resize(sqst);
	for(int m=0; m<sqst; m++){
		distribution[m] = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*9); // Allocate the distribution ove the extended domain
		for(int n=0; n<extended_elements; n++) distribution[m][n][0]=0, distribution[m][n][1]=0; // Initialise the distriubtion to zero
		growth[m] = new double*[1]; // Allocate the growth array for the (non-extended) domain
		growth[m][0] = new double[1];
		growth[m][0][0] = growth_rate[m]; // Assign the growth array the average intrinsic growth rate
	}
	p_forward.resize(sqst);
	p_backward.resize(sqst);
}

// Copy Constructor
Region::Region(const Region& Other)
{
    stages = Other.stages, sqst = stages*stages, fine = 0, width=1, elements=1, extended_elements=9, extended_width=3, used = 0; // Start with coarse mesh (fine=0), a single elements and 9 extended elements (to allow dispersal to neighbours)
    growth_rate.resize(sqst);
	OldAve.resize(sqst);
	for(int s=0; s<sqst; s++) growth_rate[s] = 1, OldAve[s] = 0;
	distribution.resize(sqst);
	fourier_distribution.resize(sqst);
	growth.resize(sqst);
	for(int m=0; m<sqst; m++){
		distribution[m] = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*9); // Allocate the distribution ove the extended domain
		for(int n=0; n<extended_elements; n++) distribution[m][n][0]=0, distribution[m][n][1]=0; // Initialise the distriubtion to zero
		growth[m] = new double*[1]; // Allocate the growth array for the (non-extended) domain
		growth[m][0] = new double[1];
		growth[m][0][0] = growth_rate[m]; // Assign the growth array the average intrinsic growth rate
	}
	p_forward.resize(sqst);
	p_backward.resize(sqst);
}

// Destructor
Region::~Region()
{
	for(int m=0; m<sqst; m++){	
		if(fine==0){ // If the mesh is coarse, we have only assigned the distribution and growth array, if it is fine we also need to free memory from the fourier_distribution and the fftw plans.
			delete[](growth[m][0]); // Free memory from each row of the growth array
		}
		else{
			fftw_free(distribution[m]); // Free memory from distribution
			fftw_free(fourier_distribution[m]); // Free memory from the fourier_distribution
			for(int n=0; n<width; n++){
				delete[](growth[m][n]); // Free memory from each row of the growth array
			}
		}
		delete[](growth[m]); // Free memory from the growth array
	}
}

// Find the average intrinsic growth rate by reading map data from an external file (map_file) and generate the local growth array (if fine = 1)
void Region::GrowthRate(double* intrinsic_growth_rate, const char * map_file, const char * map_type, int locx, int locy, int fine_width, int map_length)
{
    typedef bool Typ;

    if(strcmp(map_type,"bool")==0) typedef bool Typ;
    else if(strcmp(map_type,"double")==0) typedef double Typ;
    else cout << "Error!";
    fstream file (map_file, ios::in|ios::binary);
    Typ read;
    double proportion = 0; // proportion is the proportion of suitable habitat in the region (initially zero)
    for(int s=0; s<sqst; s++){
		for(int n=0; n<fine_width; n++){ // Get each row of map data individually
			int Na = fine_width*locx+map_length*pow(1.0*fine_width,2)*locy+n*fine_width*map_length; // location of first element of data
			file.seekg(Na*(sizeof read));
			for(int m=0; m<fine_width; m++){
				file.read(reinterpret_cast<char*>( &read ), sizeof read );
				proportion += int(read)/pow(double(fine_width),2); // For each element, if it is suitable habitat, then increment the proportion of good habitat by the relative size of that element
				if(fine == 1) growth[s][n][m] = intrinsic_growth_rate[s]*read;
			}
		}
		growth_rate[s] = intrinsic_growth_rate[s]*proportion; // Calculate and allocate the average intrinsic growth rate
		if(fine == 0) growth[s][0][0]=growth_rate[s]; // If the mesh is coares then assign the single growth array element the value of the average intrinsic growth rate
	}
    file.close();
}

// Growth phase for the region
void Region::Grow(){
    // High Resolution
	for(int s=0; s<sqst; s++) OldAve[s] = Average(s);
    if(fine == 1){
		for(int s=0; s<sqst; s++){
			for(int n=0; n<extended_width; n++){
				for(int m=0; m<extended_width; m++){
					if(n<2*width && n>width-1 && m <2*width && m>width-1);
					else distribution[s][n+extended_width*m][0]=0; // If the cell is not an element of the region, but a member of the extended region, then set its value to zero
					distribution[s][n+extended_width*m][1] = 0; // Make imaginary part of the population distribution zero
					if(distribution[s][n+extended_width*m][0]<1e-13) distribution[s][n+extended_width*m][0] = 0; // Remove numerical noise resulting from fftw in dispersal phase
				}
			}
		}
		for(int n=0; n<width; n++){
			for(int m=0; m<width; m++){
				double y = 0;
				for(int s2=0; s2<stages; s2++) for(int s1=0; s1<stages; s1++) y += growth[s1*stages+s2][n][m]*distribution[s2*(stages+1)][width+n+extended_width*(width+m)][0];
				if(y<1){
					for(int s2=0; s2<stages; s2++){
						double x = distribution[s2*(stages+1)][width+n+extended_width*(width+m)][0];
						for(int s1=0; s1<stages; s1++){
							distribution[s1*stages+s2][width+n+extended_width*(width+m)][0] = growth[s1*stages+s2][n][m] * x; // Population growth: if(r*u < 1) u' = r*u
						}
					}
				}
				else{
					for(int s2=0; s2<stages; s2++){
						double x = distribution[s2*(stages+1)][width+n+extended_width*(width+m)][0];
						for(int s1=0; s1<stages; s1++){
							distribution[s1*stages+s2][width+n+extended_width*(width+m)][0] = (1./y) * growth[s1*stages+s2][n][m] * x; // Population growth: if(r*u < 1) u' = r*u
						}
					}
				}
			}
		}
    }
    // Low Resolution, Wave Front
    else{
		for(int n=0; n<extended_width; n++){
			for(int m=0; m<extended_width; m++){
				if(n==1 && m==1){
					for(int s=0; s<sqst; s++) if(distribution[s][n+extended_width*m][0]<1e-13) distribution[s][n+extended_width*m][0] = 0; // Remove numerical noise resulting from fftw in dispersal phase
					double y = 0;
					for(int s2=0; s2<stages; s2++) for(int s1=0; s1<stages; s1++) y += growth[s1*stages+s2][0][0]*distribution[s2*(stages+1)][n+extended_width*m][0];
					if(y<1){
						for(int s2=0; s2<stages; s2++){
							double x = distribution[s2*(stages+1)][n+extended_width*m][0];
							for(int s1=0; s1<stages; s1++){
								distribution[s1*stages+s2][n+extended_width*m][0] = growth[s1*stages+s2][0][0] * x;
							}
						}
					}
					else{
						for(int s2=0; s2<stages; s2++){
							double x = distribution[s2*(stages+1)][n+extended_width*m][0];
							for(int s1=0; s1<stages; s1++){
								distribution[s1*stages+s2][n+extended_width*m][0] = (1./y) * growth[s1*stages+s2][0][0] * x;
							}
						}
					}
				}
				else{
					for(int s=0; s<sqst; s++) distribution[s][n+extended_width*m][0] = 0; // If the cell is not an element of the region, but a member of the extended region, then set its value to zero
				}
				for(int s=0; s<sqst; s++) distribution[s][n+extended_width*m][1] = 0; // Make imaginary part of the population distribution zero
			}
		}
	}
}
					

/*				for(int s2=0; s2<stages; s2++){
					double x = distribution[s2*(stages+1)][n+extended_width*m][0];
					for(int s1=0; s1<stages; s1++){
						if(n==1 && m==1){
							if(distribution[s1*stages+s2][n+extended_width*m][0]<1e-13) distribution[s1*stages+s2][n+extended_width*m][0] = 0; // Remove numerical noise resulting from fftw in dispersal phase
							if(distribution[0][n+extended_width*m][0]<1) distribution[s1*stages+s2][n+extended_width*m][0] = growth[s1*stages+s2][0][0] * x; // Population growth: if(r*u < 1) u' = r*u 
						}
						else distribution[s1*stages+s2][n+extended_width*m][0] = 0; // If the cell is not an element of the region, but a member of the extended region, then set its value to zero
						distribution[s1*stages+s2][n+extended_width*m][1] = 0; // Make imaginary part of the population distribution zero
					}
				}
			}
		}
    }
    // Low Resolution, Wave Back
    else{
		for(int s=0; s<sqst; s++){
			  for(int n=0; n<extended_width; n++){
				for(int m=0; m<extended_width; m++){
				  if(n==1 && m==1){
					  //distribution[n+extended_width*m][0] = final_popn; // Population remains constant in the wave-back
				  }
				  else distribution[s][n+extended_width*m][0] = 0; // If the cell is not an element of the region, but a member of the extended region, then set its value to zero
				  distribution[s][n+extended_width*m][1] = 0; // Make imaginary part of the population distribution zero
				}
			  }
		}
    }
}*/

// Local dispersal phase for regions with a fine mesh (fine = 1)
void Region::Disperse(vector<fftw_complex*> kernel, Parameters P){
  // Fine
  if(fine != 1); // Make sure that region has a fine mesh (otherwise argument should be double*)
  else{
	for(int s=0; s<sqst; s++){
		if(P.dispersal_on[s]==1){
			p_forward[s] = fftw_plan_dft_1d(extended_elements,distribution[s],fourier_distribution[s],1,FFTW_ESTIMATE); // Create plan for fft
			fftw_execute(p_forward[s]); // Execute fft
			// Convolution - elementwise multiplication of the fourier distribution with the fft of the dispersal kernel
			for(int n=0; n<extended_elements; n++){
				 double a = fourier_distribution[s][n][0];
				 double b = fourier_distribution[s][n][1];
				 double c = kernel[s][n][0];
				 double d = kernel[s][n][1];
				 fourier_distribution[s][n][0] =  ((a*c)-(b*d))/double(extended_elements); // Divide by the size of the arrays to keep population constant
				 fourier_distribution[s][n][1] = ((a*d)+(b*c))/double(extended_elements); // Need to do imaginary and real parts
			}
			p_backward[s] = fftw_plan_dft_1d(extended_elements,fourier_distribution[s],distribution[s],-1,FFTW_ESTIMATE); // Create plan for inverse fft
			fftw_execute(p_backward[s]); // Execute inverse fft
			fftw_destroy_plan(p_backward[s]);
			fftw_destroy_plan(p_forward[s]);
		}
	}
  }
}

// Local dispersal phase for regions with a coarse mesh (fine = 0)
void Region::Disperse(vector<double*> kernel, Parameters P){
  // Coarse
  if(fine != 0); // Make sure that region has a coarse mesh (otherwise argument should be fftw_complex*)
  else{
	  for(int s=0; s<sqst; s++){
		  if(P.dispersal_on[s]==1){
			double x = 0; // x is the total density of propagules arriving at the region (from its nearest neighbours)
			for(int n=0; n<extended_elements; n++) x += kernel[s][n]*distribution[s][n][0]; // Increment x for the contribution of each of itself and its 8 nearest neighbours
			distribution[s][extended_width+1][1] = 0; // Set the imaginary part of the population to 0
			distribution[s][extended_width+1][0] = x; // Set the population to x (the total density of propagules arriving from the nearest neighbours)
		  }
		}
  }
}

void Region::Sum(){
    // High Resolution
    if(fine == 1){
		for(int n=0; n<width; n++){
			for(int m=0; m<width; m++){
				for(int s1=0; s1<stages; s1++){
					double x = 0;
					for(int s2=0; s2<stages; s2++){ 
						x += distribution[s1*stages+s2][width+n+extended_width*(width+m)][0];
						distribution[s1*stages+s2][width+n+extended_width*(width+m)][0] = 0;
					}
					distribution[s1*(stages+1)][width+n+extended_width*(width+m)][0] = x;
				}
			}
		}
	}
	else{
		for(int s1=0; s1<stages; s1++){
			double x = 0;
			for(int s2=0; s2<stages; s2++){
				x += distribution[s1*stages+s2][4][0];
				distribution[s1*stages+s2][4][0] = 0;
			}
			distribution[s1*(stages+1)][4][0] = x;
		}
	}
}

// (n,m)th section of extended distribution is added to *pNeighbour's distribution
void Region::Share(int n, int m, Region* pNeighbour){
    // Fine Destination
    if((*pNeighbour).fine == 1){
		for(int s=0; s<sqst; s++){
			for(int i1=0; i1<(*pNeighbour).width; i1++){
				for(int i2=0; i2<(*pNeighbour).width; i2++){
					// Fine Origin
					if(fine == 1){
						(*pNeighbour).distribution[s][(1+n)*(*pNeighbour).width+i1+(*pNeighbour).extended_width*((1+m)*(*pNeighbour).width+i2)][0] = distribution[s][width+i1+extended_width*(width+i2)][0]; // Add population from corresponding cell in *pNeighbour to the element of distribution
					}
					// Coarse Origin
					else{
						(*pNeighbour).distribution[s][(1+n)*(*pNeighbour).width+i1+(*pNeighbour).extended_width*((1+m)*(*pNeighbour).width+i2)][0] = distribution[s][1+extended_width][0]; // Add non spatially structured average population of *pNeighbour to the element of distribution
					}
				}
			}
		}
    }
    // Coarse Destination
    else{
		for(int s=0; s<sqst; s++){
			// Fine Origin
			if(fine == 1){
				(*pNeighbour).distribution[s][1+n+(*pNeighbour).extended_width*(1+m)][0] = Average(s); // Add average population of fine neighbour to coarse (non spatially structured) distribution
			}
			// Coarse Origin
			else{
				(*pNeighbour).distribution[s][1+n+(*pNeighbour).extended_width*(1+m)][0] = distribution[s][1+extended_width][0]; // Add non spatially structured neighbour's population to non spatiall structured distribution
			}
		}
    }
}

// Adds constant x to each element of the (population) distribution
void Region::Add(double x, int stage_no){
  for(int n=0; n<width; n++){
        for(int m=0; m<width; m++){
                distribution[stage_no][width+n+extended_width*(width+m)][0]+=x; // Add x to each element of the distribution
        }
  }
}

// Returns the average of the (population) distribution
double Region::Average(int stage_no){
  double aa=0; // Initialise the average aa as zero
  for(int n=0; n<width; n++){
        for(int m=0; m<width; m++){
                aa += distribution[stage_no][width+n+extended_width*(width+m)][0]/double(elements); // Increments aa by the contribution to the average from each element in the (non-extended) domain
        }
  }
  return aa;
}

// Returns the minimum of the (population) distribution
double Region::Minimum(){
	double aa=0;
	for(int s=0; s<sqst; s++) aa += distribution[s][width+extended_width*width][0]; // Initialise the minimum aa as the first element in the (non-extended) domain
    for(int n=0; n<width; n++){
        for(int m=0; m<width; m++){
            double bb = 0; 
			for(int s=0; s<sqst; s++) bb += distribution[s][width+n+extended_width*(width+m)][0]; // bb is the (n,m)th element in the (non-extended) domain
            if(bb<aa) aa = bb; // If bb is smaller than aa then change aa to bb
        }
    }
    return aa;
}

// Changes the region's mesh to fine (from coarse)
void Region::Fine(Parameters P, int N, const char* map_file, const char* map_type, int locx, int locy, int map_length, double* intrinsic_growth_rate){
    if(fine==0){
        typedef bool Typ; // Programme reads map_file as Typ. This is set to double by default
        if(map_type == "bool") typedef bool Typ; // Set Typ to bool if map_type is bool
        else if(map_type == "double") typedef double Typ; // Set Typ to double if map_type is bool
        used = 1; // Change used to 1 (this means the region has entered the travelling wave)
		int coarse_extended_width = extended_width;
        fine = 1, width=N, elements=N*N, extended_elements=16*N*N, extended_width=4*N; // Update width, elements etc
		double x;
		for(int s=0; s<sqst; s++){ 
			x = distribution[s][coarse_extended_width + 1][0]; // the population x is the middle element in the 3x3 population distribution for the cell and its neighbours
			fftw_free(distribution[s]); // Free LR distribution, so a HR one can be allocated
			distribution[s] = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*extended_elements); // Allocate a HR distribution
			fourier_distribution[s] = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*extended_elements); // Allocate a HR fourier-distribution
			for(int n=0; n<extended_width; n++) for(int m=0; m<extended_width; m++) distribution[s][n+extended_width*m][0]=x, distribution[s][n+extended_width*m][1]=0; // All elements in the HR distribution are assumed to initially have the same density
			delete[](growth[s][0]); // Delete 1x1 growth array
			delete[](growth[s]); // Delete 1x1 growth array
			growth[s] = new double*[width]; // Allocate new HR growth array
			for(int n=0; n<width; n++) growth[s][n] = new double[width]; // Allocate new HR growth array
			fstream file (map_file, ios::in|ios::binary); // Create fstream for map_file
			Typ read; // Create variable of type Typ (see first few lines of this function Region::Fine)
			for(int n=0; n<width; n++){ // Get each row of map data individually
				int Na = width*locx+map_length*width*width*locy+n*width*map_length; // location of first element of data
				file.seekg(Na*(sizeof (Typ))); // Find first element in map_file
				for(int m=0; m<width; m++){ // For each element in the row of map data
					file.read(reinterpret_cast<char*>( &read ), sizeof read ); // Read each element of the map file
					growth[s][m][n] = intrinsic_growth_rate[s]*read; // (m,n)-th element of the region's array of growth rates is given by the intrinsic_growth_rate multiplied by the element of map_file
					//cout << read << "	";
				}
			}
			file.close(); // close the map_file
		}
    }
    else cout << "Error trying to fine a fine!"; // If the region is already fine, output error
}

// Changes the region's mesh to coarse (from fine)
void Region::Coarse(int locx, int locy, Parameters P){
	// Save to file
	ofstream abcdef;
	const char* file = "temporary_file.bin";
	abcdef.open(file, ios::out | ios::app | ios::binary);
	abcdef << locx << "	" << locy << "	"; // << endl;
	for(int m1=0; m1<P.region_width; m1+=P.print_size){
		for(int m2=0; m2<P.region_width; m2+=P.print_size){
			abcdef << distribution[0][P.region_width+m1+4*P.region_width*(m2+P.region_width)][0] << "	";
		}
		//abcdef << endl;
	}
	abcdef << endl;
	abcdef.close();
	final_popn.resize(sqst);
	vector<double> x;
	x.resize(sqst);
	for(int s=0; s<sqst; s++) x[s] = Average(s); // Use Region::Average to get the average density of the file (this will be the new density of the region)
	fine = 0, width=1, elements=1, extended_elements=9, extended_width=3; // Update width, elements etc
	for(int s=0; s<sqst; s++){
		for(int n=0; n<width; n++) delete[](growth[s][n]); // Delete the HR array of growth rates
		delete[](growth[s]); // Delete the HR array of growth rates
		fftw_free(distribution[s]); // Free the HR distribution
		fftw_free(fourier_distribution[s]); // Free the HR fourier-distribution
		distribution[s] = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*9); // Allocate 3x3 LR distribution for the region and its 8 nearest neighbours
		for(int n=0; n<extended_elements; n++) distribution[s][n][0]=0, distribution[s][n][1]=0; // Assign zero values to neighbours in distribution
		distribution[s][extended_width+1][0]=x[s]; // Assign average density x to the region's denisty in distribution
		cout << x[s] << "	";
		growth[s] = new double*[1]; // Allocate new 1x1 grotwh array
		growth[s][0] = new double[1]; // Allocate new 1x1 grotwh array
		growth[s][0][0] = growth_rate[s]; // Assign the average growth rate as growth-rate
		final_popn[s]=x[s]; // Final population is the average density as the region enters the wave-back (becomes coarse)
	}
}
#endif
