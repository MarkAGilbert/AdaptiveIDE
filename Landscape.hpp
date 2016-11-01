#ifndef LANDSCAPEHEADERDEF
#define LANDSCAPEHEADERDEF

#include<vector>
#include "Region.hpp"
#include "Kernel.hpp"
#include "Parameters.hpp"
#include<iostream>
#include<fstream>

using namespace std;

class Landscape
{

public:

// Constructor
Landscape(Parameters Param); // Constructor (no need for copy-constructor, as only one Landscape)

// Destructor
~Landscape(); // Destructor for the Landscape

// Functions
void Growth(); // Growth phase for the Landscape
void Dispersal(); // Dispersal phase for the Landscape
void Sum(); // Sums up stage-pairs after dispersal
void Remesh(); // Remesh phase for the Landscape
void Print(int t); // Prints distribution to file

// Variables
Parameters P; // Simulation Parameters
int t; // Number of generations that have elapsed
vector<vector<Region> > Popn; // Array of regions (each region contains its own population distribution, growth rates etc)

vector<fftw_complex*> kernel; //
vector<fftw_plan> p_gforward, p_gbackward;
vector<fftw_complex*> gkernel;
vector<double*> lckernel;
vector<fftw_complex*> RegionalAverages; // Averages in physical space
vector<double> OldAve;
vector<fftw_complex*> RA; // Averages in Fourier space
double** initial_cells;
double** initial_regions;
double** initial_regions_neighbours;
};

Landscape::Landscape(Parameters Param){

  // Set landscape parameters
  P = Param;
  t = 0;

  // Build Population (Popn) Array
  Popn.resize(P.landscape_length);
  for(int n1=0; n1<P.landscape_length; n1++){
     Popn[n1].resize(P.landscape_width,Region(P.stages));
     for(int n2=0; n2<P.landscape_width; n2++){
       Popn[n1][n2].GrowthRate(P.intrinsic_growth_rate, P.map_file, P.map_type, n1, n2, P.region_width, P.landscape_length);
     }
  }

  ofstream abcdef;
  abcdef.open("temporary_file.bin", ios::out |  ios::binary);
  abcdef.close();
  
  if(P.initial_popn_map == 0){
	  // Set Initial Values
	  int ind2 = P.initial_popn_loc_x/P.region_width;
	  int ind1 = P.initial_popn_loc_y/P.region_width;
	  int loc2 = P.initial_popn_loc_x % P.region_width;
	  int loc1 = P.initial_popn_loc_y % P.region_width;
	  for(int k1 = max(0,ind1-1); k1<min(P.landscape_width,ind1+2); k1++){
		  for(int k2 = max(0,ind2-1); k2<min(P.landscape_length,ind2+2); k2++){
			  Popn[k1][k2].Fine(P, P.region_width, P.map_file, "bool", k1, k2, P.landscape_length, P.intrinsic_growth_rate);
		  }
	  }
	  Popn[ind1][ind2].distribution[0][4*P.region_width*(P.region_width+loc2)+P.region_width+loc1][0] = .01;
	  cout << loc1 << "	" << loc2 << endl << endl;
  }
  if(P.initial_popn_map == 1){
	  initial_cells = new double*[P.landscape_length*P.region_width];
	  for(int n=0; n<P.landscape_length*P.region_width; n++) initial_cells[n] = new double[P.landscape_width*P.region_width];
	  for(int n=0; n<P.landscape_length*P.region_width; n++) for(int m=0; m<P.landscape_width*P.region_width; m++) initial_cells[n][m]=0;
	  initial_regions = new double*[P.landscape_length];
	  for(int n=0; n<P.landscape_length; n++) initial_regions[n] = new double[P.landscape_width];
	  for(int n=0; n<P.landscape_length; n++) for(int m=0; m<P.landscape_width; m++) initial_regions[n][m] = 0;
	  initial_regions_neighbours = new double*[P.landscape_length];
	  for(int n=0; n<P.landscape_length; n++) initial_regions_neighbours[n] = new double[P.landscape_width];
	  for(int n=0; n<P.landscape_length; n++) for(int m=0; m<P.landscape_width; m++) initial_regions_neighbours[n][m] = 0;

	
	  fstream file (P.map_initial_file, ios::in|ios::binary);
	  bool read;
	  file.seekg(0);
	  for(int n=0; n<P.landscape_length*P.region_width; n++){
		  for(int m=0; m<P.landscape_width*P.region_width; m++){
			  file.read(reinterpret_cast<char*>( &read ), sizeof read );
			  initial_cells[n][m] = int(read);
			  if(read==1){
				  int reg_n = n/P.region_width;
				  int reg_m = m/P.region_width;
				  initial_regions[reg_n][reg_m] = 1;
			  }
		  }
	  }
	  file.close();
	  
	  for(int n=0; n<P.landscape_length; n++){
		  for(int m=0; m<P.landscape_width; m++){
			  if(initial_regions[n][m]==1){
				  for(int k1 = max(0,n-1); k1<min(P.landscape_length,n+2); k1++){
					  for(int k2 = max(0,m-1); k2<min(P.landscape_length,m+2); k2++){
						  initial_regions_neighbours[k1][k2] = 1;
					  }
				  }
			  }
		  }
	  }
	  for(int n=0; n<P.landscape_length; n++){
		  for(int m=0; m<P.landscape_width; m++){
			  if(initial_regions_neighbours[n][m] == 1){
				  Popn[n][m].Fine(P, P.region_width, P.map_file, "bool", n, m, P.landscape_length, P.intrinsic_growth_rate);
			  }
			  if(initial_regions[n][m]==1){
				  for(int k1=0; k1<P.region_width; k1++){
					  for(int k2=0; k2<P.region_width; k2++){
						  if(initial_cells[P.region_width*n + k1][P.region_width*m + k2]==1){
							  Popn[n][m].distribution[0][4*P.region_width*(P.region_width+k2)+P.region_width+k1][0] = 1.0;
						  }
					  }
				  }
			  }
		  }
	  }
  }
  
  RegionalAverages.resize(P.stsq);
  RA.resize(P.stsq);
  for(int n1=0; n1<P.stsq; n1++){
     RegionalAverages[n1] = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*4*P.landscape_length*P.landscape_width);
	 RA[n1] = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*4*P.landscape_length*P.landscape_width);
	 for(int n2=0; n2<4*P.landscape_length*P.landscape_width; n2++) RegionalAverages[n1][n2][0]=0.0, RegionalAverages[n1][n2][1]=0.0;
  }
  

  // Build Dispersal Kernels
  kernel.resize(P.stsq);
  gkernel.resize(P.stsq);
  lckernel.resize(P.stsq);
  for(int n=0; n<P.stsq; n++){
	kernel[n] = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*pow(4.0*P.region_width,2));
	gkernel[n] = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*4*P.landscape_length*P.landscape_width);
	lckernel[n] = new double[9];
  }
  Local_Kernel(P, kernel, lckernel);
  Global_Kernel(P, gkernel);
  p_gforward.resize(P.stsq);
  p_gbackward.resize(P.stsq);
 
}

Landscape::~Landscape()
{
	for(int s=0; s<P.stsq; s++){
	   fftw_free(kernel[s]);
	   fftw_free(gkernel[s]);
	   fftw_free(RegionalAverages[s]);
	   fftw_free(RA[s]);
	   delete[](lckernel[s]);
	}
}

void Landscape::Growth(){
  for(int n1=0; n1<P.landscape_length; n1++){
    for(int n2=0; n2<P.landscape_width; n2++){
      Popn[n1][n2].Grow();
    }
  }
}


void Landscape::Dispersal(){
	
	int n[2], N[2];
	N[0]=P.landscape_length, N[1] = P.landscape_width;

	for(int s=0; s<P.stsq; s++){
		if(P.dispersal_on[s]==1){
	  
		  // Global Dispersal 1 - Get The Average Populations from each Region
		  for(int n=0; n<4*N[0]*N[1]; n++) RegionalAverages[s][n][0] = 0, RegionalAverages[s][n][1] = 0;
		  for(int n0=0; n0<N[0]; n0++){
			  for(int n1=0; n1<N[1]; n1++){
				  RegionalAverages[s][n0+2*P.landscape_length*n1][0] = Popn[n0][n1].Average(s);
				  RegionalAverages[s][n0+2*P.landscape_length*n1][1] = 0;
			  }
		  }

		  // Global Dispersal 2 - Do the Global Dispersal
		  p_gforward[s]=fftw_plan_dft_1d(4*P.landscape_length*P.landscape_width,RegionalAverages[s],RA[s],1,FFTW_ESTIMATE);
		  fftw_execute(p_gforward[s]);
		  for(int m=0; m<4*P.landscape_length*P.landscape_width;m++){
			double a = RA[s][m][0]/(4*P.landscape_length*P.landscape_width);
			double b = RA[s][m][1]/(4*P.landscape_length*P.landscape_width);
			RA[s][m][0] = a*gkernel[s][m][0] - b*gkernel[s][m][1];
			RA[s][m][1] = a*gkernel[s][m][1] + b*gkernel[s][m][0];
		  }
		  p_gbackward[s] = fftw_plan_dft_1d(4*P.landscape_length*P.landscape_width,RA[s],RegionalAverages[s],-1,FFTW_ESTIMATE);
		  fftw_execute(p_gbackward[s]);
		  fftw_destroy_plan(p_gforward[s]);
		  fftw_destroy_plan(p_gbackward[s]);
		}
	}
	
	// Local Dispersal 1
	int a[2], ns[2], ne[2];
	for(a[0]=-1; a[0]<2; a[0]++){
		for(a[1]=-1; a[1]<2; a[1]++){
			for(int j=0; j<2; j++){
				if(a[j]==-1) ns[j] = 0, ne[j] = N[j]-1;
				else if(a[j]==0) ns[j] = 0, ne[j] = N[j];
				else ns[j] = 1, ne[j] = N[j];
			}
			if(a[0] != 0 || a[1] != 0){
				for(int n0 =ns[0]; n0 < ne[0]; n0++){
					for(int n1 =ns[1]; n1 <ne[1]; n1++){
						Popn[n0][n1].Share(a[0],a[1],&Popn[n0-a[0]][n1-a[1]]);
					}
				}
			}
		}
	}
	
	// Local Dispersal 2 - Do the Local Dispersal
	for(n[0]=0; n[0]<N[0]; n[0]++){
		for(n[1]=0; n[1]<N[1]; n[1]++){
			if(Popn[n[0]][n[1]].fine==1) Popn[n[0]][n[1]].Disperse(kernel,P);
			else Popn[n[0]][n[1]].Disperse(lckernel,P);
		}
	}
	
	for(int s=0; s<P.stsq; s++){
	  // Add global to local
	  int m[2];
	  for(m[0]=0; m[0]<P.landscape_length; m[0]++){
		for(m[1]=0; m[1]<P.landscape_width; m[1]++){
			Popn[m[0]][m[1]].Add(RegionalAverages[s][m[0]+2*P.landscape_length*m[1]][0],s);
		}
	  }
	}
}

void Landscape::Sum(){
  for(int n1=0; n1<P.landscape_length; n1++){
    for(int n2=0; n2<P.landscape_width; n2++){
      Popn[n1][n2].Sum();
    }
  }
}

void Landscape::Remesh(){
    int m[2];
    for(m[0]=0; m[0]<P.landscape_length; m[0]++){
        for(m[1]=0; m[1]<P.landscape_width; m[1]++){
            if(Popn[m[0]][m[1]].used ==0 && Popn[m[0]][m[1]].fine ==0 && Popn[m[0]][m[1]].Average(0)<=P.threshcoarse){
                for(int a=-1; a<2; a++){
                    for(int b=-1; b<2; b++){
                        if(m[0]+a >= 0 && m[0]+a < P.landscape_length && m[1]+b >=0 && m[1]+b < P.landscape_width){
                            if(Popn[m[0]][m[1]].fine ==0 && Popn[m[0]+a][m[1]+b].Average(0)>P.threshfine) Popn[m[0]][m[1]].Fine(P, P.region_width, P.map_file, P.map_type, m[0], m[1], P.landscape_length, P.intrinsic_growth_rate);
                        }
                    }
                }
            }
            if(Popn[m[0]][m[1]].fine == 1 && Popn[m[0]][m[1]].Minimum()>P.threshcoarse) Popn[m[0]][m[1]].Coarse(m[0],m[1],P);
            if(Popn[m[0]][m[1]].fine == 1 && t>0 && fabs(Popn[m[0]][m[1]].Average(0)-Popn[m[0]][m[1]].OldAve[0])<P.threshcoarse2 && Popn[m[0]][m[1]].Average(0)>P.threshcoarse3){ 
				Popn[m[0]][m[1]].Coarse(m[0],m[1],P);
				cout << endl << endl <<  "Wooht" << endl << endl;
			}
			for(int s=0; s<P.stsq; s++) Popn[m[0]][m[1]].OldAve[s] = Popn[m[0]][m[1]].Average(s);
        }
    }
    t++;
}

void Landscape::Print(int t){
    if(P.print == 1){
        // Open Output File
        ofstream abcdef;
        const char* file = P.file;
		if(t==0){
			abcdef.open(file, ios::out | ios::binary);
		}
		else abcdef.open(file, ios::out | ios::app | ios::binary);
		if((P.T-1-t)%P.print_frequency == 0){
			cout << t << endl;
			for(int i=0; i<P.landscape_length; i++){
				for(int n1=0; n1<P.region_width; n1+=P.print_size){
					for(int s=0; s<P.stages; s++){
						for(int j=0; j<P.landscape_width; j++){
							for(int n2=0; n2<P.region_width; n2+=P.print_size){
								if(Popn[i][j].fine!=0){
									double a = 0;
									for(int m1=0; m1<P.print_size; m1++) for(int m2=0; m2<P.print_size; m2++) a+= Popn[i][j].distribution[s*(P.stages+1)][P.region_width+n1+m1+4*P.region_width*(n2+m2+P.region_width)][0]/pow(1.0*P.print_size,2);
									int m1=0;
									int m2=0;
									//cout << Popn[i][j].distribution[0][P.region_width+n1+m1+4*P.region_width*(n2+m2+P.region_width)][0]/pow(1.0*P.print_size,2) << "	" << Popn[i][j].distribution[0][P.region_width+n1+m1+4*P.region_width*(n2+m2+P.region_width)][0] << endl;
									abcdef << a << "    ";
								}
								//else if(Popn[i][j].used==1) abcdef << -1 << "	";
								else abcdef << Popn[i][j].distribution[s*(P.stages+1)][4][0] << "    ";
							}
						}
					}
				abcdef << endl; // New Line For Each 'Row' of the Landscape
				}
			}
		}
        abcdef << endl; // Additional New Line For Each Time-Step
		abcdef.close();
    }
}

#endif

