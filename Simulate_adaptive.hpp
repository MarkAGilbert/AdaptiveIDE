#include<iostream>
#include<fstream>
#include "fftw3.h"
#include "Region.hpp"
#include "Landscape.hpp"
#include "Parameters.hpp"


using namespace std;

void Simulate_adaptive(Parameters P)
{
    // Initialise Landscape
    Landscape Britain(P);

    // Run Simulation
    for(int t=0; t<P.T; t++){
		
        Britain.Print(t); // Print Population Distribution Data to Output File (Every P.print_frequency generations)
		
        Britain.Growth(); // Growth

        Britain.Dispersal(); // Dispersal

		Britain.Sum(); // Sums up stage-pairs after dispersal

        Britain.Remesh(); // Remesh Landscape
		
    }

}
