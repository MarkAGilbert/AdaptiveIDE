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
    Simulate_adaptive(P);
    if(P.smooth_output==1) Add_Fine_Detail(P);

}
