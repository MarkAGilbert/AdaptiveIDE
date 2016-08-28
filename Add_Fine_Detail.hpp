#include<fstream>
#include <string>
#include <sstream>
#include "Parameters.hpp"

using namespace std;

void Add_Fine_Detail(Parameters P){
	//ifstream f("temporary_file.bin");
    //if (f.good()) {

	double** map;
	map = new double*[P.landscape_length*P.region_width*P.T/P.print_size]; // Allocate the growth array for the (non-extended) domain
	for(int n=0; n<P.landscape_length*P.region_width*P.T/P.print_size; n++) map[n] = new double[P.landscape_width*P.region_width/P.print_size];

	fstream file (P.file, ios::in|ios::binary); // Create fstream for map_file
	file.seekg (0, file.end);
	int length = file.tellg();
	file.seekg (0, file.beg);
	char * buffer = new char [length];
	file.read(buffer,length);
	stringstream stream(buffer);
	for(int ind1=0; ind1<P.landscape_length*P.region_width*P.T/P.print_size; ind1++){
		for(int ind2=0; ind2<P.landscape_width*P.region_width/P.print_size; ind2++){
			double n;
			stream >> n;
			map[ind1][ind2] = n;
		}
	}
	file.close();

/*		double** fine_to_coarse_map;
		fine_to_coarse_map = new double*[P.landscape_length*P.region_width/P.print_size];
		for(int n=0; n<P.landscape_length*P.region_width/P.print_size; n++) map[n] = new double[P.landscape_width*P.region_width/P.print_size];
		fstream file2 ("temporary_file.bin", ios::in|ios::binary); // Create fstream for map_file
		file2.seekg (0, file2.end);
		length = file2.tellg();
		file2.seekg (0, file2.beg);
		char * buffer2 = new char [length];
		file2.read(buffer2,length);
		stringstream stream2(buffer2);
		file2.close();
		while(stream2.peek() != EOF){
			int n1, n2;
			stream >> n1;
			stream >> n2;
			for(int n=0; 
				*/



	for(int n=P.landscape_length*P.region_width/P.print_size; n<P.landscape_length*P.region_width*P.T/P.print_size; n++){
		for(int m=0; m<P.landscape_width*P.region_width/P.print_size; m++){
			if(map[n][m]<-.5){
				map[n][m] = map[n-P.landscape_length*P.region_width/P.print_size][m];
			}
		}
	}

	ofstream abcdef;
	abcdef.open(P.file, ios::out | ios::binary);
	for(int n=0; n<P.landscape_length*P.region_width*P.T/P.print_size; n++){
		for(int m=0; m<P.landscape_width*P.region_width/P.print_size; m++){
			abcdef << map[n][m] << "	";
		}
		abcdef << endl;
	}
		
	/*	ofstream abcde;
		abcde.open("temporary_file.bin", ios::out | ios::binary);
	}*/
}