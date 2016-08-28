#ifndef PARAMETERSHEADERDEF
#define PARAMETERSHEADERDEF

using namespace std;

class Parameters
{
    public:

    Parameters();
    Parameters(const Parameters& Other); // Copy Constructor

	bool adaptive;
    double* dispersal_distance;
	int* dispersal_on;
    int landscape_width;
    int landscape_length;
    int region_width;
    double* intrinsic_growth_rate;
	double resolution;
    int kernel_no;
	int stages;
	int stsq;
    int T;
    double threshfine;
    double threshcoarse;
    double threshcoarse2;
    double threshcoarse3;
    double period;
    double length;
	int initial_popn_map;
	int initial_popn_loc_x;
	int initial_popn_loc_y;
	const char * map_initial_file;
    const char * map_type;
    const char * map_file; // Location of map

    // Outputs
    int print;
    int print_frequency;
    bool dist_all;
    bool area_all;
    bool time_all;
    bool dist_end;
    bool area_end;
    bool time_end;
    int print_size; // Should divide the region_width
    const char * file; // Where to print the distributions to
	const char * file_additional;
	int smooth_output;
	
    void Values();
    void Default();
};

Parameters::Parameters(){
    char lastword[200];
    char oneword[200];
    lastword[0] = 'X';
    int c;
    double x,y;
    FILE *fp1;
	fopen_s(&fp1,"parameters.txt", "r");
    do {
        c = fscanf_s(fp1,"%s",oneword,200); 
		if(strcmp(lastword,"initial_popn_map:")==0) initial_popn_map = atoi(oneword);
		if(strcmp(lastword,"initial_popn_loc_pixel_x:")==0) initial_popn_loc_x = atoi(oneword)-1;
		if(strcmp(lastword,"initial_popn_loc_pixel_y:")==0) initial_popn_loc_y = atoi(oneword)-1;
        if(strcmp(lastword,"kernel_no:")==0) kernel_no = atoi(oneword);
        if(strcmp(lastword,"resolution:")==0) resolution = atof(oneword);
        if(strcmp(lastword,"print:")==0) print = atoi(oneword);
		if(strcmp(lastword,"stages:")==0){ 
			stages = atoi(oneword);
			stsq = stages*stages;
			intrinsic_growth_rate = new double[stsq];
			dispersal_distance = new double[stsq];
			dispersal_on = new int[stsq];
		}
		if(strcmp(lastword,"print_frequency:")==0) print_frequency = atoi(oneword);
        if(strcmp(lastword,"output_resolution:")==0) print_size = atoi(oneword);
        if(strcmp(lastword,"T:")==0) T = atoi(oneword);
        if(strcmp(lastword,"threshcoarse2:")==0) threshcoarse2 = atof(oneword);
        if(strcmp(lastword,"threshcoarse3:")==0) threshcoarse3 = atof(oneword);
		if(strcmp(lastword,"region_width_in_pixels:")==0) region_width = atoi(oneword);
        if(strcmp(lastword,"landscape_width_in_pixels:")==0) landscape_width = atoi(oneword)/region_width;
        if(strcmp(lastword,"landscape_length_in_pixels:")==0) landscape_length = atoi(oneword)/region_width;
        if(strcmp(lastword,"threshcoarse:")==0) threshcoarse = atof(oneword);
        if(strcmp(lastword,"threshfine:")==0) threshfine = atof(oneword);
        if(strcmp(lastword,"intrinsic_growth_rate:")==0){ 
			for(int n=0; n<stsq; n++){ 
				if(n>0) c = fscanf_s(fp1,"%s",oneword,200); 
				intrinsic_growth_rate[n] = atof(oneword);
			}
		}
		if(strcmp(lastword,"dispersal_distance:")==0){
			for(int n=0; n<stsq; n++){
				if(n>0) c = fscanf_s(fp1,"%s",oneword,200);
				dispersal_distance[n] = atof(oneword);
			}
		}
		if(strcmp(lastword,"dispersal_on:")==0){
			for(int n=0; n<stsq; n++){
				if(n>0) c = fscanf_s(fp1,"%s",oneword,200);
				dispersal_on[n] = atoi(oneword);
			}
		}
		if(strcmp(lastword,"smooth_output:")==0) smooth_output = atoi(oneword);
        file = "output.bin";
        map_file = "new_forest.bin";
		map_initial_file = "line.bin";
		file_additional = "temporary_file.bin";
        map_type =  "bool";
        strcpy_s(lastword,oneword);
    } while (c != EOF);
	length = landscape_width * region_width * resolution;
}

// Copy Constructor
Parameters::Parameters(const Parameters& Other){
	initial_popn_loc_x = Other.initial_popn_loc_x;
	initial_popn_loc_y = Other.initial_popn_loc_y;
	initial_popn_map = Other.initial_popn_map;
	map_initial_file = Other.map_initial_file;
    adaptive = Other.adaptive;
	stages = Other.stages;
	stsq = Other.stsq;
    kernel_no = Other.kernel_no;
    length = Other.length;
    period = Other.period;
    print = Other.print;
    file = Other.file;
	intrinsic_growth_rate = new double[stsq];
    for(int n=0; n<stsq; n++) intrinsic_growth_rate[n] = Other.intrinsic_growth_rate[n];
	dispersal_distance  = new double[stsq]; 
	for(int n=0; n<stsq; n++) dispersal_distance[n] = Other.dispersal_distance[n];
	dispersal_on = new int[stsq];
	for(int n=0; n<stsq; n++) dispersal_on[n] = Other.dispersal_on[n];
    map_file = Other.map_file;
	resolution = Other.resolution;
    map_type = Other.map_type;
    print_frequency = Other.print_frequency;
    print_size = Other.print_size;
    T = Other.T;
    threshcoarse2 = Other.threshcoarse2;
    time_all = Other.time_all;
    threshcoarse3 = Other.threshcoarse3;
    landscape_width = Other.landscape_width;
    landscape_length = Other.landscape_length;
    region_width = Other.region_width;
    threshcoarse = Other.threshcoarse;
    threshfine = Other.threshfine;
	print_frequency = Other.print_frequency;
}

#endif
