#include <iostream>
#include <fstream>
#include <string>

#include "fullSolver.hpp"

void readDataFile( string path , vector<double>* VGS , vector<double>* VDS , vector<double>* IDS );
double sum_of_squares( vector<double>* s_model , vector<double>* s_measured );
void search_values( vector<double>* s_measured , vector<double>* VGS , vector<double>* VDS , vector<double>* IDS , double vgs , double vds );