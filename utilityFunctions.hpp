#include <iostream>
#include <fstream>
#include <string>
#include <cmath>

#include "fullSolver.hpp"

void readDataFile( string path , vector<double>* VGS , vector<double>* VDS , vector<double>* IDS );
double sumSquares( vector<double>* s_model , vector<double>* s_measured );
void searchValues( vector<double>* s_measured , vector<double>* VGS , vector<double>* VDS , vector<double>* IDS , double vgs , double vds );
double calculateIds( double vgs , double vds , double kappa , double vth , double is );
void modelIds( vector<double>* IDS_model , vector<double>* VGS , vector<double>* VDS , double kappa , double vth , double is );
void add_vectors( vector<double>* x,  vector<double>* dx,  vector<double>* sum );
double delta_norm_2( vector<double>* delta_a, vector<double>* a );
