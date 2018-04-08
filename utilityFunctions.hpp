#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <cstdlib>

#include "fullSolver.hpp"
using namespace std;

void readDataFile( string path , vector<double>* VGS , vector<double>* VDS , vector<double>* IDS );
double sumSquares( vector<double>* s_model , vector<double>* s_measured , bool n = false );
void searchValues( vector<double>* s_measured , vector<double>* VGS , vector<double>* VDS , vector<double>* IDS , double vgs , double vds );
double calculateIds( double vgs , double vds , double kappa , double vth , double is );
void modelIds( vector<double>* IDS_model , vector<double>* VGS , vector<double>* VDS , double kappa , double vth , double is );
void add_vectors( vector<double>* x,  vector<double>* dx,  vector<double>* sum );
double delta_norm_rel( vector<double>* delta_a, vector<double>* a );
double delta_norm_abs( vector<double>* delta_a);
double parameterSensitivity( vector<double>* parameters , double pertubation , int which ,
    vector<double>* VGS , vector<double>* VDS , vector<double>* IDS );
void scaleVector( double scalar, vector<double>* a,  vector<double>* result);