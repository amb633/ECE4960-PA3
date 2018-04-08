//
//  main.cpp
//  PA3
//
//  Created by Ariana Bruno on 4/7/18.
//  Copyright © 2018 Ariana Bruno. All rights reserved.
//

#include <iostream>
#include <fstream>
#include <string>
#include <ctime>
#include <cstdlib>
#include "fullSolver.hpp"
#include "generateSample.hpp"
#include "utilityFunctions.hpp"
#include "quasiNewtonMethod.hpp"

vector<double>* VGS = new vector<double>;
vector<double>* VDS = new vector<double>;
vector<double>* IDS = new vector<double>;

string path = "/Users/arianabruno/Desktop/ECE4960/ProgrammingAssignments/ECE4960-PA3/outputNMOS.txt";
//string path = "C:/Users/Haritha/Documents/ECE4960-PAs/ECE4960-PA3/outputNMOS.txt";


int main(int argc, const char * argv[]) {
    // insert code here...
    
    readDataFile( path , VGS ,  VDS ,  IDS );
    
    vector<double>* x_samples  = new vector<double>;
    vector<double>* y_samples  = new vector<double>;
    vector<double>* noise_samples  = new vector<double>;
    vector<double>* y_noisey_samples  = new vector<double>;
    
    randomSamples(x_samples, y_samples, noise_samples, y_noisey_samples);
    
    vector<vector<double>*>* H_matrix = new vector<vector<double>*>;
    vector<double>* RHS = new vector<double>;
    
//    find the H matrix and RHS vector for this problem for parameter extraction of  m and c0 for the S_measured data
    LinearLSF(x_samples, y_noisey_samples, H_matrix, RHS);
    vector<double>* solution = new vector<double>;

//    solve for the parameters given H and RHS
    fullSolver( solution , H_matrix , RHS );
    
    cout << "This is the solution for m and log(c0)" << endl;
    printMatrix( solution );
    cout << "m = " << (*solution)[0] << endl;
    cout << "c0 = " << exp((*solution)[1]) << endl;
    
    cout << endl << "Quasi Newton Parameter Extraction" << endl << endl;
    double k_0 = 1.0;
    double Vth_0 = 1.0;
    double Is_0 = 1e-7;
    
    vector<double>* current_parameters = new vector<double>;
    (*current_parameters).push_back(k_0); (*current_parameters).push_back(Vth_0); (*current_parameters).push_back(Is_0);
    vector<double>* updated_parameters = new vector<double>;
    
    double norm_V, norm_delta;
    
    quasiNetwon_itr(VGS, VDS, IDS, current_parameters, updated_parameters, &norm_V, &norm_delta);
    
    if( norm_V > 1e-4){
        (*current_parameters) = (*updated_parameters);
        vector<double>* updated_parameters = new vector<double>;
        quasiNetwon_itr(VGS, VDS, IDS, current_parameters, updated_parameters, &norm_V, &norm_delta);
    }
    
    
    return 0;
}
