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



int main(int argc, const char * argv[]) {
    // insert code here...

    vector<double>* VGS = new vector<double>;
    vector<double>* VDS = new vector<double>;
    vector<double>* IDS = new vector<double>;   
    
    //string path = "/Users/arianabruno/Desktop/ECE4960/ProgrammingAssignments/ECE4960-PA3/outputNMOS.txt";
	string path = "C:/Users/Haritha/Documents/ECE4960-PAs/ECE4960-PA3/outputNMOS.txt";
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
    return 0;
}
