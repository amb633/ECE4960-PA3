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
#include "fullSolver.hpp"
#include "utilityFunctions.hpp"

int main(int argc, const char * argv[]) {
    // insert code here...
    fullSolver();

    vector<double>* VGS = new vector<double>;
    vector<double>* VDS = new vector<double>;
    vector<double>* IDS = new vector<double>;   
    
    string path = "/Users/arianabruno/Desktop/ECE4960/ProgrammingAssignments/ECE4960-PA3/outputNMOS.txt";
   // string path = "C:/Users/Haritha/Documents/ECE4960-PAs/ECE4960-PA3/outputNMOS.txt";
    readDataFile( path , VGS ,  VDS ,  IDS );

    vector<double>* s_measured = new vector<double>;
    search_values( s_measured , VGS , VDS , IDS , 1.0 , 1.95 );
    printMatrix( s_measured);
    return 0;
}