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

void readDataFile( string path , vector<double>* VGS , vector<double>* VDS , vector<double>* IDS );

int main(int argc, const char * argv[]) {
    // insert code here...
    fullSolver();

    vector<double>* VGS = new vector<double>;
    vector<double>* VDS = new vector<double>;
    vector<double>* IDS = new vector<double>;   
    
    //string path = "/Users/arianabruno/Desktop/ECE4960/ProgrammingAssignments/ECE4960-PA3/outputNMOS.txt";
    string path = "C:/Users/Haritha/Documents/ECE4960-PAs/ECE4960-PA3/outputNMOS.txt";
    readDataFile( path , VGS ,  VDS ,  IDS );
   
    return 0;
}

void readDataFile( string path , vector<double>* VGS , vector<double>* VDS , vector<double>* IDS ){
    // function to read and store file
    ifstream myfile( path );
    if (myfile.is_open()){
        myfile.ignore(2048 , '\n');
        while( !myfile.eof() ){
            double vgs , vds , ids;
            myfile >> vgs >> vds >> ids;
            (*VGS).push_back(vgs);
            (*VDS).push_back(vds);
            (*IDS).push_back(ids);
        }
        myfile.close();
        
        // note that the last entry is read twice when using eof
        // so use pop_back to get rid of extra entry
        (*VGS).pop_back();
        (*VDS).pop_back();
        (*IDS).pop_back();
    }
    else cout << "Unable to open file " << endl;
    return;
}