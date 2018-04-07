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
#include "utilityFunctions.hpp"

void LinearLSF(vector<double>* x, vector<double>* y, vector<vector<double>*>* H, vector<double>* RHS){
    zeroMatrix(H, 2);
    double sum_xi = 0.0, sum_xi_sq = 0.0, sum_bs = 0.0, sum_yi = 0.0, sum_xi_yi = 0.0;
    if( x->size() != y->size() ){
        cout << "The x and y pairs are not the same size... cannot calculate LSF ";
    } else{
        for(int i = 0; i<x->size(); i++){
            sum_xi = sum_xi + (*x)[i];
            sum_bs = sum_bs + 1;
            sum_xi_sq = sum_xi_sq + pow((*x)[i],2);
            sum_yi = sum_yi + (*y)[i];
            sum_xi_yi = sum_xi_yi + ((*x)[i]*(*y)[i]);
        }
        
        (*RHS).push_back(sum_xi_yi);
        (*RHS).push_back(sum_yi);
        
        (*(*H)[0])[0] = sum_xi_sq;
        (*(*H)[0])[1] = sum_xi;
        (*(*H)[1])[0] = sum_xi;
        (*(*H)[1])[1] = sum_bs;
    }
    
    
    
}

int main(int argc, const char * argv[]) {
    // insert code here...

    vector<double>* VGS = new vector<double>;
    vector<double>* VDS = new vector<double>;
    vector<double>* IDS = new vector<double>;   
    
    //string path = "/Users/arianabruno/Desktop/ECE4960/ProgrammingAssignments/ECE4960-PA3/outputNMOS.txt";
    string path = "C:/Users/Haritha/Documents/ECE4960-PAs/ECE4960-PA3/outputNMOS.txt";
    readDataFile( path , VGS ,  VDS ,  IDS );

    vector<double>* s_measured = new vector<double>;
    search_values( s_measured , VGS , VDS , IDS , 1.0 , 1.95 );
    //printMatrix( s_measured);
    
    double c_0 = 10.0;
    double m = -0.5;
    vector<double>* x_samples  = new vector<double>;
    vector<double>* y_samples  = new vector<double>;
    vector<double>* noise_samples  = new vector<double>;
    vector<double>* y_noisey_samples  = new vector<double>;
    
//    srand(time(0));
    for( int i = 0; i<10; i++){
        
        double x = rand()/100000.0;
        (*x_samples).push_back(x);
        double y = c_0*pow(x,m);
        (*y_samples).push_back(y);
        double rand_noise = (rand()%100000)/100000.0;
        if (rand_noise < 0.50){
            rand_noise = -(rand_noise/5.0) - 0.1;
        }else{
            rand_noise = (rand_noise/10.0) + 0.1;
        }
        (*noise_samples).push_back(rand_noise);
        (*y_noisey_samples).push_back(y + rand_noise*y);
    }
    

    printMatrix(x_samples);
    printMatrix(y_samples);
    printMatrix(y_noisey_samples);
    printMatrix(noise_samples);
    
    vector<vector<double>*>* H_matrix = new vector<vector<double>*>;
    vector<double>* RHS = new vector<double>;
    LinearLSF(x_samples, y_noisey_samples, H_matrix, RHS);
    vector<double>* solution = new vector<double>;

    fullSolver( solution , H_matrix , RHS );
    printMatrix( solution );
    return 0;
}
