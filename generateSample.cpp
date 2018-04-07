//
//  generateSample.cpp
//  PA3
//
//  Created by Ariana Bruno on 4/7/18.
//  Copyright © 2018 Ariana Bruno. All rights reserved.
//

#include "generateSample.hpp"

void randomSamples(vector<double>* x, vector<double>* y, vector<double>* noise, vector<double>* y_n){
    //    the model parameters
    double c_0 = 10.0;
    double m = -0.5;
    
    //    For generating different random samples each time function is run
    //    srand(time(0));
    for( int i = 0; i<10; i++){
        
        //        generating random values for x input
        double x_value = rand()/100000.0;
        (*x).push_back(x_value);
        
        //        calculating the actual y value for the power law model
        double y_value = c_0*pow(x_value,m);
        (*y).push_back(y_value);
        
        //        calculating the noise to add to the y model result
        double rand_noise = (rand()%100000)/100000.0;    //calculate random values from 0 to 1
        
        if (rand_noise < 0.50){
            //        values less than 0.5 will be negative 0.1 to 0.2 noise
            rand_noise = -(rand_noise/5.0) - 0.1;
        }else{
            //        values greater than or equal to 0.5 will be positive 0.1 to 0.2 noise
            rand_noise = (rand_noise/10.0) + 0.1;
        }
        (*noise).push_back(rand_noise);
        
        //        adding noise to the y model result
        (*y_n).push_back(y_value + rand_noise*y_value);
    }
    
    //    print the values for verification that function is working properly
    //    printMatrix(x_samples);
    //    printMatrix(y_samples);
    //    printMatrix(y_noisey_samples);
    //    printMatrix(noise_samples);
    
}

void LinearLSF(vector<double>* x, vector<double>* y, vector<vector<double>*>* H, vector<double>* RHS){
    //    initialize the Hessian matrix
    zeroMatrix(H, 2);
    
    //    initialize each term for the Hessian matrix and RHS matrix
    double sum_xi = 0.0, sum_xi_sq = 0.0, sum_bs = 0.0, sum_yi = 0.0, sum_xi_yi = 0.0;
    
    //    check that the size of the x and y vectors are the same (ensure each value has a pair)
    if( x->size() != y->size() ){
        cout << "The x and y pairs are not the same size... cannot calculate LSF ";
    } else{
        
        //        calculate the summed values for H and RHS
        for(int i = 0; i<x->size(); i++){
            // need to take the log of x and y because we are treating the power law as a linear function
            sum_xi = sum_xi + log((*x)[i]);
            sum_bs = sum_bs + 1;
            sum_xi_sq = sum_xi_sq + pow(log((*x)[i]),2);
            sum_yi = sum_yi + log((*y)[i]);
            sum_xi_yi = sum_xi_yi + (log((*x)[i])*log((*y)[i]));
        }
        
        (*RHS).push_back(sum_xi_yi);
        (*RHS).push_back(sum_yi);
        
        (*(*H)[0])[0] = sum_xi_sq;
        (*(*H)[0])[1] = sum_xi;
        (*(*H)[1])[0] = sum_xi;
        (*(*H)[1])[1] = sum_bs;
    }
    
}
