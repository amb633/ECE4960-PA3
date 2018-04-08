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
//#include "utilityFunctions.hpp"
#include "quasiNewtonMethod.hpp"
#include "secant.hpp"

vector<double>* VGS = new vector<double>;
vector<double>* VDS = new vector<double>;
vector<double>* IDS = new vector<double>;

string path = "/Users/arianabruno/Desktop/ECE4960/ProgrammingAssignments/ECE4960-PA3/outputNMOS.txt";
//string path = "C:/Users/Haritha/Documents/ECE4960-PAs/ECE4960-PA3/outputNMOS.txt";


int main(int argc, const char * argv[]) {
    // insert code here...
   
    readDataFile( path , VGS ,  VDS ,  IDS );
    
    /*vector<double>* x_samples  = new vector<double>;
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
    cout << "c0 = " << exp((*solution)[1]) << endl; */
    
    cout << endl << " -------- Task 4A : Quasi Newton Parameter Extraction-------- " << endl << endl;
    
    double k_0 = 1.0;
    double Vth_0 = 1.0;
    double Is_0 = 1e-7;
    
    vector<double>* current_parameters = new vector<double>;
    (*current_parameters).push_back(k_0); (*current_parameters).push_back(Vth_0); (*current_parameters).push_back(Is_0);
    vector<double>* updated_parameters = new vector<double>;
    
    double norm_V, norm_delta_rel, norm_delta_abs;
    
    quasiNetwon_itr(VGS, VDS, IDS, current_parameters, updated_parameters, &norm_V, &norm_delta_rel, &norm_delta_abs, false);
    (*current_parameters) = (*updated_parameters);
    
    int itr = 0;
    while( norm_V > 1e-9 && norm_delta_rel > 1e-9){
        vector<double>* new_parameters = new vector<double>;
        quasiNetwon_itr(VGS, VDS, IDS, current_parameters, new_parameters, &norm_V, &norm_delta_rel, &norm_delta_abs, false);
        (*current_parameters) = (*new_parameters);
        itr++;
    }
    
    cout << "the converged solutions after " << itr << " iterations are: " << endl;
    cout << "kappa = " << (*current_parameters)[0] << endl;
    cout << "V_th = " << (*current_parameters)[1] << endl;
    cout << "Is = " << (*current_parameters)[2] << endl;
    cout << "absolute residual error = " << norm_delta_abs << endl;
    cout << "relative residual error = " << norm_delta_rel << endl;
    cout << "least squares = " << norm_V << endl;
    cout << endl;
    
    cout << endl << " -------- Task 5A : Quasi Newton Parameter Extraction-------- " << endl << endl;
    
    cout << "for normalized soltions: " << endl << endl;
    
    vector<double>* current_parameters_N = new vector<double>;
    (*current_parameters_N).push_back(k_0); (*current_parameters_N).push_back(Vth_0); (*current_parameters_N).push_back(Is_0);
    vector<double>* updated_parameters_N = new vector<double>;
    
    double norm_V_N, norm_delta_rel_N, norm_delta_abs_N;
    
    quasiNetwon_itr(VGS, VDS, IDS, current_parameters_N, updated_parameters_N, &norm_V_N, &norm_delta_rel_N, &norm_delta_abs_N, true);
    (*current_parameters_N) = (*updated_parameters_N);
    
    int itr_N = 0;
    while( norm_V_N > 1e-9 && norm_delta_rel_N > 1e-9){
        vector<double>* new_parameters_N = new vector<double>;
        quasiNetwon_itr(VGS, VDS, IDS, current_parameters_N, new_parameters_N, &norm_V_N, &norm_delta_rel_N, &norm_delta_abs_N, true);
        (*current_parameters_N) = (*new_parameters_N);
        itr_N++;
    }
    
    cout << "the converged solutions after " << itr_N << " iterations are: " << endl;
    cout << "kappa = " << (*current_parameters_N)[0] << endl;
    cout << "V_th = " << (*current_parameters_N)[1] << endl;
    cout << "Is = " << (*current_parameters_N)[2] << endl;
    cout << "absolute residual error = " << norm_delta_abs_N << endl;
    cout << "relative residual error = " << norm_delta_rel_N << endl;
    cout << "least squares = " << norm_V_N << endl;
    cout << endl;
    

    cout << endl << endl;
    
    cout << " --------------- Task 4B : Secant Convergence --------------- " << endl;
    // two initial guesses
    // parameters in the following order: kappa , vth , is
    vector<double> guess_0 = { 1.25 , 1.1 , 2.5e-7 };
    vector<double> guess_1 = { 1.0 , 1.0 , 1.0e-7 };
    // third guess using recurrence relation
    vector<double> guess_2;
    recurrenceRelation( &guess_2 , &guess_1 , &guess_0 ,VGS , VDS , IDS );
    //printMatrix( &guess_2 );

    // declare some variables for storing convergence answers
    int iterations; 
    vector<double> parameter_solutions; 
    double relative_residual , absolute_residual , least_squares; 

    // call the simple secant method for convergence
    /* input parameters are as follows:
     * iterations -> where the number of interations it took for the solution to converge will be stored
     * parameter_solutions -> to store the converged solutions
     * absolute , relative residuals -> to store the values of the residuals
     * least squares -> to store the least squares value after convergence
     * guess_0 -> oldest guess
     * guess_1 -> next guess
     * guess_2 -> third data point from the recurrence relation
     * VGS, VDS, IDS -> data from the file
     */
    secantConvergence( iterations , &parameter_solutions , 
        absolute_residual , relative_residual , least_squares ,
        &guess_0 , &guess_1 , &guess_2 , VGS , VDS , IDS );

    cout << endl;
    cout << " the converged solutions after " << iterations << " iterations are : " << endl;
    cout << " kappa = " << parameter_solutions[0] << endl;
    cout << " vth = " << parameter_solutions[1] << endl;
    cout << " is = " << parameter_solutions[2] << endl;
    cout << " absolute residual error = " << absolute_residual << endl;
    cout << " relative residual error = " << relative_residual << endl;
    cout << " least squares = " << least_squares << endl;

    // parameter sensitivity analysis
    /* input parameters are as follows:
     * parameter_solutions -> solution from convergence above
     * pertubation -> amount by which to disturb the parameter for analysis
     * which -> which parameter to perturb: 0 is kappa , 1 is vth and 2 is is
     * VGS, VDS, IDS -> data from the file
     */
    double kappa_sensitivity = parameterSensitivity( &parameter_solutions , 0.1 , 0 , VGS , VDS , IDS );
    double vth_sensitivity = parameterSensitivity ( &parameter_solutions , 0.1 , 1 , VGS , VDS , IDS );
    double is_sensitivity = parameterSensitivity ( &parameter_solutions , 0.1 , 2 , VGS , VDS , IDS );
    cout << " sensitivity with respect to kappa = " << kappa_sensitivity << endl;
    cout << " sensitivity with respect to vth   = " << vth_sensitivity << endl;
    cout << " sensitivity with respect to is    = " << is_sensitivity << endl;
    
    cout << endl << endl;
    cout << " --------------- Task 5B : Normalized Secant Convergence --------------- " << endl;
    // use same initial guesses as above
    //guess_2.erase(guess_2.begin(),guess_2.end());
    //recurrenceRelation( &guess_2 , &guess_1 , &guess_0 ,VGS , VDS , IDS , true );
    // reset some variables 
    iterations = 0;
    relative_residual = 1.0;
    absolute_residual = 1.0;
    least_squares =  1.0;
    parameter_solutions.erase(parameter_solutions.begin(),parameter_solutions.end());

    secantConvergence( iterations , &parameter_solutions , 
        absolute_residual , relative_residual , least_squares ,
        &guess_0 , &guess_1 , &guess_2 , VGS , VDS , IDS , true );

    cout << endl;
    cout << " the converged solutions after " << iterations << " iterations are : " << endl;
    cout << " kappa = " << parameter_solutions[0] << endl;
    cout << " vth = " << parameter_solutions[1] << endl;
    cout << " is = " << parameter_solutions[2] << endl;
    cout << " absolute residual error = " << absolute_residual << endl;
    cout << " relative residual error = " << relative_residual << endl;
    cout << " least squares = " << least_squares << endl;

    kappa_sensitivity = parameterSensitivity( &parameter_solutions , 0.1 , 0 , VGS , VDS , IDS );
    vth_sensitivity = parameterSensitivity ( &parameter_solutions , 0.1 , 1 , VGS , VDS , IDS );
    is_sensitivity = parameterSensitivity ( &parameter_solutions , 0.1 , 2 , VGS , VDS , IDS );
    cout << " sensitivity with respect to kappa = " << kappa_sensitivity << endl;
    cout << " sensitivity with respect to vth   = " << vth_sensitivity << endl;
    cout << " sensitivity with respect to is    = " << is_sensitivity << endl;

    cout << endl;

    return 0;
}
