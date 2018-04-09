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

vector<double> VGS;
vector<double> VDS;
vector<double> IDS;

//string path = "/Users/arianabruno/Desktop/ECE4960/ProgrammingAssignments/ECE4960-PA3/outputNMOS.txt";
string path = "C:/Users/Haritha/Documents/ECE4960-PAs/ECE4960-PA3/outputNMOS.txt";


int main(int argc, const char * argv[]) {
    cout << endl << " -------- Task 2 : Parameter Extraction for Power Law -------- " << endl << endl;
   
    readDataFile( path , &VGS ,  &VDS ,  &IDS );
    
    vector<double> x_samples;
    vector<double> y_samples;
    vector<double> noise_samples;
    vector<double> y_noisey_samples;
    
    randomSamples(&x_samples, &y_samples, &noise_samples, &y_noisey_samples);
    
    vector<vector<double>> H_matrix;
    vector<double> RHS;
    
//    find the H matrix and RHS vector for this problem for parameter extraction of  m and c0 for the S_measured data
    LinearLSF(&x_samples, &y_noisey_samples, &H_matrix, &RHS);
    vector<double> solution;

//    solve for the parameters given H and RHS
    fullSolver( &solution , &H_matrix , &RHS );
    
    cout << " This is the solution for m and log(c0)" << endl;
    printMatrix( &solution );
    cout << " m = " << (solution)[0] << endl;
    cout << " c0 = " << exp((solution)[1]) << endl;
    
    cout << endl << " -------- Task 4A : Quasi Newton Parameter Extraction-------- " << endl << endl;
    
    double k_0 = 1.0;
    double Vth_0 = 1.0;
    double Is_0 = 1e-7;
    
    vector<double> current_parameters;
    (current_parameters).push_back(k_0); (current_parameters).push_back(Vth_0); (current_parameters).push_back(Is_0);
    vector<double> updated_parameters;
    
    double norm_V, norm_delta_rel, norm_delta_abs;
    
    quasiNetwon_itr(&VGS, &VDS, &IDS, &current_parameters, &updated_parameters, &norm_V, &norm_delta_rel, &norm_delta_abs, false);
    (current_parameters) = (updated_parameters);
    
    int itr = 0;
    while( norm_V > 1e-9 && norm_delta_rel > 1e-9){
        vector<double> new_parameters;
        quasiNetwon_itr(&VGS, &VDS, &IDS, &current_parameters, &new_parameters, &norm_V, &norm_delta_rel, &norm_delta_abs, false);
        (current_parameters) = (new_parameters);
        itr++;
    }
    
    cout << " the converged solutions after " << itr << " iterations are: " << endl;
    cout << " kappa = " << (current_parameters)[0] << endl;
    cout << " V_th = " << (current_parameters)[1] << endl;
    cout << " Is = " << (current_parameters)[2] << endl;
    cout << " absolute residual error = " << norm_delta_abs << endl;
    cout << " relative residual error = " << norm_delta_rel << endl;
    cout << " least squares = " << norm_V << endl;
    
    double k_sensitivity = parameterSensitivity( &current_parameters , 0.1, 0 , &VGS, &VDS, &IDS );
    double Vth_sensitivity = parameterSensitivity ( &current_parameters , 0.1 , 1 , &VGS , &VDS , &IDS );
    double Is_sensitivity = parameterSensitivity ( &current_parameters , 0.1 , 2 , &VGS , &VDS , &IDS );
    cout << " sensitivity with respect to kappa = " << k_sensitivity << endl;
    cout << " sensitivity with respect to vth   = " << Vth_sensitivity << endl;
    cout << " sensitivity with respect to is    = " << Is_sensitivity << endl;
    
    cout << endl;
    
    cout << endl << " -------- Task 5A : Normalized Quasi Newton Parameter Extraction-------- " << endl << endl;
    
    cout << " for normalized soltions: " << endl << endl;
    
    vector<double> current_parameters_N;
    (current_parameters_N).push_back(k_0); (current_parameters_N).push_back(Vth_0); (current_parameters_N).push_back(Is_0);
    vector<double> updated_parameters_N;
    
    double norm_V_N, norm_delta_rel_N, norm_delta_abs_N;
    
    quasiNetwon_itr(&VGS, &VDS, &IDS, &current_parameters_N, &updated_parameters_N, &norm_V_N, &norm_delta_rel_N, &norm_delta_abs_N, true);
    (current_parameters_N) = (updated_parameters_N);
    
    int itr_N = 0;
    while( norm_V_N > 1e-9 && norm_delta_rel_N > 1e-9){
        vector<double> new_parameters_N;
        quasiNetwon_itr(&VGS, &VDS, &IDS, &current_parameters_N, &new_parameters_N, &norm_V_N, &norm_delta_rel_N, &norm_delta_abs_N, true);
        (current_parameters_N) = (new_parameters_N);
        itr_N++;
    }
    
    cout << " the converged solutions after " << itr_N << " iterations are: " << endl;
    cout << " kappa = " << (current_parameters_N)[0] << endl;
    cout << " V_th = " << (current_parameters_N)[1] << endl;
    cout << " Is = " << (current_parameters_N)[2] << endl;
    cout << " absolute residual error = " << norm_delta_abs_N << endl;
    cout << " relative residual error = " << norm_delta_rel_N << endl;
    cout << " least squares = " << norm_V_N << endl;
    cout << endl;
    double k_sensitivity_N = parameterSensitivity( &current_parameters_N , 0.1 , 0 , &VGS , &VDS , &IDS );
    double Vth_sensitivity_N = parameterSensitivity ( &current_parameters_N , 0.1 , 1 , &VGS , &VDS , &IDS );
    double Is_sensitivity_N = parameterSensitivity ( &current_parameters_N , 0.1 , 2 , &VGS , &VDS , &IDS );
    cout << " sensitivity with respect to kappa = " << k_sensitivity_N << endl;
    cout << " sensitivity with respect to vth   = " << Vth_sensitivity_N << endl;
    cout << " sensitivity with respect to is    = " << Is_sensitivity_N << endl;
    

    cout << endl << endl;
    
    cout << " --------------- Task 4B : Secant Convergence --------------- " << endl;
    // two initial guesses
    // parameters in the following order: kappa , vth , is
    vector<double> guess_0 = { 1.25 , 1.1 , 2.5e-7 };
    vector<double> guess_1 = { 1.0 , 1.0 , 1.0e-7 };
    // third guess using recurrence relation
    vector<double> guess_2;
    recurrenceRelation( &guess_2 , &guess_1 , &guess_0 , &VGS , &VDS , &IDS );
    //printMatrix( &guess_2 );

    // declare some variables for storing convergence answers
    int iterations; 
    vector<double> parameter_solutions; 
    double relative_residual , absolute_residual , least_squares; 

    // call the simple secant method for convergence
    //* input parameters are as follows:
    //* iterations -> where the number of interations it took for the solution to converge will be stored
    //* parameter_solutions -> to store the converged solutions
    //* absolute , relative residuals -> to store the values of the residuals
    //* least squares -> to store the least squares value after convergence
    //* guess_0 -> oldest guess
    //* guess_1 -> next guess
    //* guess_2 -> third data point from the recurrence relation
    //* VGS, VDS, IDS -> data from the file
    
    secantConvergence( iterations , &parameter_solutions ,
        absolute_residual , relative_residual , least_squares ,
        &guess_0 , &guess_1 , &guess_2 , &VGS , &VDS , &IDS );

    cout << endl;
    cout << " the converged solutions after " << iterations << " iterations are : " << endl;
    cout << " kappa = " << parameter_solutions[0] << endl;
    cout << " vth = " << parameter_solutions[1] << endl;
    cout << " is = " << parameter_solutions[2] << endl;
    cout << " absolute residual error = " << absolute_residual << endl;
    cout << " relative residual error = " << relative_residual << endl;
    cout << " least squares = " << least_squares << endl;

    // parameter sensitivity analysis
    //* input parameters are as follows:
    //* parameter_solutions -> solution from convergence above
    //* pertubation -> amount by which to disturb the parameter for analysis
    //* which -> which parameter to perturb: 0 is kappa , 1 is vth and 2 is is
    //* VGS, VDS, IDS -> data from the file
     
    double kappa_sensitivity = parameterSensitivity( &parameter_solutions , 0.1 , 0 , &VGS , &VDS , &IDS );
    double vth_sensitivity = parameterSensitivity ( &parameter_solutions , 0.1 , 1 , &VGS , &VDS , &IDS );
    double is_sensitivity = parameterSensitivity ( &parameter_solutions , 0.1 , 2 , &VGS , &VDS , &IDS );
    cout << " sensitivity with respect to kappa = " << kappa_sensitivity << endl;
    cout << " sensitivity with respect to vth   = " << vth_sensitivity << endl;
    cout << " sensitivity with respect to is    = " << is_sensitivity << endl;
    cout << endl << endl;

    cout << " --------------- Task 5B : Normalized Secant Convergence --------------- " << endl;
    // use same initial guesses as above
    // reset some variables 
    iterations = 0;
    relative_residual = 1.0;
    absolute_residual = 1.0;
    least_squares =  1.0;
    parameter_solutions.erase(parameter_solutions.begin(),parameter_solutions.end());

    secantConvergence( iterations , &parameter_solutions , 
        absolute_residual , relative_residual , least_squares ,
        &guess_0 , &guess_1 , &guess_2 , &VGS , &VDS , &IDS , true );

    cout << endl;
    cout << " the converged solutions after " << iterations << " iterations are : " << endl;
    cout << " kappa = " << parameter_solutions[0] << endl;
    cout << " vth = " << parameter_solutions[1] << endl;
    cout << " is = " << parameter_solutions[2] << endl;
    cout << " absolute residual error = " << absolute_residual << endl;
    cout << " relative residual error = " << relative_residual << endl;
    cout << " least squares = " << least_squares << endl;

    kappa_sensitivity = parameterSensitivity( &parameter_solutions , 0.1 , 0 , &VGS , &VDS , &IDS );
    vth_sensitivity = parameterSensitivity ( &parameter_solutions , 0.1 , 1 , &VGS , &VDS , &IDS );
    is_sensitivity = parameterSensitivity ( &parameter_solutions , 0.1 , 2 , &VGS , &VDS , &IDS );
    cout << " sensitivity with respect to kappa = " << kappa_sensitivity << endl;
    cout << " sensitivity with respect to vth   = " << vth_sensitivity << endl;
    cout << " sensitivity with respect to is    = " << is_sensitivity << endl;
    cout << endl << endl;
    
    cout << " --------------- Task 6A : Unnormalized Quasi Newton Convergence with Different Starting Points --------------- " << endl;
    
    vector<double> kappa_initial_points = { 0.5 , 0.6 , 0.7 , 0.8 , 0.9 , 1.0 };
    vector<double> vth_initial_points = { 0.8 , 0.9 , 1.0 , 1.1 , 1.2 , 1.3 , 1.4 , 1.5 , 1.6 , 1.7 , 1.8 , 1.9 , 2.0 };
    vector<double> is_initial_points = { 1e-8 , 3e-8 , 1e-7 , 3e-7 , 1e-6 , 3e-6 , 1e-5 , 3e-5 };
    
    vector<double> current_parameters_itr;
    vector<double> updated_parameters_itr;

    double min_least_squares_uQuasi = 10000.0;
    double min_absolute_residual_uQuasi = 0.0;
    double min_relative_residual_uQuasi = 0.0;
    double min_kappa_uQuasi = 0.0;
    double min_vth_uQuasi = 0.0;
    double min_is_uQuasi = 0.0;
    double min_kappa_initial_uQuasi = 0.0;
    double min_vth_initial_uQuasi = 0.0;
    double min_is_initial_uQuasi = 0.0;
    double min_kappa_sensitivity_uQuasi = 0.0;
    double min_vth_sensitivity_uQuasi = 0.0;
    double min_is_sensitivity_uQuasi = 0.0;
    
    for ( int itk = 0 ; itk < kappa_initial_points.size() ; itk++ ){
        for ( int itv = 0 ; itv < vth_initial_points.size() ; itv++ ){
            for ( int its = 0 ; its < is_initial_points.size() ; its++ ){
                k_0 = kappa_initial_points[itk];
                Vth_0 = vth_initial_points[itv];
                Is_0 = is_initial_points[its];
                
                vector<double> current_parameters_itr;
                vector<double> updated_parameters_itr;
                (current_parameters_itr).push_back(k_0); (current_parameters_itr).push_back(Vth_0);(current_parameters_itr).push_back(Is_0);
                
                double norm_V, norm_delta_rel, norm_delta_abs;
                
                quasiNetwon_itr(&VGS, &VDS, &IDS, &current_parameters_itr, &updated_parameters_itr, &norm_V, &norm_delta_rel, &norm_delta_abs, false);
                (current_parameters_itr) = (updated_parameters_itr);
                
                int itr = 0;
                while( norm_V > 1e-9 && norm_delta_rel > 1e-9 && itr < 100 ){
                    vector<double> new_parameters_itr;
                    quasiNetwon_itr(&VGS, &VDS, &IDS, &current_parameters_itr, &new_parameters_itr, &norm_V, &norm_delta_rel, &norm_delta_abs, false);
                    (current_parameters_itr) = (new_parameters_itr);
                    itr++;
                    if( itr >= 100 ){
                        cout << "did not converge" << endl;
                    }
                }
                
                cout << " the converged solutions after " << itr << " iterations are: " << endl;
                cout << " kappa = " << (current_parameters_itr)[0] << endl;
                cout << " V_th = " << (current_parameters_itr)[1] << endl;
                cout << " Is = " << (current_parameters_itr)[2] << endl;
                cout << " absolute residual error = " << norm_delta_abs << endl;
                cout << " relative residual error = " << norm_delta_rel << endl;
                cout << " least squares = " << norm_V << endl;
                
                double k_sensitivity = parameterSensitivity( &current_parameters_itr , 0.1 , 0 , &VGS , &VDS , &IDS );
                double Vth_sensitivity = parameterSensitivity ( &current_parameters_itr , 0.1 , 1 , &VGS , &VDS , &IDS );
                double Is_sensitivity = parameterSensitivity ( &current_parameters_itr , 0.1 , 2 , &VGS , &VDS , &IDS );
                cout << " sensitivity with respect to kappa = " << k_sensitivity << endl;
                cout << " sensitivity with respect to vth   = " << Vth_sensitivity << endl;
                cout << " sensitivity with respect to is    = " << Is_sensitivity << endl;
                
                if ( norm_V < min_least_squares_uQuasi ){
                	min_least_squares_uQuasi = norm_V;
                	min_absolute_residual_uQuasi = norm_delta_abs;
                	min_relative_residual_uQuasi = norm_delta_rel;
                	min_kappa_uQuasi = (current_parameters_itr)[0];
                	min_vth_uQuasi = (current_parameters_itr)[1];
                	min_is_uQuasi = (current_parameters_itr)[2];
                	min_kappa_initial_uQuasi = k_0;
                	min_vth_initial_uQuasi = Vth_0;
                	min_is_initial_uQuasi = Is_0;
                	min_kappa_sensitivity_uQuasi = k_sensitivity;
                	min_vth_sensitivity_uQuasi = Vth_sensitivity;
                	min_is_sensitivity_uQuasi = Is_sensitivity;
                }

                cout << endl;
                
            }
        }
    }
    
    cout << " --------------- Task 6B : Unnormalized Secant Convergence with Different Starting Points --------------- " << endl;
    
    double min_least_squares_uSecant = 10000.0;
    double min_absolute_residual_uSecant = 0.0;
    double min_relative_residual_uSecant = 0.0;
    double min_kappa_uSecant = 0.0;
    double min_vth_uSecant = 0.0;
    double min_is_uSecant = 0.0;
    double min_kappa_initial_uSecant = 0.0;
    double min_vth_initial_uSecant = 0.0;
    double min_is_initial_uSecant = 0.0;
    double min_kappa_sensitivity_uSecant = 0.0;
    double min_vth_sensitivity_uSecant = 0.0;
    double min_is_sensitivity_uSecant = 0.0;

    for ( int itk = 0 ; itk < kappa_initial_points.size() ; itk++ ){
        for ( int itv = 0 ; itv < vth_initial_points.size() ; itv++ ){
            for ( int its = 0 ; its < is_initial_points.size() ; its++ ){
                // clear some old variables
                guess_0.erase(guess_0.begin(),guess_0.end());
                guess_1.erase(guess_1.begin(),guess_1.end());
                guess_2.erase(guess_2.begin(),guess_2.end());
                parameter_solutions.erase(parameter_solutions.begin(),parameter_solutions.end());
                iterations = 0;
                relative_residual = 1.0;
                absolute_residual = 1.0;
                least_squares = -1.0;
                kappa_sensitivity = 0.0;
                vth_sensitivity = 0.0;
                is_sensitivity = 0.0;
                
                // inital guess
                guess_1 = { kappa_initial_points[itk] , vth_initial_points[itv] , is_initial_points[its] };
                guess_0 = { 1.0 , 1.0 , 1.0e-7 };

                recurrenceRelation( &guess_2 , &guess_1 , & guess_0 , &VGS , &VDS , &IDS );
                secantConvergence( iterations , &parameter_solutions , 
                    absolute_residual , relative_residual , least_squares ,
                    &guess_0 , &guess_1 , &guess_2 , &VGS , &VDS , &IDS );
                cout << endl;
                cout << " parameters being tested : " ; printMatrix(&guess_1);
                cout << " the converged solutions after " << iterations << " iterations are : " << endl;
                cout << " kappa = " << parameter_solutions[0] << endl;
                cout << " vth = " << parameter_solutions[1] << endl;
                cout << " is = " << parameter_solutions[2] << endl;
                cout << " absolute residual error = " << absolute_residual << endl;
                cout << " relative residual error = " << relative_residual << endl;
                cout << " least squares = " << least_squares << endl;

                kappa_sensitivity = parameterSensitivity( &parameter_solutions , 0.1 , 0 , &VGS , &VDS , &IDS );
                vth_sensitivity = parameterSensitivity ( &parameter_solutions , 0.1 , 1 , &VGS , &VDS , &IDS );
                is_sensitivity = parameterSensitivity ( &parameter_solutions , 0.1 , 2 , &VGS , &VDS , &IDS );
                cout << " sensitivity with respect to kappa = " << kappa_sensitivity << endl;
                cout << " sensitivity with respect to vth   = " << vth_sensitivity << endl;
                cout << " sensitivity with respect to is    = " << is_sensitivity << endl;
                cout << endl << endl;

                if ( least_squares < min_least_squares_uSecant ){
                	min_least_squares_uSecant = least_squares;
                	min_absolute_residual_uSecant = absolute_residual;
                	min_relative_residual_uSecant = relative_residual;
                	min_kappa_uSecant = parameter_solutions[0];
                	min_vth_uSecant = parameter_solutions[1];
                	min_is_uSecant = parameter_solutions[2];
                	min_kappa_initial_uSecant = kappa_initial_points[itk];
                	min_vth_initial_uSecant = vth_initial_points[itv];
                	min_is_initial_uSecant = is_initial_points[its];
                	min_kappa_sensitivity_uSecant = kappa_sensitivity;
                	min_vth_sensitivity_uSecant = vth_sensitivity;
                	min_is_sensitivity_uSecant = is_sensitivity;
                }


            }
        }
    }

    cout << " --------------- Task 6D : Normalized Secant Convergence with Different Starting Points --------------- " << endl;
    
    double min_least_squares_Secant = 10000.0;
    double min_absolute_residual_Secant = 0.0;
    double min_relative_residual_Secant = 0.0;
    double min_kappa_Secant = 0.0;
    double min_vth_Secant = 0.0;
    double min_is_Secant = 0.0;
    double min_kappa_initial_Secant = 0.0;
    double min_vth_initial_Secant = 0.0;
    double min_is_initial_Secant = 0.0;
    double min_kappa_sensitivity_Secant = 0.0;
    double min_vth_sensitivity_Secant = 0.0;
    double min_is_sensitivity_Secant = 0.0;

    for ( int itk = 0 ; itk < kappa_initial_points.size() ; itk++ ){
        for ( int itv = 0 ; itv < vth_initial_points.size() ; itv++ ){
            for ( int its = 0 ; its < is_initial_points.size() ; its++ ){
                // clear some old variables
                guess_0.erase(guess_0.begin(),guess_0.end());
                guess_1.erase(guess_1.begin(),guess_1.end());
                guess_2.erase(guess_2.begin(),guess_2.end());
                parameter_solutions.erase(parameter_solutions.begin(),parameter_solutions.end());
                iterations = 0;
                relative_residual = 1.0;
                absolute_residual = 1.0;
                least_squares = -1.0;
                kappa_sensitivity = 0.0;
                vth_sensitivity = 0.0;
                is_sensitivity = 0.0;
                
                // inital guess
                guess_1 = { kappa_initial_points[itk] , vth_initial_points[itv] , is_initial_points[its] };
                guess_0 = { 1.0 , 1.0 , 1.0e-7 };

                recurrenceRelation( &guess_2 , &guess_1 , & guess_0 , &VGS , &VDS , &IDS );
                secantConvergence( iterations , &parameter_solutions , 
                    absolute_residual , relative_residual , least_squares ,
                    &guess_0 , &guess_1 , &guess_2 , &VGS , &VDS , &IDS );
                cout << endl;
                cout << " parameters being tested : " ; printMatrix(&guess_1);
                cout << " the converged solutions after " << iterations << " iterations are : " << endl;
                cout << " kappa = " << parameter_solutions[0] << endl;
                cout << " vth = " << parameter_solutions[1] << endl;
                cout << " is = " << parameter_solutions[2] << endl;
                cout << " absolute residual error = " << absolute_residual << endl;
                cout << " relative residual error = " << relative_residual << endl;
                cout << " least squares = " << least_squares << endl;

                kappa_sensitivity = parameterSensitivity( &parameter_solutions , 0.1 , 0 , &VGS , &VDS , &IDS );
                vth_sensitivity = parameterSensitivity ( &parameter_solutions , 0.1 , 1 , &VGS , &VDS , &IDS );
                is_sensitivity = parameterSensitivity ( &parameter_solutions , 0.1 , 2 , &VGS , &VDS , &IDS );
                cout << " sensitivity with respect to kappa = " << kappa_sensitivity << endl;
                cout << " sensitivity with respect to vth   = " << vth_sensitivity << endl;
                cout << " sensitivity with respect to is    = " << is_sensitivity << endl;
                cout << endl << endl;

                if ( least_squares < min_least_squares_Secant ){
                	min_least_squares_Secant = least_squares;
                	min_absolute_residual_Secant = absolute_residual;
                	min_relative_residual_Secant = relative_residual;
                	min_kappa_Secant = parameter_solutions[0];
                	min_vth_Secant = parameter_solutions[1];
                	min_is_Secant = parameter_solutions[2];
                	min_kappa_initial_Secant = kappa_initial_points[itk];
                	min_vth_initial_Secant = vth_initial_points[itv];
                	min_is_initial_Secant = is_initial_points[its];
                	min_kappa_sensitivity_Secant = kappa_sensitivity;
                	min_vth_sensitivity_Secant = vth_sensitivity;
                	min_is_sensitivity_Secant = is_sensitivity;
                }


            }
        }
    }

    cout << endl;
    return 0;
}
