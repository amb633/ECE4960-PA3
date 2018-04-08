//
//  quasiNewtonMethod.hpp
//  PA3
//
//  Created by Ariana Bruno on 4/7/18.
//  Copyright © 2018 Ariana Bruno. All rights reserved.
//

#ifndef quasiNewtonMethod_hpp
#define quasiNewtonMethod_hpp

#include <stdio.h>
#include <vector>
#include "fullSolver.hpp"
#include "utilityFunctions.hpp"

using namespace std;

void quasiNetwon_dx(vector<double>* VGS , vector<double>* VDS, vector<double>* IDS, double k, double Vth, double Is, vector<double>* delta_x);
double t_adjusted_sum_sq( vector<double>* VGS, vector<double>* VDS, vector<double>* IDS, double t, vector<double>* delta, vector<double>* paramters);
double linear_search( vector<double>* VGS , vector<double>* VDS, vector<double>* IDS, vector<double>* paramters, vector<double>* delta, double t_min = 0.0, double t_max = 1.0 );
void quasiNetwon_itr( vector<double>* VGS , vector<double>* VDS, vector<double>* IDS, vector<double>* current_paramters,  vector<double>* new_parameters, double* norm_V, double* norm_delta_rel,  double* norm_delta_abs);
#endif /* quasiNewtonMethod_hpp */
