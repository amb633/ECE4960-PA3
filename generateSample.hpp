//
//  generateSample.hpp
//  PA3
//
//  Created by Ariana Bruno on 4/7/18.
//  Copyright © 2018 Ariana Bruno. All rights reserved.
//

#ifndef generateSample_hpp
#define generateSample_hpp

#include <stdio.h>
#include <cmath>
#include <vector>
#include "fullSolver.hpp"

using namespace std;

void randomSamples(vector<double>* x, vector<double>* y, vector<double>* noise, vector<double>* y_n);
void LinearLSF(vector<double>* x, vector<double>* y, vector<vector<double>>* H, vector<double>* RHS);

#endif /* generateSample_hpp */
