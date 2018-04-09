//
//  fullSolver2.hpp
//  PA3
//
//  Created by Ariana Bruno on 4/7/18.
//  Copyright © 2018 Ariana Bruno. All rights reserved.
//

#ifndef fullSolver_hpp
#define fullSolver_hpp

#include <stdio.h>
#include <iostream>
#include <cstdlib>
#include <cmath>
#include <vector>

using namespace std;

void fullSolver( vector<double>* solution , vector<vector<double>>* matrix , vector<double>* b );
void conditionMatrix( vector<vector<double>>* matrix );
void printMatrix ( vector<vector<double>>* matrix );
void printMatrix( vector<double>* vector );
void identityMatrix( vector<vector<double>>* m , int rank );
void zeroMatrix( vector<vector<double>>* m , int rank );
void calcAtomicVector( vector<vector<double>>* m , vector<vector<double>>* matrix, int row );
void copyMatrix( vector<vector<double>>* copy , vector<vector<double>>* original );
void matrixProduct( vector<vector<double>>* result , vector<vector<double>>* matrix_1 , vector<vector<double>>* matrix_2 );
void vectorProduct( vector<double>* y , vector<vector<double>>* matrix , vector<double>* x );
void backwardSubstitution( vector<double>* x , vector<vector<double>>* U , vector<double>* y );


#endif /* fullSolver_hpp */
