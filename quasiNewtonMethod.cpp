//
//  quasiNewtonMethod.cpp
//  PA3
//
//  Created by Ariana Bruno on 4/7/18.
//  Copyright © 2018 Ariana Bruno. All rights reserved.
//

#include "quasiNewtonMethod.hpp"

void quasiNetwon_dx(vector<double>* VGS , vector<double>* VDS, vector<double>* IDS, double k, double Vth, double Is, vector<double>* delta_x){
    
    //initializes the Hessian matrix
    vector<vector<double>*>* H = new vector<vector<double>*>;
    zeroMatrix(H, 3);

    //Ids_model which will be filled in by modelIds for different paramter pertubations
    vector<double>* Ids_current = new vector<double>;
    vector<double>* Ids_2dk = new vector<double>;
    vector<double>* Ids_dk = new vector<double>;
    vector<double>* Ids_2dVth = new vector<double>;
    vector<double>* Ids_dVth = new vector<double>;
    vector<double>* Ids_2dIs = new vector<double>;
    vector<double>* Ids_dIs = new vector<double>;
    vector<double>* Ids_dkdVth = new vector<double>;
    vector<double>* Ids_dkdIs = new vector<double>;
    vector<double>* Ids_dVthdIs = new vector<double>;
    
    //calculating the pertubations of the parameters
    double pertubation = 1e-4;
    double dk = k*pertubation;
    double dVth = Vth*pertubation;
    double dIs = Is*pertubation;
    
    //calculating various valuse of Ids_model for different pertubations of the parameters
    //needed for the V gradient vector and Hessian matrix calculations
    modelIds(Ids_current, VGS, VDS, k, Vth, Is);
    modelIds(Ids_2dk, VGS, VDS, (k + 2*dk), Vth, Is);
    modelIds(Ids_dk, VGS, VDS, (k + dk), Vth, Is);
    modelIds(Ids_2dVth, VGS, VDS, k, (Vth + 2*dVth), Is);
    modelIds(Ids_dVth, VGS, VDS, k, (Vth + dVth), Is);
    modelIds(Ids_2dIs, VGS, VDS, k, Vth, (Is + 2*dIs));
    modelIds(Ids_dIs, VGS, VDS, k, Vth, (Is + dIs));
    modelIds(Ids_dkdVth, VGS, VDS, (k + dk), (Vth + dVth), Is);
    modelIds(Ids_dkdIs, VGS, VDS, (k + dk), Vth, (Is + dIs));
    modelIds(Ids_dVthdIs, VGS, VDS, k, (Vth + dVth), (Is + dIs));
    
    //calculating the current least squares result for the given parameters(unchanged)
    double current_sum_sqs = sumSquares(Ids_current, IDS);
    
    //calculating the V gradient vector terms
    vector<double>* v_grad = new vector<double>;
    (*v_grad).push_back((sumSquares(Ids_dk, IDS)-current_sum_sqs)/dk);
    (*v_grad).push_back((sumSquares(Ids_dVth, IDS)-current_sum_sqs)/dVth);
    (*v_grad).push_back((sumSquares(Ids_dIs, IDS)-current_sum_sqs)/dIs);
    
    //calculating the Hessian matrix terms
    (*(*H)[0])[0] = (sumSquares(Ids_2dk, IDS) - 2*((*v_grad)[0]) + current_sum_sqs)/(dk*dk);
    (*(*H)[0])[1] = (sumSquares(Ids_dkdVth, IDS) - (*v_grad)[0] - (*v_grad)[1] + current_sum_sqs)/(dk*dVth);
    (*(*H)[0])[2] = (sumSquares(Ids_dkdIs, IDS) - (*v_grad)[0] - (*v_grad)[2] + current_sum_sqs)/(dk*dIs);
    
    (*(*H)[1])[0] = (*(*H)[0])[1];
    (*(*H)[1])[1] = (sumSquares(Ids_2dVth, IDS) - 2*((*v_grad)[1]) + current_sum_sqs)/(dVth*dVth);
    (*(*H)[1])[2] = (sumSquares(Ids_dVthdIs, IDS) - (*v_grad)[1] - (*v_grad)[2] + current_sum_sqs)/(dVth*dIs);
    
    (*(*H)[2])[0] = (*(*H)[0])[2];
    (*(*H)[2])[1] = (*(*H)[1])[2];
    (*(*H)[2])[2] = (sumSquares(Ids_2dIs, IDS) - 2*((*v_grad)[2]) + current_sum_sqs)/(dIs*dIs);
    
    //solving for the delta of the parameters, given the Hessian matrix and V gradient vector
    fullSolver( delta_x , H , v_grad );
    for(int i = 0; i<(*delta_x).size(); i++){
        (*delta_x)[i] = (-1.0)* (*delta_x)[i];
    }

}

double t_adjusted_sum_sq( vector<double>* VGS, vector<double>* VDS, vector<double>* IDS, double t, vector<double>* delta, vector<double>* paramters){
    
    //calculating the least squares results for the scaled delta_paramters by t
    vector<double>* delta_adjusted = new vector<double>;
    scaleVector(t, delta, delta_adjusted);
    
    vector<double>* Ids_t_adjusted = new vector<double>;
    vector<double>* new_adjusted_paramters = new vector<double>;
    
    add_vectors(paramters, delta_adjusted, new_adjusted_paramters);
    modelIds(Ids_t_adjusted, VGS, VDS, (*new_adjusted_paramters)[0], (*new_adjusted_paramters)[1], (*new_adjusted_paramters)[2]);
    
    return sumSquares(Ids_t_adjusted, IDS);
}

double linear_search( vector<double>* VGS , vector<double>* VDS, vector<double>* IDS, vector<double>* paramters, vector<double>* delta, double t_min, double t_max){
    
    //finding the optimal t using binary search to get the best t for the smallest Least Squares result
    double t_mid = (t_min+t_max)/2.0;
    double t_mid_adj_sum_sqs = t_adjusted_sum_sq( VGS, VDS, IDS, t_mid, delta, paramters);
    double t_min_adj_sum_sqs = t_adjusted_sum_sq( VGS, VDS, IDS, t_min, delta, paramters);
    double t_max_adj_sum_sqs = t_adjusted_sum_sq( VGS, VDS, IDS, t_max, delta, paramters);
    
    double t_opt = t_mid;
    double sum_sq_opt = t_mid_adj_sum_sqs;
    
    if( isnan(t_max_adj_sum_sqs) || isnan(t_mid_adj_sum_sqs)){
        t_opt = linear_search( VGS , VDS, IDS, paramters, delta, t_min, t_mid);
    }
    if( t_min_adj_sum_sqs < sum_sq_opt ){
        t_opt = linear_search( VGS , VDS, IDS, paramters, delta, t_min, t_mid);
    }
    if( t_max_adj_sum_sqs < sum_sq_opt ){
        t_opt = linear_search( VGS , VDS, IDS, paramters, delta, t_mid, t_max);
    }
    return t_opt;
}

void quasiNetwon_itr( vector<double>* VGS , vector<double>* VDS, vector<double>* IDS, vector<double>* current_parameters,  vector<double>* new_parameters, double* norm_V, double* norm_delta_rel,  double* norm_delta_abs){

    double k = (*current_parameters)[0];
    double Vth = (*current_parameters)[1];
    double Is = (*current_parameters)[2];
    
    //calculating current sum squares for the parameters at the start of this iteration
    vector<double>* Ids_current = new vector<double>;
    modelIds(Ids_current, VGS, VDS, k, Vth, Is);
    double current_sum_sqs = sumSquares(Ids_current, IDS);
    
    //calculate the delta for the parameters in this iteration
    vector<double>* delta_parameters = new vector<double>;
    quasiNetwon_dx(VGS, VDS, IDS, k, Vth, Is, delta_parameters);
    
    //find the t to scale the delta for updating the parameters in this iterations
    double t = linear_search(VGS, VDS, IDS, current_parameters, delta_parameters);
    
    //Calculating the scaled detla for the parameters using t from the linear search
    vector<double>* t_delta_parameters = new vector<double> ;
    scaleVector(t, delta_parameters, t_delta_parameters);
    
    //updating the current parameters with t*delta to get the new parameters in this iteration
    add_vectors(current_parameters, t_delta_parameters, new_parameters);

    //calculating the absolute and relative residual error for the parameters
    *norm_delta_rel = delta_norm_rel(t_delta_parameters, new_parameters);
    *norm_delta_abs = delta_norm_abs(t_delta_parameters);
    
    vector<double>* Ids_new = new vector<double>;
    modelIds(Ids_new, VGS, VDS, (*new_parameters)[0], (*new_parameters)[1], (*new_parameters)[2]);
    double new_sum_sqs = sumSquares(Ids_new, IDS);
    
    //calculate the least squares for this iteration
    *norm_V = new_sum_sqs;
    
}

