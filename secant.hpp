//#include "utilityFunctions.hpp"
#include "quasiNewtonMethod.hpp"

void secantConvergence( int& iterations , vector<double>* parameter_solutions , 
    double& absolute_residual , double& relative_residual , double& least_squares ,
    vector<double>* guess_0 , vector<double>* guess_1 , vector<double>* guess_2 , 
    vector<double>* VGS , vector<double>* VDS , vector<double>* IDS ,
    bool normalize = false);
void recurrenceRelation( vector<double>* guess_2 , vector<double>* guess_1 , vector<double>* guess_0 ,
    vector<double>* VGS , vector<double>* VDS , vector<double>* IDS , bool normalize = false );
void secantGradient( vector<double>* gradients, vector<double>* kappa_history , vector<double>* vth_history , vector<double>* is_history , vector<double>* v_history );
void secantHessian( vector<vector<double>>* hessians , vector<double>* kappa_history , vector<double>* vth_history , vector<double>* is_history , vector<double>* v_history );
