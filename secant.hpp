#include "utilityFunctions.hpp"

void secantConvergence( int& iterations , vector<double>* parameter_solutions , double& relative_residual ,
	vector<double>* guess_0 , vector<double>* guess_1 , vector<double>* guess_2 , 
	vector<double>* VGS , vector<double>* VDS , vector<double>* IDS );
void secantGradient( vector<double>* gradients, vector<double>* kappa_history , vector<double>* vth_history , vector<double>* is_history , vector<double>* v_history );
void secantHessian( vector<vector<double>*>* hessians , vector<double>* kappa_history , vector<double>* vth_history , vector<double>* is_history , vector<double>* v_history );