#include "secant.hpp"

void secantGradient( vector<double>* gradients, vector<double>* kappa_history , vector<double>* vth_history , vector<double>* is_history , vector<double>* v_history ){
	// variable_history is a vector of the form : variable(k-2) , variable(k-1) , variable(k)

	// extract all the current and previous values
	double v_0 = (*v_history)[2];
	double v_1 = (*v_history)[1];

	double kappa_0 = (*kappa_history)[2];
	double kappa_1 = (*kappa_history)[1];

	double vth_0 = (*vth_history)[2];
	double vth_1 = (*vth_history)[1]; 

	double is_0 = (*is_history)[2];
	double is_1 = (*is_history)[1];

	// calculate the gradients using the secant method
	double gradientKappa = ( v_0 - v_1 )/( kappa_0 - kappa_1 );
	double gradientVth = ( v_0 - v_1 )/( vth_0 - vth_1 );
	double gradientIs = ( v_0 - v_1 )/( is_0 - is_1 );

	(*gradients).push_back(gradientKappa);
	(*gradients).push_back(gradientVth);
	(*gradients).push_back(gradientIs);
	return;
};

void secantHessian( vector<vector<double>*>* hessians , vector<double>* kappa_history , vector<double>* vth_history , vector<double>* is_history , vector<double>* v_history ) {

	double kappa_0 = (*kappa_history)[2];
	double kappa_1 = (*kappa_history)[1]; 
	double kappa_2 = (*kappa_history)[0];

	double vth_0 = (*vth_history)[2];
	double vth_1 = (*vth_history)[1]; 
	double vth_2 = (*vth_history)[0]; 

	double is_0 = (*is_history)[2];
	double is_1 = (*is_history)[1]; 
	double is_2 = (*is_history)[0];

	double v_0 = (*v_history)[2];
	double v_1 = (*v_history)[1]; 
	double v_2 = (*v_history)[0];
	
	/*double h_kappa_kappa = (((v_0 - v_1)/(kappa_0 - kappa_1)) - (( v_1 - v_2 )/(kappa_1 - kappa_2)))/(kappa_0 - kappa_1);
	double h_kappa_vth   = (((v_0 - v_1)/(kappa_0 - kappa_1)) - (( v_1 - v_2 )/(kappa_1 - kappa_2)))/(vth_0 - vth_1);
	double h_kappa_is   = (((v_0 - v_1)/(kappa_0 - kappa_1)) - (( v_1 - v_2 )/(kappa_1 - kappa_2)))/(is_0 - is_1);

	double h_vth_kappa = ((( v_0 - v_1)/(vth_0 - vth_1)) - (( v_1 - v_2 )/(vth_1 - vth_2)))/(kappa_0 - kappa_1);
	double h_vth_vth   = ((( v_0 - v_1)/(vth_0 - vth_1)) - (( v_1 - v_2 )/(vth_1 - vth_2)))/(vth_0 - vth_1);
	double h_vth_is    = ((( v_0 - v_1)/(vth_0 - vth_1)) - (( v_1 - v_2 )/(vth_1 - vth_2)))/(is_0 - is_1);

	double h_is_kappa = ((( v_0 - v_1)/(is_0 - is_1)) - (( v_1 - v_2 )/(is_1 - is_2)))/(kappa_0 - kappa_1);
	double h_is_vth   = ((( v_0 - v_1)/(is_0 - is_1)) - (( v_1 - v_2 )/(is_1 - is_2)))/(vth_0 - vth_1);
	double h_is_is    = ((( v_0 - v_1)/(is_0 - is_1)) - (( v_1 - v_2 )/(is_1 - is_2)))/(is_0 - is_1);

	zeroMatrix( hessians , 3 );
	(*(*hessians)[0])[0] = h_kappa_kappa;
	(*(*hessians)[1])[1] = h_vth_vth;
	(*(*hessians)[2])[2] = h_is_is;

	(*(*hessians)[0])[1] = h_kappa_vth;
	(*(*hessians)[0])[2] = h_kappa_is;

	(*(*hessians)[1])[0] = h_vth_kappa;
	(*(*hessians)[1])[1] = h_vth_is;

	(*(*hessians)[2])[0] = h_is_kappa;
	(*(*hessians)[2])[1] = h_is_vth;*/

	double h_kappa = (v_0 - v_1)/ ((kappa_0 - kappa_1)*(kappa_0 - kappa_1));
	double h_vth = (v_0 - v_1)/ ((vth_0 - vth_1)*(vth_0 - vth_1));
	double h_is = (v_0 - v_1)/ ((is_0 - is_1)*(is_0 - is_1));

	zeroMatrix( hessians , 3 );
	(*(*hessians)[0])[0] = h_kappa;
	(*(*hessians)[1])[1] = h_vth;
	(*(*hessians)[2])[2] = h_is;	

}