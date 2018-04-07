#include <iostream>
#include <cstdlib>
#include <cmath>
#include <experimental/random>

using namespace std;

double binarySearch( double start , double end );
double newtonRhapson( double guess );
double exp_fcn( double x );
double exp_der( double x );
double newtonWithLine( double guess );
double quasiNewton( double guess , double pert );
double secantMethod( double guess_1 , double guess_2 , double& f_val_1 );

int quasiNewtonOptimization( double& delta_1 , double& delta_2 , double x1 , double x2 , double pert );
int steepDescent( double& delta_1 , double& delta_2 , double x1 , double x2 , double pert );
double fcn_v( double x1 , double x2 );
int fcn_del_V( double& j1 , double& j2 , double x1 , double x2 , double pert );
int fcn_hessian( double& h11 , double& h22 , double& h12 , double& h21 , double x1 , double x2 , double pert );
int main (){

	cout.precision(6);

	cout << endl << " ---------- Hacker Practice 5.1 : Bisection Search Method ---------- " << endl;
	double ans = binarySearch( -5 , 10 );
	cout << " solution to binary search = " << ans << endl;

	cout << endl << " ---------- Hacker Practice 5.2 : Newton-Rhapson Method ---------- " << endl;
	double residual = 1.0;
	double iterate = 1.0;
	double delta = 0.0;
	int counter = 0;
	while( residual > 10e-10 ){
		iterate += delta;
		delta = newtonRhapson( iterate );
		residual = abs( delta/iterate );
		counter++;
	}
	cout << " solution to Newton's Method with inital guess of 1.0 = " << iterate << ", found after " << counter << " iterations" << endl;

	residual = 1.0;
	iterate = 10.0;
	delta = 0.0;
	counter = 0;
	while( residual > 10e-10 ){
		iterate += delta;
		delta = newtonRhapson( iterate );
		residual = abs( delta/iterate );
		counter++;
	}
	cout << " solution to Newton's Method with inital guess of 10.0 = " << iterate << ", found after " << counter << " iterations" << endl;

	cout << endl << " ---------- Hacker Practice 5.3 : Newton-Rhapson Method with Line Search ---------- " << endl;
	residual = 1.0;
	iterate = 10.0;
	delta = 0.0;
	counter = 0;
	while( residual > 10e-10 ){
		iterate += delta;
		delta = newtonWithLine( iterate );
		residual = abs( delta/iterate );
		counter++;
	}
	cout << " solution to Newton with Line Search and inital guess of 10.0 = " << iterate << ", found after " << counter << " iterations" << endl;

	cout << endl << " ---------- Hacker Practice 5.4 : Quasi-Newton Method ---------- " << endl;
	residual = 1.0;
	iterate = 1.0;
	delta = 0.0;
	counter = 0;
	double delta_prev = 10.0;
	double pert = 1e-4;
	while( residual > 10e-10 ){
		iterate += delta;
		delta = quasiNewton( iterate , pert );
		residual = abs( delta/iterate );
		counter++;
	}
	cout << " solution to quasi newton method = " << iterate << " , found after " << counter << " iterations" << endl;

	residual = 1.0;
	double iterate_1 = 1.1;
	double iterate_2 = 1.0;
	double f_val_1 = exp_fcn( iterate_1 );
	delta = 0.0;
	counter = 0;
	while( residual > 10e-10 ){
		iterate_2 += delta;
		delta = secantMethod( iterate_1 , iterate_2 , f_val_1 );
		//cout << delta << endl;
		residual = abs( delta/iterate_2 );
		iterate_1 = iterate_2;
		counter++;
	}
	cout << " solution to secant method = " << iterate_2 << " , found after " << counter << " iterations" << endl;

	cout << endl << " ---------- Hacker Practice 5.5 : Optimization using Newton Method ---------- " << endl;
	double x1 = 1.1;
	double x2 = 1.1;
	double delta_1 = 0.0;
	double delta_2 = 0.0;
	double residual_1 = 1.0;
	double residual_2 = 1.0;
	pert = 1e-4;
	counter = 0;
	while ( residual_1 > 1e-9 || residual_2 > 1e-9 ) {
	//for ( int i = 0 ; i < 100 ; i ++ ){
		x1 += delta_1;
		x2 += delta_2;
		quasiNewtonOptimization( delta_1 , delta_2 , x1 , x2 , pert );
		residual_1 = abs( delta_1/x1 );
		residual_2 = abs( delta_2/x2 );
		counter++;
	}
	cout << " solution to newton optimization method = " << x1 << " , " << x2 << " , found after " << counter << " iterations" << endl;
	
	cout << endl << " ---------- Hacker Practice 5.6 : Steepest Descent Search ---------- " << endl;
	x1 = 0.5;
	x2 = 0.5;
	delta_1 = -1.0;
	delta_2 = -1.0;
	residual_1 = 1.0;
	residual_2 = 1.0;
	pert = 1e-4;
	counter = 0;
	//while ( residual_1 > 1e-9 || residual_2 > 1e-9 ){
		x1 += delta_1;
		x2 += delta_2;
		steepDescent( delta_1 , delta_2 , x1 , x2 , pert );
		residual_1 = abs( delta_1/x1 );
		residual_2 = abs( delta_2/x2 );
		counter++;
	//}	

	cout << endl << " ---------- Hacker Practice 5.7 : Parameter Extraction ---------- " << endl;

	cout << endl << endl;
	return 0;
}

double binarySearch( double start , double end ){
	double midpoint = (end + start)/2.0;
	double test = exp(midpoint);
	if ( abs(test - 1.0 ) == 0.0 ) return midpoint;
	else if ( (test - 1.0) > 0 ) return binarySearch( start , midpoint );
	else if ( (test - 1.0) < 0 ) return binarySearch( midpoint , end );
}

double newtonRhapson( double guess ){
	double f_val = exp_fcn( guess );
	double f_der = exp_der( guess );
	double delta = -f_val/f_der;
	return delta;
}

double newtonWithLine( double guess ){
	double f_val = exp_fcn( guess );
	double f_der = exp_der( guess );
	double delta = -f_val/f_der;
	double test = exp_fcn( guess + delta );
	while( test < 0.0 ){
		delta /= 2.0;
		test = exp_fcn( guess + delta );
	}
	return delta;
}

double quasiNewton( double guess , double pert ){
	double f_val = exp_fcn( guess );
	double f_der = ( exp_fcn((1.0+pert)*guess) - exp_fcn (guess) ) / (pert*guess);
	double delta = -f_val/f_der;
	if ( isinf(abs(delta))) return 0.0;
	double test = -1.0;
	return delta;
}

double secantMethod( double guess_1 , double guess_2 , double& f_val_1 ){
	double f_val_2 = exp_fcn( guess_2 );
	double f_der = ( f_val_2 - f_val_1 ) / ( guess_2 - guess_1 );
	double delta = -f_val_1/f_der;
	if( isinf(abs(delta))) return 0.0;
	f_val_1 = f_val_2;
	double test = exp_fcn( guess_2 + delta );
	while( test < 0.0 ){
		delta /= 2.0;
		test = exp_fcn( guess_2 + delta );
	}
	return delta;
}

int quasiNewtonOptimization( double& delta_1 , double& delta_2 , double x1 , double x2 , double pert ){
	double j1 , j2 , h11 , h12 , h21 , h22;
	fcn_del_V( j1 , j2 , x1 , x2 , pert );
	fcn_hessian( h11 , h22 , h12 , h21 , x1 , x2 , pert );
	// invert the hessian
	double det = h11*h22 - h12*h21;
	double temp = h11;
	h11 = h22/det;
	h22 = temp/det;
	h12 = -h12/det;
	h21 = -h21/det;

	delta_1 = -(h11*j1 + h12*j2);
	delta_2 = -(h21*j1 + h22*j2);

	//cout << j1 << " , " << j2 << " , " << h11 << " , " << h12 << " , " << h21 << " , " << h22 << endl;
	return 0;

}

int steepDescent( double& delta_1 , double& delta_2 , double x1 , double x2 , double pert ){

}

double exp_fcn( double x ){
	return exp( 50.0*x ) - 1.0;
}

double exp_der( double x ){
	return 50.0*exp( 50.0*x );
}

double fcn_v ( double x1 , double x2 ){
	double v1 = 3.0*pow(x1,2.0) + x2 - 4.0;
	double v2 = pow(x1,2.0) -3.0*x2 + 2.0;
	return ( pow(v1,2) + pow(v2,2) );
}

int fcn_del_V( double& j1 , double& j2 , double x1 , double x2 , double pert ){
	j1 = 40.0*pow(x1,3.0) - 40.0*x1; //
	j1 = ( fcn_v( x1 + pert , x2 ) - fcn_v( x1 , x2 ) ) / ( pert );
	j2 = 20.0*x2 - 20.0; //
	j2 = ( fcn_v( x1 , x2 + pert ) - fcn_v( x1 , x2 ) ) / ( pert );
	return 0;
}

int fcn_hessian( double& h11 , double& h22 , double& h12 , double& h21 , double x1 , double x2 , double pert ){
	h11 = ( fcn_v( (x1+2*pert) , x2 ) - 2*fcn_v( (x1+pert) , x2 ) + fcn_v( x1 , x2 ) ) / ( pert* pert );
	h22 = ( fcn_v( x1 , (x2+2*pert) ) - 2*fcn_v( x1 , (x2+pert) ) + fcn_v( x1 , x2 ) ) / ( pert* pert );
	h12 = ( fcn_v( (x1+pert) , (x2+pert) ) - fcn_v( (x1+pert) , (x2) ) - fcn_v( x1 , (x2+pert) ) + fcn_v( x1 , x2 ) ) / ( pert*pert );
	h21 = h12;//( fcn_v( (x1+pert) , (x2+pert) ) - 2*fcn_v( (x2+pert) , x2 ) + fcn_v( x1 , x2 ) ) / ( pert*pert );
	return 0;
}