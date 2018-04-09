#include "secant.hpp"

void recurrenceRelation( vector<double>* guess_2 , vector<double>* guess_1 , vector<double>* guess_0 ,
   vector<double>* VGS , vector<double>* VDS , vector<double>* IDS , bool normalize ){
   
   vector<double> IDS_model;
   modelIds( &IDS_model , VGS , VDS , (*guess_1)[0] , (*guess_1)[1] , (*guess_1)[2]);
   double v_1 = sumSquares( &IDS_model , IDS , normalize );
   
   IDS_model.erase( IDS_model.begin(), IDS_model.end()); 
   modelIds( &IDS_model , VGS , VDS , (*guess_0)[0] , (*guess_0)[1] , (*guess_0)[2]);
   double v_0 = sumSquares( &IDS_model , IDS , normalize );
   
   double kappa_2 = (((*guess_0)[0])*v_1 - ((*guess_1)[0])*v_0) / (v_1 - v_0);
   double vth_2 = (((*guess_0)[1])*v_1 - ((*guess_1)[1])*v_0) / (v_1 - v_0);
   double is_2 = (((*guess_0)[2])*v_1 - ((*guess_1)[2])*v_0) / (v_1 - v_0);
   
   (*guess_2).push_back(kappa_2);
   (*guess_2).push_back(vth_2);
   (*guess_2).push_back(is_2);
   return;
}

void secantConvergence( int& iterations , vector<double>* parameter_solutions , 
   double& absolute_residual , double& relative_residual , double& least_squares ,
   vector<double>* guess_0 , vector<double>* guess_1 , vector<double>* guess_2 , 
   vector<double>* VGS , vector<double>* VDS , vector<double>* IDS ,
   bool normalize ){

   // initialize vectors required for storing data
   vector<double> kappa_history;
   vector<double> vth_history;
   vector<double> is_history;
   vector<double> v_history;
   vector<double> IDS_model;

   // initialize some parameters
   double kappa , vth , is , v;

   // parse first guess and calculate model
   kappa = (*guess_0)[0];
   vth = (*guess_0)[1];
   is = (*guess_0)[2];
   modelIds( &IDS_model , VGS , VDS , kappa , vth , is );
   v = sumSquares( &IDS_model , IDS , normalize );
   // store values from the first guess
   kappa_history.push_back(kappa);
   vth_history.push_back(vth);
   is_history.push_back(is);
   v_history.push_back(v);
   // clear memory for second iteration
   IDS_model.erase(IDS_model.begin(), IDS_model.end()); 

   // repeat for second guess
   kappa = (*guess_1)[0];
   vth = (*guess_1)[1];
   is = (*guess_1)[2];
   modelIds( &IDS_model , VGS , VDS , kappa , vth , is );
   v = sumSquares( &IDS_model , IDS , normalize );

   kappa_history.push_back(kappa);
   vth_history.push_back(vth);
   is_history.push_back(is);
   v_history.push_back(v);
   IDS_model.erase(IDS_model.begin(), IDS_model.end()); 

   // repeat for third guess
   kappa = (*guess_2)[0];
   vth = (*guess_2)[1];
   is = (*guess_2)[2];
   modelIds( &IDS_model , VGS , VDS , kappa , vth , is );
   v = sumSquares( &IDS_model , IDS , normalize );
   kappa_history.push_back(kappa);
   vth_history.push_back(vth);
   is_history.push_back(is);
   v_history.push_back(v);
   IDS_model.erase(IDS_model.begin(), IDS_model.end()); 

   // initialize things for the first update
   vector<double> g;
   vector<vector<double>> h;
   vector<double> d;
   // calculate the gradient (delV)
   secantGradient( &g , &kappa_history , &vth_history , &is_history , &v_history );
   // calculate the hessians (jacobian)
   secantHessian( &h , &kappa_history , &vth_history , &is_history , &v_history );
   // call the solver to find the delta
   fullSolver( &d , &h , &g );
   // update the values for next iterate
   kappa -= d[0];
   vth -= d[1];
   is -= d[2];

   // now iterate until convergence
   int counter = 0;
   relative_residual = 1.0;
   //for ( int i = 0 ; i < 100 ; i++ ){
   while ( relative_residual > 1e-9 ){
       // calculate the current v
       modelIds( &IDS_model , VGS , VDS , kappa , vth , is );
       v = sumSquares( &IDS_model , IDS , normalize );

       // erase the "oldest" position
       kappa_history.erase(kappa_history.begin());
       vth_history.erase(vth_history.begin());
       is_history.erase(is_history.begin());
       v_history.erase(v_history.begin());
       IDS_model.erase(IDS_model.begin(), IDS_model.end());

       // push in the current position
       kappa_history.push_back(kappa);
       vth_history.push_back(vth);
       is_history.push_back(is);
       v_history.push_back(v);

       // calculate the update 
       vector<double> gradients;
       vector<vector<double>> hessians;
       vector<double> delta;

       secantGradient( &gradients , &kappa_history , &vth_history , &is_history , &v_history );
       secantHessian( &hessians , &kappa_history , &vth_history , &is_history , &v_history ); 
       
       fullSolver( &delta , &hessians , &gradients );
         
       // covergence criteria taken to be the relative residual
       absolute_residual = (delta[0])*(delta[0]) + (delta[1])*(delta[1]) + (delta[2])*(delta[2]);

       relative_residual = ((delta[0])*(delta[0]))/(kappa*kappa) + ((delta[1])*(delta[1]))/(vth*vth) + ((delta[2])*(delta[2]))/(is*is);

       // update variables for next iteration 
       kappa -= delta[0];
       vth -= delta[1];
       is -= delta[2];

       //cout << counter << " : " << v << " : " << absolute_residual << endl;
       counter++;
       if ( counter > 15000 || vth > 15 ) break;
   }
   iterations = counter;
   modelIds( &IDS_model , VGS , VDS , kappa , vth , is );
   least_squares = sumSquares( &IDS_model , IDS , normalize );
   (*parameter_solutions).push_back(kappa);
   (*parameter_solutions).push_back(vth);
   (*parameter_solutions).push_back(is);
   return;

}

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

void secantHessian( vector<vector<double>>* hessians , vector<double>* kappa_history , vector<double>* vth_history , vector<double>* is_history , vector<double>* v_history ) {

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

   double num = (v_0 - 2.0*v_1 + v_2);

   double h_kappa = num/((kappa_0 - kappa_1)*(kappa_1 - kappa_2));
   double h_vth = num/((vth_0 - vth_1)*(vth_1 - vth_2));
   double h_is = num/((is_0 - is_1)*(is_1 - is_2));
   zeroMatrix( hessians , 3 );
   ((*hessians)[0])[0] = h_kappa;
   ((*hessians)[1])[1] = h_vth;
   ((*hessians)[2])[2] = h_is;

   return;
}

