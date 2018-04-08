#include "utilityFunctions.hpp"

void readDataFile( string path , vector<double>* VGS , vector<double>* VDS , vector<double>* IDS ){
    // function to read and store file
    ifstream myfile( path );
    if (myfile.is_open()){
        myfile.ignore(2048 , '\n');
        while( !myfile.eof() ){
            double vgs , vds , ids;
            myfile >> vgs >> vds >> ids;
            (*VGS).push_back(vgs);
            (*VDS).push_back(vds);
            (*IDS).push_back(ids);
        }
        myfile.close();
        
        // note that the last entry is read twice when using eof
        // so use pop_back to get rid of extra entry
        (*VGS).pop_back();
        (*VDS).pop_back();
        (*IDS).pop_back();
    }
    else cout << "Unable to open file " << endl;
    return;
}

double sumSquares( vector<double>* s_model , vector<double>* s_measured , bool n ){
    // function to calculate the sum of squares
    // note: if you want the norm, you need to get the square root of sum
    if ( (*s_model).size() != (*s_measured).size() ){
        cout << "error is dimensions of s_vectors" << endl;
        return -1.0;
    }
    double sum = 0.0;
    if ( !n ){
        for ( int i = 0 ; i < (*s_model).size() ; i++ ){
            double temp = (*s_model)[i] - (*s_measured)[i];
            sum += temp*temp;
        }
    }
    else if ( n ){
        for ( int i = 0 ; i < (*s_model).size() ; i++ ){
            double temp = 1.0-((*s_model)[i])/((*s_measured)[i]);
            sum += temp*temp;
        }
    }
    return sum;
}

double calculateIds( double vgs , double vds , double kappa , double vth , double is ){
    // function to calculate the modelled ids values
    double vt = 26e-3;
    double term_1 = log(1.0 + exp( ( kappa*( vgs -  vth ) ) / ( 2*vt ) ));
    term_1 = term_1*term_1;
    double term_2 = log(1.0 + exp( ( kappa*(vgs-vth)-vds ) / (2*vt) ));
    term_2 = term_2*term_2;
    return is*(term_1 - term_2);
}

void modelIds( vector<double>* IDS_model , vector<double>* VGS , vector<double>* VDS , double kappa , double vth , double is ){
    for ( int i = 0 ; i < (*VGS).size() ; i++ ){
        double vgs = (*VGS)[i];
        double vds = (*VDS)[i];
        double value = calculateIds( vgs , vds , kappa , vth , is );
        (*IDS_model).push_back(value);
    }
    return;
}

void searchValues( vector<double>* s_measured , vector<double>* VGS , vector<double>* VDS , vector<double>* IDS , double vgs , double vds ){
    // ensure search values are in acceptable region
    if ( vgs < 0.5 || vgs > 5 ) {
        cout << "error in vgs search parameters" << endl;
        return;
    }
    if ( vds < 0 || vds > 5 ){
        cout << "error in vds search parameters" << endl;
        return;
    }

    // find which region of VGS
    int start = 0;
    while ( (*VGS)[start] < vgs ){
        start++;
    }
    int end = 0;
    while ( (*VGS)[end] <= vgs ){
        end++;
    }
    // then search smaller region for vds
    int entry;
    for ( entry = start ; entry < end ; entry ++ ){
        if ( (*VDS)[entry] == vds ) break;
    }
    // then retrieve the corresponding ids
    double ids = (*IDS)[entry];
    (*s_measured).push_back(vgs);
    (*s_measured).push_back(vds);
    (*s_measured).push_back(ids);
    return;
}

void add_vectors( vector<double>* x,  vector<double>* dx,  vector<double>* sum ){
    if((*x).size() != (*dx).size()){
        cout << "vectors have different dimensions... cannot add them togther! " << endl;
    }else{
        for( int i = 0; i < (*x).size(); i++){
            (*sum).push_back((*x)[i] + (*dx)[i]);
        }
    }
}

double delta_norm_rel( vector<double>* delta_a, vector<double>* a ){
    double norm = 0.0;
    for(int i = 0; i<(*a).size(); i++){
        norm = norm + (pow((*delta_a)[i], 2)/pow((*a)[i], 2));
    }
    return norm;
}

double delta_norm_abs( vector<double>* delta_a){
    double norm = 0.0;
    for(int i = 0; i<(*delta_a).size(); i++){
        norm = norm + (pow((*delta_a)[i], 2));
    }
    return norm;
}

double parameterSensitivity( vector<double>* parameters , double pertubation , int which ,
    vector<double>* VGS , vector<double>* VDS , vector<double>* IDS ){
    vector<double>* IDS_model = new vector<double>;
    
    double kappa = (*parameters)[0];
    double vth = (*parameters)[1];
    double is = (*parameters)[2];

    modelIds( IDS_model , VGS , VDS , kappa , vth , is );
    double v = sumSquares( IDS_model , IDS ); 

    double kappa_new = kappa;
    double vth_new = vth;
    double is_new = is;

    switch( which ){
        case 0: kappa_new *= (1.0+pertubation); break;
        case 1: vth_new *= (1.0+pertubation); break;
        case 2: is_new *= (1.0+pertubation); break;
        default: cout << "no valid parameter chosen for sensitivity analysis" << endl; return 0.0;
    }

    (*IDS_model).erase((*IDS_model).begin(), (*IDS_model).end() );
    modelIds( IDS_model , VGS , VDS , kappa_new , vth_new , is_new );
    double v_new = sumSquares( IDS_model , IDS ); 

    double num = v_new/v;
    double den = kappa_new/kappa;
    return num/den;
}

void scaleVector( double scalar, vector<double>* a,  vector<double>* result){
    for( int i = 0; i<(*a).size(); i++){
        (*result).push_back((*a)[i]*scalar);
    }
}