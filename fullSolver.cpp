#include <iostream>
#include <cstdlib>
#include <cmath>
#include <vector>

using namespace std;

int conditionMatrix( vector<vector<double>*>* matrix );
int printMatrix ( vector<vector<double>*>* matrix );
int identityMatrix( vector<vector<double>*>* m , int rank );
int zeroMatrix( vector<vector<double>*>* m , int rank );
int calcAtomicVector( vector<vector<double>*>* m , vector<vector<double>*>* matrix, int row );
int copyMatrix( vector<vector<double>*>* copy , vector<vector<double>*>* original );
int matrixProduct( vector<vector<double>*>* result , vector<vector<double>*>* matrix_1 , vector<vector<double>*>* matrix_2 );
int vectorProduct( vector<double>* y , vector<vector<double>*>* matrix , vector<double>* x );
// remove if not required
int invertAtomicMatrix( vector<vector<double>*>* inverse , vector<vector<double>*>* matrix , int col );
// remove if not required
int forwardSubstitution( vector<double>* y , vector<vector<double>*>* L , vector<double>* b );
int backwardSubstitution( vector<double>* x , vector<vector<double>*>* U , vector<double>* y );

int main(int argc, char const *argv[])
{
	// using vector of vectors in heap to define full matrix
	vector< vector<double>*>* matrix = new vector< vector<double>*>;
	vector<double>* row = new vector<double>;
	(*row) = { 4 , 4 , 2 };
	(*matrix).push_back(row);
	row = new vector<double>;
	(*row) = { 2 , 3 , 3 };
	(*matrix).push_back(row);
	row = new vector<double>;
	(*row) = { 4 , 5 , 3 };
	(*matrix).push_back(row);

	// get some information about the matrix
	int rank = (*matrix).size();
	
	// condition the matrix through partial row pivoting
	conditionMatrix( matrix );

	// create a dummy copy of matrix
	vector< vector<double>*>* dummy = new vector< vector<double>*>;
	copyMatrix( dummy , matrix );

	// number of M matrices required = rank - 1;
	vector<vector< vector<double>*>*> M_matrices;
	for ( int i = 0 ; i < (rank - 1) ; i++ ){
		
		// first create an identity matrix
		vector< vector<double>*>* m = new vector< vector<double>*>;
		identityMatrix( m , rank );

		// second figure out the off-diagonal elements 
		calcAtomicVector( m , dummy , i );
		M_matrices.push_back(m); 

		// update dummy matrix for next iteration
		matrixProduct( dummy , m , dummy );

	}

	/*printMatrix(M_matrices[0]);
	cout << endl;
	printMatrix(M_matrices[1]);
	cout << endl;*/

	// multiply all the m matrices in reverse order
	vector<vector<double>*>* M = new vector< vector<double>*>;
	identityMatrix( M , rank );
	
	for ( int i = M_matrices.size() - 1 ; i >= 0 ; i-- ){
		matrixProduct( M , M , M_matrices[i] );
	}
	vector<double> b = {2 , 3 , 5 };
	vector<double> y2;
	vectorProduct( &y2 , M , &b );
	
	// define upper triangular matrix
	vector<vector<double>*>* U = new vector< vector<double>*>;
	matrixProduct( U , M , matrix );
	cout << "upper triangular matrix:" << endl;
	printMatrix(U);
	cout << endl;
	// define lower triangular matrix
	/*vector<vector<double>*>* L = new vector< vector<double>*>;
	identityMatrix( L , rank );
	for ( int i = 0 ; i < M_matrices.size() ; i++ ){
		invertAtomicMatrix( M_matrices[i] , M_matrices[i] , i );
		matrixProduct( L , L , M_matrices[i]);
	}
	printMatrix(L);
	cout << endl;
	
	// check if L*U returns the original matrix
	vector<vector<double>*>* check = new vector< vector<double>*>;
	matrixProduct( check , L , U );
	printMatrix(check);
	cout << endl;
	// forward substitution -> find Ux = L\b
	
	forwardSubstitution( &y , L , &b );
	for ( int i = 0 ; i < rank ; i++ ){
		cout << y[i] << "   ";
	}
	cout << endl;*/
	// backward substitution -> x = U\y
	vector<double> x;
	backwardSubstitution( &x , U , &y2 );
	for ( int i = 0 ; i < rank ; i++ ){
		cout << x[i] << "   ";
	}

	cout << endl;
	return 0;
}

int conditionMatrix( vector<vector<double>*>* matrix ){
// function to condition matrix through partial row pivoting
	// assume square matrix
	int rank = (*matrix).size();

	// iterate over each column
	for ( int col = 0 ; col < rank ; col++ ){

		double max = 0.0;
		int max_row = col;

		// search for the max value along current column
		for ( int i = col ; i < rank ; i++ ){
			if ( abs((*(*matrix)[i])[col]) > max ) {
				max = (*(*matrix)[i])[col];
				max_row = i;
			}
		}

		// swap row with max_row
		vector<double> temp = (*(*matrix)[col]);
		(*(*matrix)[col]) = (*(*matrix)[max_row]);
		(*(*matrix)[max_row]) = temp;
	}
	return 0;
}

int printMatrix ( vector<vector<double>*>* matrix ){
	int rank = (*matrix).size();
	for ( int i = 0 ; i < rank ; i++ ){
		for ( int j = 0 ; j < rank ; j++ ){
			cout << (*(*matrix)[i])[j] << "   ";
		}
		cout << endl;
	}
	return 0;
}

int identityMatrix( vector<vector<double>*>* m , int rank ){
	vector<double>* row;
	for ( int c = 0 ; c < rank ; c++ ){
		row = new vector<double>;
		for ( int r = 0 ; r < rank ; r++ ){
			if ( r==c ) (*row).push_back(1.0);
			else (*row).push_back(0.0);
		}
		(*m).push_back(row);
	}
	return 0;
}

int zeroMatrix( vector<vector<double>*>* matrix , int rank ){
	vector<double>*row;
	for( int c = 0 ; c < rank ; c++ ){
		row = new vector<double>;
		for ( int r = 0 ; r < rank ; r++ ){
			(*row).push_back(0.0);
		}
		(*matrix).push_back(row);
	}
	return 0;
}

int calcAtomicVector( vector<vector<double>*>* m , vector<vector<double>*>* matrix, int row ){
	int rank = (*matrix).size();
	double pivot = (*(*matrix)[row])[row];
	for ( int i = (row + 1) ; i < rank ; i++ ){
		(*(*m)[i])[row] = -(*(*matrix)[i])[row]/pivot;
	}
	return 0;
}

int copyMatrix( vector<vector<double>*>* copy , vector<vector<double>*>* original ){
	int rank = (*original).size();
	zeroMatrix( copy , rank );
	for ( int i = 0 ; i < rank ; i++ ){
		for ( int j = 0 ; j < rank ; j++ ){
			(*(*copy)[i])[j] = (*(*original)[i])[j];
		}
	}
	return 0;
}

int matrixProduct( vector<vector<double>*>* result , vector<vector<double>*>* matrix_1 , vector<vector<double>*>* matrix_2 ){
	// calculates result = matrix_1*matrix*2
	int rank = (*matrix_1).size();
	
	vector<vector<double>*>* temp = new vector<vector<double>*>;
	
	zeroMatrix( temp , rank );
	for ( int i = 0 ; i < rank ; i++ ){
		for ( int j = 0 ; j < rank ; j++ ){
			for ( int k = 0 ; k < rank ; k++ ){
				(*(*temp)[i])[j] += ((*(*matrix_1)[i])[k])*((*(*matrix_2)[k])[j]);
			}
		}
	}
	(*result) = (*temp);
	return 0;
}

int vectorProduct( vector<double>* y , vector<vector<double>*>* matrix , vector<double>* x ){
	// calculates y = matrix*x
	int rank = (*x).size();
	for ( int i = 0 ; i < rank ; i++ ){
		(*y).push_back(0.0);
	}
	for ( int i = 0 ; i < rank ; i++ ){
		for ( int j = 0 ; j < rank ; j++ ){
			(*y)[i] += ((*(*matrix)[i])[j])*( (*x)[j] );
		}
	}
	return 0;
}

int invertAtomicMatrix( vector<vector<double>*>* inverse , vector<vector<double>*>* matrix , int col ){
	// only works for unit atomic matrices
	int rank = (*matrix).size();
	vector<vector<double>*>* temp = new vector<vector<double>*>;
	zeroMatrix( temp , rank );
	// first set the diagonals to 1.0
	for ( int i = 0 ; i < rank ; i++ ){
		(*(*temp)[i])[i] = 1.0;
	}
	// then change the size of the atomic vector
	for ( int i = (col + 1) ; i < rank ; i++ ){
		(*(*temp)[i])[col] = -(*(*matrix)[i])[col]; 
	}
	(*inverse) = (*temp);
	return 0;
}

int forwardSubstitution( vector<double>* y , vector<vector<double>*>* L , vector<double>* b ){
	int rank = (*L).size();
	double temp = (*b)[0] /(*(*L)[0])[0];
	(*y).push_back(temp);

	double sum;

	for ( int i = 1 ; i < rank ; i++ ){
		sum = 0.0;
		for ( int j = 0 ; j < i ; j++ ){
			sum += (*y)[j]*(*(*L)[i])[j];
		}
		temp = ( (*b)[i] - sum )/(*(*L)[i])[i];
		(*y).push_back(temp);
	}
	return 0;
}

int backwardSubstitution( vector<double>* x , vector<vector<double>*>* U , vector<double>* y ){
	int rank = (*y).size();
	for ( int i = 0 ; i < rank ; i++ ){
		(*x).push_back(0.0);
	}
	int end = rank - 1;
	double temp = (*y)[end] /(*(*U)[end])[end];
	(*x)[end] = temp;
	double sum;

		for ( int i = end-1 ; i >= 0 ; i-- ){
		sum = 0.0;
		for ( int j = end ; j > i ; j-- ){
			sum += (*x)[j]*(*(*U)[i])[j];
		}
		temp = ( (*y)[i] - sum )/(*(*U)[i])[i];
		(*x)[i] = temp;
	}

	return 0;
}