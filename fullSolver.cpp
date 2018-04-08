#include "fullSolver.hpp"

void fullSolver( vector<double>* solution , vector< vector<double>*>* matrix , vector<double>* b ){
	
	// get some information about the matrix
	int rank = (*matrix).size();
	
	// condition the matrix through partial row pivoting
	conditionMatrix( matrix );
	//printMatrix( matrix );

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

	// multiply all the m matrices in reverse order
	vector<vector<double>*>* M = new vector< vector<double>*>;
	identityMatrix( M , rank );
	
	for ( int i = M_matrices.size() - 1 ; i >= 0 ; i-- ){
		matrixProduct( M , M , M_matrices[i] );
	}
	
	// no need to explicity calculate lower triangular matrix
		// Ux = (L^-1)*b where L = (M1^-1)(M2^-1)(M3^-1)...
		// so (L^-1) = ...(M3)(M2)(M1) -> in reverse order (we're already calculated this for the U matrix)
		// so Ux = y = M*b
	vector<double>* y = new vector<double>;
	vectorProduct( y , M , b );
	
	// define upper triangular matrix
	vector<vector<double>*>* U = new vector< vector<double>*>;
	matrixProduct( U , M , matrix );

	// backward substitution -> x = U\y
	backwardSubstitution( solution , U , y );

	return;
}

void conditionMatrix( vector<vector<double>*>* matrix ){
// function to condition matrix through partial row pivoting
	// assume square matrix
	int rank = (*matrix).size();

	// iterate over each column
	for ( int col = 0 ; col < rank ; col++ ){
		//printMatrix(matrix);

		double max = 0.0;
		int max_row = col;

		// search for the max value along current column
		for ( int i = col ; i < rank ; i++ ){
			double v = abs((*(*matrix)[i])[col]);
			if ( v > max ) {
				max = v;
				max_row = i;
			}
		}
		// swap row with max_row
		vector<double> temp = (*(*matrix)[col]);
		(*(*matrix)[col]) = (*(*matrix)[max_row]);
		(*(*matrix)[max_row]) = temp;
	}
//    return 0;
}

void printMatrix ( vector<vector<double>*>* matrix ){
	// overloaded function to print matrix
	int rank = (*matrix).size();
	for ( int i = 0 ; i < rank ; i++ ){
		for ( int j = 0 ; j < rank ; j++ ){
			cout << (*(*matrix)[i])[j] << "   ";
		}
		cout << endl;
	}
	cout << endl;
}

void printMatrix( vector<double>* vector ){
	// overloaded function to print vector
	int rank = (*vector).size();
	for ( int i = 0 ; i < rank ; i++ ){
		cout << (*vector)[i] << endl;
	}
	cout << endl;
//    return 0;
}

void identityMatrix( vector<vector<double>*>* m , int rank ){
	// function to create a identity matrix of size rank
	vector<double>* row;
	for ( int c = 0 ; c < rank ; c++ ){
		row = new vector<double>;
		for ( int r = 0 ; r < rank ; r++ ){
			if ( r==c ) (*row).push_back(1.0);
			else (*row).push_back(0.0);
		}
		(*m).push_back(row);
	}
//    return 0;
}

void zeroMatrix( vector<vector<double>*>* matrix , int rank ){
	// function to create a zero matrix of size rank
	vector<double>*row;
	for( int c = 0 ; c < rank ; c++ ){
		row = new vector<double>;
		for ( int r = 0 ; r < rank ; r++ ){
			(*row).push_back(0.0);
		}
		(*matrix).push_back(row);
	}
//    return 0;
}

void calcAtomicVector( vector<vector<double>*>* m , vector<vector<double>*>* matrix, int row ){
	// calculates the off diagonal elements for gaussian elimination where the current pivot is matrix[row][row]
	// the coefficients are stored in m
	int rank = (*matrix).size();
	double pivot = (*(*matrix)[row])[row];
	for ( int i = (row + 1) ; i < rank ; i++ ){
		(*(*m)[i])[row] = -(*(*matrix)[i])[row]/pivot;
	}
//    return 0;
}

void copyMatrix( vector<vector<double>*>* copy , vector<vector<double>*>* original ){
	// creates a deep copy of a matrix
	int rank = (*original).size();
	zeroMatrix( copy , rank );
	for ( int i = 0 ; i < rank ; i++ ){
		for ( int j = 0 ; j < rank ; j++ ){
			(*(*copy)[i])[j] = (*(*original)[i])[j];
		}
	}
//    return 0;
}

void matrixProduct( vector<vector<double>*>* result , vector<vector<double>*>* matrix_1 , vector<vector<double>*>* matrix_2 ){
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
//    return 0;
}

void vectorProduct( vector<double>* y , vector<vector<double>*>* matrix , vector<double>* x ){
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
//    return 0;
}

void backwardSubstitution( vector<double>* x , vector<vector<double>*>* U , vector<double>* y ){
	// backward sub to find x in Ux = y where U is upper triangular and y = (L^-1)b
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

//    return 0;
}
