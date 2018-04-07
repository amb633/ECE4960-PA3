#include <iostream>
#include <cstdlib>
#include <cmath>
#include <vector>

using namespace std;

void matrixTests(int argc, char const *argv[])
{
	vector< vector<double> > matrix;
	vector<double> row1 = {1 , 0 , 0};
	vector<double> row2 = {0 , 1 , 0};
	matrix.push_back(row1);
	matrix.push_back(row2);
	matrix.push_back({0,0,1});

	for ( int i = 0 ; i < 3 ; i++ ){
		for ( int j = 0 ; j < 3 ; j++ ){
			cout << matrix[i][j] << "   ";
		}
		cout << endl;
	}
	cout << endl;

	vector<vector<double>*> matrix2;
	vector<double>* r1 = new vector<double>;
	(*r1).push_back(1);
	(*r1).push_back(0);
	(*r1).push_back(0);
	matrix2.push_back(r1);

	for ( int i = 0 ; i < 3 ; i++ ){
		cout << (*matrix2[0])[i] << "   ";
	}
	cout << endl << endl;

	vector<vector<double>*>* matrix3 = new vector<vector<double>*>;
	(*matrix3).push_back(r1);
	r1 = new vector<double>;
	(*r1).push_back(0);
	(*r1).push_back(1);
	(*r1).push_back(0);
	(*matrix3).push_back(r1);
	r1 = new vector<double>;
	(*r1).push_back(0);
	(*r1).push_back(0);
	(*r1).push_back(1);
	(*matrix3).push_back(r1);

	for ( int i = 0 ; i < 3 ; i++ ){
		for ( int j = 0 ; j < 3 ; j++ ){
			cout << (*(*matrix3)[i])[j] << "   ";
		}
		cout << endl;
	}
	cout << endl;
//    return 0;
}
