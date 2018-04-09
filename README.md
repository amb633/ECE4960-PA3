# ECE4960-PA3 - Programming Assignment 3 Parameter Extraction from Least Squares Fitting

- *output_log.txt* contains the log of all checks, tests and the convergence runs
- *optimal_parameter_summary.txt* contains the summary of the convergence runs
- *report* contains our graphs and analysis of the results

**Overview** 
This programming assignment focuses on constructing a modular program that is able to identify optimal parameters using least-squares fitting on a large dataset. Two different convergence methods (Quasi Newton and Secant) were tested on both unnormalized and normalized datasets.

***************************************************************************
## Part 1 - Full Solver
***************************************************************************
**Overview:**
Implementation of functions for full matrix solvers using the Gaussian Elimination method.

**Documentation**
- vectors and matrices have to be created using the C++ stdlib
- vectors are `vector<double>` and matrices are `vector<vector<double>>`

The following functions are defined in the full solver files.  All input arguments are typically pointers.

- `fullSolver( x , A , b )` : parent function to solve a set of simultaneous equations of the form Ax = b where A is a martrix of coefficients, b is a vector of constants and x is a vector of unknowns
- `conditionMatrix( A , b )` : function for partial row pivoting to ensure matrix is sufficiently conditioned
- `printMatrix( A )` : overloaded utility function to print matrix or vector to console
- `identityMatrix( mat , rank )` : creates a identity matrix of size rank and stores it in mat
- `zeroMatrix( mat , rank )` : creates a zero matrix of size rank and stores it in mat
- `calAtomicVector( m , mat , row )` : calculates the pivoting matrix (m) required for gaussian elimination of element specified in mat(row,row)
- `copyMatrix( copy , original )` : creates a deep copy of the matrix specified by original
- `matrixProduct( result , matrix_1 , matrix_2 )` : calculates result = matrix_1 * matrix_2
- `vectorProduct( result , matrix , vector )` : caluclates result = matrix * vector
- `backwardSubstitution( x , upperMatrix , y )` : finds the solution of Ux = L\b where U is the upper triangular matrix and L is the lower triangular matrix from LU decomposition 


***************************************************************************
## Part 2 - Validation of Parameter Extraction Program with Power Law
***************************************************************************

***************************************************************************
## Part 3 - Visualization of Data
***************************************************************************
Graphs were generated using MATLAB. Please refer to the *report* document for the graphs.

***************************************************************************
## Part 4A - Unnormalized Quasi Newton Convergence
***************************************************************************

***************************************************************************
## Part 4B - Unnormalized Secant Convergence
***************************************************************************

***************************************************************************
## Part 5 - Normalized Quasi Newton and Secant Convergence
***************************************************************************

***************************************************************************
## Part 6 - Various Initial Guess for all Four Convergence Schemes
***************************************************************************

***************************************************************************
## Part 7 - Vizualization of Converged Datasets
***************************************************************************
Graphs were generated using MATLAB. Please refer to the *report* document for the graphs and analysis