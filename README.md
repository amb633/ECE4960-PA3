# ECE4960-PA3 - Programming Assignment 3 Parameter Extraction from Least Squares Fitting

- *output_log.txt* contains the log of all checks, tests and the convergence runs
- *optimal_parameter_summary.txt* contains the summary of the convergence runs
- *report* contains our graphs and analysis of the results

**Summary**  

This programming assignment focuses on constructing a modular program that is able to identify optimal parameters using least-squares fitting on a large dataset. Two different convergence methods (Quasi Newton and Secant) were tested on both unnormalized and normalized datasets.

***************************************************************************
### Part 0 - Utility Functions
***************************************************************************
**Overview:** Common utility functions that will be accessed and used commonly throughout the assignment

** Documentation** 
- `readDataFile()` : reads and parses *ouputNMOS.txt* file. stores the data points in individual vectors for VGS , VDS and IDS respectively
- `sumSquares( model , measured , bool )` : calculates the sum of squared difference between two vectors of the same size. `false` is default for the final bool argument. Pass in a true boolean value to calculate the nomalized sum of squared differences (as opposed to the unnormalized raw data)
- `modelIds()` : calucaltes a vector of the modelled Id given specificied parameters
- `add_vectors( x , y , sum )` : calculates element-wise additon: sum = x + y
- `delta_norm_rel( delta , a )` : calculates the relative residual 
- `delta_norm_abs( delta )` : calculates the absolute residual
- `parameterSensitivity()` : calcualtes the sensitivity of the function output to specified parameter
- `scaleVector( scalar , a , result )`: scales `a` by `scalar` and stores in `result`

***************************************************************************
### Part 1 - Full Solver
***************************************************************************
**Overview:** Implementation of functions for full matrix solvers using the Gaussian Elimination method.

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
### Part 2 - Validation of Parameter Extraction Program with Power Law
***************************************************************************
**Overview:** Functions for validating parameter extraction for the Power Law using random samples with 10-20% noise

**Documentation**
- `randomSamples(x, y, noise, y_noisey)`: calculates the x and y sample pairs with y returning a value that has 10-20% noise randomly added to it ;
- `LinearLSF(x, y_noisey, H_matrix, RHS)` : calculates the Hessaian matrix and RHS vector for the Power Law problem that has been altered to be a linear problem by taking the log of it. The Hessaian matrix and RHS vector are then used in the fullSolver function to achieve the parameter extraction of  m and c0 for the S_measured data.
***************************************************************************
### Part 3 - Visualization of Data
***************************************************************************
Graphs were generated using MATLAB. Please refer to the * ReportAndAnalysis.pdf * for the graphs under section "Task 3".

***************************************************************************
### Part 4A - Unnormalized Quasi Newton Convergence
***************************************************************************
**Overview:** Functions for parameter fitting using the Quasi Newton Convergence Method

**Documentation**
- `quasiNetwon_itr()` : calculates the new parameters for each Quasi Newton iteration. The function take the current parameter values and adds a delta value to each parameter to update them and return the new parameter values.
- `quasiNetwon_dx()` : calculates the delta for each parameter by calculating the Hessian Matrix and the gradient of the least squares function.
- `linear_search()` : calculate the optimal "t" to scale the deltas of the parameters. By scaling the deltas by "t", then adding the scaled deltas to the parameters, it will achieve the optimal least sqaures value for that iteration.
- `t_adjusted_sum_sq()` : calculates the least sqaure sum for each "t" value being evaluated in the line search to decide is a "t" is better than the previous "t" value.
***************************************************************************
### Part 4B - Unnormalized Secant Convergence
***************************************************************************
**Overview:** Functions for parameter fitting using the Secant Convergence Method

**Documentation**
- `secantConvergence()` : parent function for calling the secant convergence iteration given sets of inital guesses for the parameters. pass in `false` (default) for calculating unnormalized convergence, and `true` for calculating normalized convergence
- `recurrenceRelation()` : calculates third data point given two initial guesses using the recurrence relation
- `secantGradient()` : calculates the partial-derivative gradients using the Secant method
- `secantHession()` : calculates the partial-derivatibe hessians using the Secant method

***************************************************************************
### Part 5 - Normalized Quasi Newton and Secant Convergence
***************************************************************************
**Overview:** All functions are identical to Parts 4A and 4B, with the exception that the functions will habe to be called with a boolean `true` argument to indicate that the normalized values have to be calculated. Ensures maximum code reuse. 

***************************************************************************
### Part 6 - Various Initial Guesses for all Four Convergence Schemes
***************************************************************************
**Overview:** All functions are identical to Parts 4A and 4B. Initial gueses are passed in through three nested loops to test every combination possible. The loops are repeated four times for each of the following: unnormalized Quasi-Newton, unnormalized Secant , normalized Quasi-Newton and normalized Secant.

***************************************************************************
### Part 7 - Vizualization of Converged Datasets
***************************************************************************
Graphs were generated using MATLAB. Please refer to the * ReportAndAnalysis.pdf * document for the graphs and analysis under "Task 7"
