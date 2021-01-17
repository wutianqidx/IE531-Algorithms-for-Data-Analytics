// Illustration of concepts from Linear Algebra using the NEWMAT library
// Written by Prof. R.S. Sreenivas for IE523: Financial Computation
// We want to check if there is (or, if there is no) solution to Ax = y
// Notice, in my example A is singular, and therefore has no inverse.

#include <iostream>
#include <iomanip>
#include <cmath>
#include "/Users/alexyu/Desktop/Welcome_to_IE531__Algorithms_for_Data_Analytics/newmat11/newmatap.h"
#include "/Users/alexyu/Desktop/Welcome_to_IE531__Algorithms_for_Data_Analytics/newmat11/newmat.h"
#include "//Users/alexyu/Desktop/Welcome_to_IE531__Algorithms_for_Data_Analytics/newmat11/newmatio.h"


#define EPSILON 0.0000001

// This routine computes the rank of a matrix by counting the number of
// non-zero singular-values of the matrix.
// See http://en.wikipedia.org/wiki/Singular_value_decomposition for details
// regarding the singular-value-decomposition of a matrix

int Rank(Matrix A, int no_of_rows_of_A, int no_of_cols_of_A)
{
    // The SVD routine in NEWMAT assume #rows >= #cols in A
    // if #rows < #cols, then we compute the SVD of the transpose of A
    if (no_of_rows_of_A >= no_of_cols_of_A)
    {
        DiagonalMatrix D(no_of_cols_of_A);
        SVD(A,D);
        int rank = 0;
        for (int i = 1; i <= no_of_cols_of_A; i++)
            if (abs(D(i)) > EPSILON)
                rank++;
        return (rank);
    }
    else {
        DiagonalMatrix D(no_of_rows_of_A);
        SVD(A.t(),D);
        int rank = 0;
        for (int i = 1; i <= no_of_rows_of_A; i++)
            if (abs(D(i)) > EPSILON)
                rank++;
        return (rank);
    }
    
}

int main (int argc, char* argv[])
{
    Matrix A(6,5);
    
    A(1,1) = 2; A(1,2) =  0; A(1,3) =  0; A(1,4) =  6; A(1,5) = 2;
    A(2,1) = 1; A(2,2) = -1; A(2,3) = -1; A(2,4) =  4; A(2,5) = 0;
    A(3,1) = 2; A(3,2) = -2; A(3,3) =  2; A(3,4) =  4; A(3,5) = 0;
    A(4,1) = 2; A(4,2) = -2; A(4,3) = 0; A(4,4) = 6; A(4,5) = 0;
    A(5,1) = -4; A(5,2) =  4; A(5,3) = -8; A(5,4) = -4; A(5,5) = 0;
    A(6,1) = 4; A(6,2) =  0; A(6,3) = -2; A(6,4) = 14; A(1,5) = 4;
    
    Matrix y(6,1);
    y(1,1) = 2; y(2,1) = 7; y(3,1) = 18; y(4,1) = 16; y(5,1) = -40; y(6,1) = 2;
    
    cout << "Illustration of concepts from Linear Algebra & NEWMAT" << endl;
    cout << "A matrix:" << endl;
    cout << setw(9) << setprecision(3) << A;
    cout << "y vector:" << endl;
    cout << setw(9) << setprecision(3) << y << endl;
    
    cout << "The determinant of A = " << A.Determinant() << endl;
    if (abs(A.Determinant()) < EPSILON)
        cout << "Matrix A does not have an inverse" << endl;
    else {
        cout << "Matrix A has an inverse" << endl;
    }
    
    cout << "Rank of A = " << Rank(A, 6, 5) << endl;
    
    Matrix B(6,6);
    B = (A | y);
    
    cout << endl << "(A | y) Matrix:" << endl;
    cout << setw(9) << setprecision(3) << B << endl;
    cout << "Rank of (A | y) = " << Rank(B, 6, 6) << endl;
    
    if (Rank(A, 6, 5) == Rank(B, 6, 6))
        cout << "There is a solution to Ax = y" << endl;
    else
        cout << "There is no solution to Ax = y" << endl;
    
    cout << "---------------------------------" << endl;
    
}


