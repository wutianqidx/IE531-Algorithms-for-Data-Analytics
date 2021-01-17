//
//  main.cpp
//  Matrix Sketch
//
//  Created by Ramavarapu Sreenivas on 3/29/17.
//  Copyright Â© 2017 Ramavarapu Sreenivas. All rights reserved.
//
//  Edo Liberty's paper is about using a streaming algorithm for
//  Matrix Sketching.  But we are going to implement/test this in
//  in a non-streaming context.  "A" is a large (n x m) data-matrix
//  where m >>> n. We want to find a smaller (n x k) "B" matrix, where
//  k = 2/epsilon, such that the Frobenius-Norm of (A*A.t() - B*B.t())
//  is <= epsilon * the Frobenius-Norm of (A*A.t())

#include <iostream>
#include <iomanip>
#include <cmath>
#include <fstream>
#include <cstdlib>
#include <algorithm>
#include "include.h"
#include "newmat.h"
#include "newmatio.h"
#include "newmatap.h"

using namespace std;

double Frobenious_Norm(Matrix Data)
{
    // write this part yourself
    Matrix sq_Data = Data*Data.t();
    double F_norm = 0;
    for (int i = 1; i <= sq_Data.nrows(); i++){
        for (int j = 1; j <= sq_Data.ncols(); j++)
        {
            F_norm = F_norm + pow(sq_Data(i,j),2);
        }
    }
    return sqrt(F_norm);
}

Matrix Matrix_Sketch(Matrix Data, double epsilon)
{
    // Edo Liberty's Matrix Sketch will have the same number of rows
    // as the original data; the #cols is ceil(2.0/epsilon)
    // write this part of the code yourself
    int cols_s = ceil(2.0/epsilon);
    if (cols_s >= Data.nrows()) cols_s = Data.nrows();
    Matrix U,Result(Data.nrows(), cols_s);
    DiagonalMatrix sigma;
    int j = 1;
    for (int i = 1; i <= Data.ncols(); i++){
        SVD(Result, sigma, U);
        if (j<=cols_s) Result.column(j)=Data.column(i);
        else{
            for (int k=1; k<=cols_s;k++){
                double newsigma = sqrt(pow(sigma(k,k),2) - pow(sigma(cols_s,cols_s),2));
                if (newsigma<=0) sigma(k,k) = 0;
                else sigma(k,k) = newsigma;
            }
            Result = U*sigma;
            j = cols_s;
            Result.column(j)=Data.column(i);
        }
        j = j+1;
    }
    return Result;
}

int main (int argc, char* argv[])
{
    int dimension, no_of_data_points;
    double epsilon;
    
    sscanf (argv[1], "%d", &dimension);
    sscanf (argv[2], "%d", &no_of_data_points);
    sscanf (argv[3], "%lf", &epsilon);
    ifstream input_file(argv[4]);
    ofstream output_file(argv[5]);
    
    Matrix Data(dimension, no_of_data_points);
    
    cout << "Edo Liberty's Matrix Sketching Algorithm" << endl;
    cout << "----------------------------------------" << endl;
    cout << "Original Data-Matrix has " << dimension << "-rows & " << no_of_data_points << "-cols" << endl;
    cout << "Epsilon = " << epsilon << " (i.e. max. of " << 100*epsilon << "% reduction of  Frobenius-Norm of the Sketch Matrix)"<< endl;
    cout << "Input File = " << argv[4] << endl;
    
    // Read the Data
    for (int i = 1; i <= dimension; i++)
        for (int j = 1; j <= no_of_data_points; j++)
        {
            double x;
            input_file >> x;
            Data(i,j) = x;
        }
    
    // Compute the Frobenius-Norm of the original Data-Matrix
    double Data_Forbenius_Norm = Frobenious_Norm(Data);
    cout << "Frobenius Norm of the (" << Data.nrows() << " x " << Data.ncols() << ") Data Matrix = ";
    cout << Data_Forbenius_Norm << endl;
    
    Matrix Sketch(dimension, min(dimension, (int) ceil(2/epsilon)));
    Sketch = Matrix_Sketch(Data, epsilon);
    double Sketch_Forbenius_Norm = Frobenious_Norm(Sketch);
    cout << "Frobenius Norm of the (" << Sketch.nrows() << " x " << Sketch.ncols() << ") Sketch Matrix = ";
    cout << Sketch_Forbenius_Norm << endl;
    cout << "Change in Frobenius-Norm between Sketch & Original  = ";
    cout << setprecision(3) << 100*(Sketch_Forbenius_Norm - Data_Forbenius_Norm)/Data_Forbenius_Norm << "%" << endl;
    
    output_file << Sketch;
    cout << "File `" << argv[5] << "' contains a (" << Sketch.nrows() << " x " << Sketch.ncols();
    cout << ") Matrix-Sketch" << endl;
    
    
}
