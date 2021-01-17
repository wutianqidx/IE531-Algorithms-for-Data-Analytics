//  main.cpp
//  HW3
//  Created by Tianqi Wu on 3/29/18.
//  Copyright © 2018 Tianqi Wu. All rights reserved.

#include <iostream>
#include <fstream>
#include <sstream>
#include <time.h>
#include <random>
#include <chrono>
#include <string.h>
#include <math.h>
#include "newmat.h"
#include "newmatio.h"
#include "newmatap.h"

using namespace std;
unsigned seed = (unsigned) std::chrono::system_clock::now().time_since_epoch().count();
default_random_engine generator (seed);

// U.I.I.D. RV generator
double get_uniform()
{
    std::uniform_real_distribution <double> distribution(0.0, 1.0);
    double number = distribution(generator);
    return (number);
}

Matrix get_matrix(int no_rows)
{
    Matrix A(no_rows, no_rows);
    for (int i = 1; i <= no_rows; i++) {
        for (int j = 1; j <= no_rows; j++) {
            if (get_uniform() < 0.5) A(i,j) = 5*get_uniform();
            else A(i,j) = -5*get_uniform();
        };
    };
    return A;
}

Matrix straight_multiplication(Matrix B, int exponent, int no_rows)
{
    IdentityMatrix I(no_rows);
    Matrix C = I;
    for (int i = 1; i<= exponent; i++){
        C = B*C;
    }
    return C;
}

Matrix repeated_squaring(Matrix B, int exponent, int no_rows)
{
    IdentityMatrix I(no_rows);
    if (exponent == 0) return I;
    if (exponent%2 != 0) return (B * repeated_squaring(B*B, (exponent-1)/2, no_rows));
    else return repeated_squaring(B*B, exponent/2, no_rows);
}

int main(int argc, const char * argv[]){
    int no_rows, exponent;
    double diff;
    clock_t time_before, time_after;
    
    sscanf (argv[1], "%d", &exponent);
    sscanf (argv[2], "%d", &no_rows);
    cout << "The number of rows/columns in the square matrix is: " << argv[2] << endl;
    cout << "The exponent is: " << argv[1] << endl;
    
    //repeated_squaring
    Matrix B = get_matrix(no_rows);
    time_before = clock();
    Matrix D = repeated_squaring(B,exponent,no_rows);
    time_after = clock();
    diff = ((double) time_after - (double) time_before);
    cout << "Repeated Squaring Result: " << endl;
    cout << "It took " << (diff/CLOCKS_PER_SEC) << " seconds to complete " << endl;
    //straight_multiplication
    time_before = clock();
    Matrix C = straight_multiplication(B,exponent,no_rows);
    time_after = clock();
    diff = ((double) time_after - (double) time_before);
    cout << "Direct Multiplication Result: " << endl;
    cout << "It took " << (diff/CLOCKS_PER_SEC) << " seconds to complete " << endl;
}
