//
//  main.cpp
//  Johnson_Lindenstrauss
//
//  Code frame created by Ramavarapu Sreenivas on 4/3/17.
//  Student work done by Tianqi Wu on 3/1/18
//  Copyright Â© 2017 Ramavarapu Sreenivas. All rights reserved.
//

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

typedef std::pair<double,double> Box_Muller_Pair;

using namespace std;

// cf http://www.cplusplus.com/reference/random/uniform_real_distribution/operator()/
// If you want to set a seed -- do it only after debug phase is completed
// otherwise errors will not be repeatable.
unsigned seed = (unsigned) std::chrono::system_clock::now().time_since_epoch().count();
default_random_engine generator (seed);

// U.I.I.D. RV generator
double get_uniform()
{
    std::uniform_real_distribution <double> distribution(0.0, 1.0);
    double number = distribution(generator);
    return (number);
}

// Using the Gaussian generator that is part of the C++ STL
double get_gaussian(double mean, double standard_deviation)
{
    std::normal_distribution<double> distribution(mean, standard_deviation);
    double number = distribution(generator);
    return (number);
}

// student code for generating random projection

Matrix Generate_Random_Projection(int k, int d)
{
    double scale = sqrt(d/k);
    Matrix A(k, d);
    
    IdentityMatrix I(d);
    ColumnVector mean(d); mean = 0.0;
    
    RowVector u1(d);
    for (int i = 1; i <= d; i++) u1(i) = get_gaussian(0.0, 1.0);
    double u1_norm = u1.norm_Frobenius();
    A.row(1) = u1 / u1_norm;
    
    for (int i = 2; i <= k; i++) {
        RowVector u(d);
        for (int j = 1; j <= d; j++) u(j) = get_gaussian(0.0, 1.0);
        u = scale * u / u1_norm;
        A.row(i) = u;
    };
    
    return A;
}

// student code ends

int main(int argc, const char * argv[])
{
    int original_dimension, reduced_dimension, no_of_cols, no_of_trials;
    double epsilon, delta, diff;
    clock_t time_before, time_after;
    
    sscanf (argv[1], "%d", &original_dimension);
    sscanf (argv[2], "%d", &reduced_dimension);
    sscanf (argv[3], "%d", &no_of_cols);
    sscanf (argv[4], "%lf", &epsilon);
    sscanf (argv[5], "%lf", &delta);
    ifstream input_file(argv[6]);
    sscanf (argv[7], "%d", &no_of_trials);
    
    cout << "Johnson-Lindenstrauss Lemma Demo" << endl;
    cout << "Reading a (" << original_dimension << " x " << no_of_cols << ") Matrix from file '";
    cout << argv[6] << "'" << endl;
    cout << "Reduced Dimension = " << reduced_dimension << endl;
    cout << "epsilon = " << epsilon << endl;
    cout << "delta = " << delta << endl;
    cout << "Reduced Dimension (i.e. " << reduced_dimension << ") should be >= " << ceil(1/(epsilon*epsilon)*log(((double) no_of_cols)/delta));
    cout << " for the bound to hold with probability " << 1 - delta << endl;
    
    time_before = clock(); // recording time before we started to read data
    
    Matrix A(original_dimension, no_of_cols);
    
    // student code for reading input csv file
    for (int j = 1; j <= original_dimension; j++) {
        string line;
        getline(input_file, line);
        
        stringstream row(line);
        string entry;
        for (int i = 1; i <= no_of_cols; i++) {
            getline(row, entry, ',');
            double x;
            stringstream value(entry);
            value >> x;
            A(j,i) = x;
        }
    }
    // student code ends
    
    time_after = clock(); // recording time after testing is complete
    diff = ((double) time_after - (double) time_before);
    cout << "It took " << (diff/CLOCKS_PER_SEC)/60.0 << " minutes to read data from file '" << argv[6] << "'" << endl;
    // testing Johnson-Lindenstrauss Lemma
    int no_of_hits = 0;
    cout << "#Trails for the testing-phase = " << no_of_trials << endl;
    
    Matrix R = Generate_Random_Projection(reduced_dimension, original_dimension);
    
    // this is the reduced-dimension representation of the x's (i.e the matrix of y's)
    Matrix C(reduced_dimension, original_dimension);
    C = R*A;
    
    time_before = clock(); // recording time before the testing starts
    
    // student code for verification of JL-Lemma
    for (int i = 1; i <= no_of_trials; i++){
        int m = (no_of_cols - 1) * get_uniform() + 1;
        int n = (no_of_cols - 1) * get_uniform() + 1;
        
        ColumnVector x1 = A.column(m);
        ColumnVector x2 = A.column(n);
        ColumnVector y1 = C.column(m);
        ColumnVector y2 = C.column(n);
        
        double lower = (1 - epsilon) * (x1 - x2).norm_Frobenius();
        double upper = (1 + epsilon) * (x1 - x2).norm_Frobenius();
        double val = (y1 - y2).norm_Frobenius();
        
        if ( (val >= lower) && (val <= upper) ) no_of_hits++;
    }
    // student code ends
    
    time_after = clock(); // recording time after testing is complete
    diff = ((double) time_after - (double) time_before);
    cout << "It took " << (diff/CLOCKS_PER_SEC)/60.0 << " minutes for testing to be completed" << endl;
    
    cout << "Johnson-Lindenstrauss Lemma is satisfied " << no_of_hits << "-many times over ";
    cout << no_of_trials << " attempts" << endl;
    cout << "Empirical Probability = " << ((double) no_of_hits/no_of_trials) <<
            " (Theory says it should be at least: " << 1-delta << ")" << endl;
    cout << " Is the JL-Lemma verified:" << ((double) no_of_hits/no_of_trials  >= 1 - delta ? "True" : "False" )<< endl;
    return 0;
    
}
