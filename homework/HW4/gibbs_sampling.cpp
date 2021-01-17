//  main.cpp
//  gibbs_sampling
//  Created by Tianqi Wu on 4/12/18.
//  Copyright © 2018 Tianqi Wu. All rights reserved.

#include <iostream>
#include <iomanip>
#include <cmath>
#include <fstream>
#include <cstdlib>
#include <random>
#include <chrono>
#include "include.h"
#include "newmat.h"
#include "newmatio.h"

#define PI 3.141592654

using namespace std;

// cf http://www.cplusplus.com/reference/random/uniform_real_distribution/operator()/
// If you want to set a seed -- do it only after debug phase is completed
// otherwise errors will not be repeatable.
unsigned seed = (unsigned) std::chrono::system_clock::now().time_since_epoch().count();
default_random_engine generator (seed);

double get_gaussian(double mean, double standard_deviation)
{
    std::normal_distribution<double> distribution(mean, standard_deviation);
    double number = distribution(generator);
    return (number);
}

// this is just a hint -- if you are confused with what I suggest, you are welcome to take a
// different approach... no compulsion to use this thought-process, per se!
ColumnVector Gibbs_sampler_for_Multivariate_Gaussian(int index, ColumnVector Previous_value, SymmetricMatrix C, ColumnVector mean)
{
    // we want the conditional distribution of the "index"-th rv, conditioned on the
    // values of the other indices (i.e. "index"-th rv is random; all others are fixed).
    // See desription and linked YouTube Video from Prof. Niemi.
    
    // Took the formulae from http://fourier.eng.hmc.edu/e161/lectures/gaussianprocess/node7.html
    // Keep in mind that the expression for the correlation-matrix for the conditional probability
    // has a minor-typo at this reference (he has transposed "i" and "j"; you could have guessed it
    // by checking for dimensional consistency of the resulting correlation-matrix.
    
    // In the following -- I am assuming the original correlation matrix (Sigma) is partioned as
    // [Sigma11 Sigma12; Sigma21 Sigma22] (in MATLAB notation).  Sigma11 is an (n-1) x (n-1) matrix
    // Sigma12 is an (n-1)-Columnvector; Sigma21 is an (n-1)-Rowvector; Sigma22 is a scalar.
    // The import being -- we are keeping all the "top" (n-1)-many variables constant, and picking
    // a value for the bottom variable (i.e conditioning the bottom-variable on the previous (n-1)
    // variables.
    
    // Since the Gibbs-Sampler is not always going to sample the bottom-value, you need to go through
    // some clerical-work to make things work.

    int n;
    n = mean.nrows();
    Matrix S_11(n-1,n-1), S_12(n-1,1), S_21(1,n-1);
    double S_22;
    ColumnVector mean_1(n-1), x_j(n-1);
    
    // construct conditional x_j, mean_1 here
    for (int i = 1; i < n; i++){
        if(i<index){
            x_j(i) = Previous_value(i);
            mean_1(i) = mean(i);
        }
        else{
            x_j(i) = Previous_value(i+1);
            mean_1(i) = mean(i+1);
        }
    }
    
    // Write code to construct Sigma11 here
    for (int i = 1; i < n; i++){
        for (int j = 1; j < n; j++){
            if(i<index) S_11(i,j) = C(i,j);
            else{
                if (j>=index) S_11(i,j) = C(i+1,j+1);
            }
        }
    }
    
    // Write code to construct Sigma12, Sigma21 here
    for (int i = 1; i < n; i++){
        if(i<index) {
            S_12(i,1) = C(i,index);     //Sigma12
            S_21(1,i) = C(index,i);     //Sigma21
        }
        else{
            S_12(i,1) = C(i+1,index);
            S_21(1,i) = C(index,i+1);
            i +=1;
        }
    }

    // Write code to construct Sigma22 here
    S_22 = C(index,index);
    
    // Write some account-keeping code here... in my "teliology" https://www.merriam-webster.com/dictionary/teleology
    // xx is a ColumnVector, which starts off with xx = Previous_value; and then the "index"-th variable is an
    // appropriately generated Univariate-Normal RV (called "x")
    double cond_mean, cond_sigma;
    ColumnVector xx(n);
    xx = Previous_value;
    cond_mean = (mean(index) + S_21*(S_11.i())*(x_j-mean_1)).as_scalar();
    cond_sigma = (S_22- S_12.t()*(S_11.i())*S_12).as_scalar();
    xx(index) = get_gaussian(cond_mean, cond_sigma);
    return xx;
}

double Theoretical_PDF(ColumnVector x, SymmetricMatrix C, ColumnVector mean)
{
    // write C++ code that computes the expression of equation 1 of the assignment description
    double a,b;
    a = 1/sqrt((pow(2*PI,x.nrows()))*C.determinant());
    b = exp(-0.5*(((x-mean).t())*(C.i())*(x-mean)).as_scalar());
    return a*b;
}

int main (int argc, char* argv[])
{
    ColumnVector y_prev, y;
    Matrix count(100,100);
    int no_of_trials, dimension;
    double diff;
    clock_t time_before, time_after;
    
    // 2D case
    dimension = 2;
    
    sscanf (argv[1], "%d", &no_of_trials);
    ofstream pdf_data(argv[2]);
    ofstream pdf_theory(argv[3]);
    
    // The covariance matrix
    SymmetricMatrix C(2);
    C(1,1) = 0.75;
    C(1,2) = 0.25;
    C(2,1) = 0.25;
    C(2,2) = 0.5;
    
    // The mean vector
    ColumnVector mean(2);
    mean(1) = 1.0;
    mean(2) = 2.0;
    
    cout << "Multivariate Gaussian Generator using Gibbs Sampling" << endl;
    cout << "Dimension = " << mean.nrows() << endl;
    cout << endl << "Mean Vector = " << endl << mean;
    cout << endl << "Covariance Matrix = " << endl << C;
    
    
    for (int i = 1; i <= 100; i++)
        for (int j = 1; j <= 100; j++)
            count(i,j) = 0.0;
    
    time_before = clock();
    
    y_prev = mean;
    for (int i = 0; i < no_of_trials; i++)
    {
        y = Gibbs_sampler_for_Multivariate_Gaussian(i%(mean.nrows())+1, y_prev, C, mean);
        for (int j = 1; j <= 100; j++) {
            for (int k = 1; k <= 100; k++) {
                if ( (y(1) >= ((double) (j-52)/10)) && (y(1) < ((float) (j-51)/10)) &&
                    (y(2) >= ((double) (k-52)/10)) && (y(2) < ((float) (k-51)/10)) )
                    count(j,k)++;
            }
        }
        y_prev = y;
    }
    
    for (int j = 1; j <= 100; j++) {
        for (int k = 1; k <= 100; k++) {
            if (k < 100)
                pdf_data << count(j,k)/((double) no_of_trials) << ", ";
            if (k == 100)
                pdf_data << count(j,k)/((double) no_of_trials) << endl;
        }
    }
    
    time_after = clock();
    diff = ((double) time_after - (double) time_before);
    cout << "practice:\t" << (diff/60/CLOCKS_PER_SEC) <<"mins" << endl;
    
    double x1, x2;
    for (int j = 1; j <= 100; j++) {
        x1 = ((double) (j-51)/10);
        for (int k = 1; k <= 100; k++) {
            x2 = ((double) (k-51)/10);
            ColumnVector x(2);
            x(1) = x1;
            x(2) = x2;
            if (k < 100)
                pdf_theory << Theoretical_PDF(x, C, mean)*0.01 << ", ";
            if (k == 100)
                pdf_theory << Theoretical_PDF(x, C, mean)*0.01 << endl;
        }
    }
}
