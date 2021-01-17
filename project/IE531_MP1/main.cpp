//  main.cpp
//  Generating Exponential and Normal RVs from an unknown discrete distribution

#include <stdio.h>

#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <random>
#include <chrono>
#include <cmath>
#include <fstream>

#define HEADS 1
#define TAILS 0
#define EVENT1 1
#define EVENT2 2
#define EVENT3 3

using namespace std;

double prob1, prob2, prob3;

int headtail;

unsigned seed = (unsigned) std::chrono::system_clock::now().time_since_epoch().count();
default_random_engine generator (seed);

double get_uniform() {
    std::uniform_real_distribution <double> distribution(0.0,1.0);
    double number = distribution(generator);
    return (number);
}

void assign_random_probabilities_to_three_events() {
    
    // write code here
    
    double a = get_uniform();
    double b = get_uniform();
    
    if (a <= b) {
        prob1 = a;
        prob2 = b - a;
        prob3 = 1 - b;
    }
    else {
        prob1 = b;
        prob2 = a - b;
        prob3 = 1 - a;
    }
}

int three_event_generator(double prob1, double prob2, double prob3) {
    
    // write code here
    
    double num = get_uniform();
    if (num <= prob1) {
        return EVENT1;
    }
    else if ((num <= prob2 + prob1) && (num > prob1)){
        return EVENT2;
    }
    else {
        return EVENT3;
    }
}

int simulate_fair_coin_toss_from_three_event_generator(double prob1, double prob2, double prob3) {
    
    // write code here
    
    // Head: [123, 132, 213]  Tail: [231, 312, 321]
    
    int trial[3];
    
    trial[0] = 0;
    trial[1] = 0;
    trial[2] = 0;
    
    while ((trial[0] == trial[1]) | (trial[0] == trial[2]) | (trial[1] == trial[2])) {
        
        trial[0] = three_event_generator(prob1, prob2, prob3);
        trial[1] = three_event_generator(prob1, prob2, prob3);
        trial[2] = three_event_generator(prob1, prob2, prob3);
    }
    
    if ((trial[0] == 1) | ((trial[0] == 2) && (trial[1] == 1))) {
        return HEADS;
    }
    else {
        return TAILS;
    }
    
}

double get_uniform_from_three_event_generator(double prob1, double prob2) {
    
    // write code here
    
    int side, num_trials = 32;
    long dec = 0, base = 1;
    
    double prob3 = 1 - prob1 - prob2;
    
    for (int i = 0; i < num_trials; i++) {
        side = simulate_fair_coin_toss_from_three_event_generator(prob1, prob2, prob3);
        dec += base * side;
        base = base * 2;
    }
    
    double result = dec / (pow(2, 32) - 1);
    return result;
}

double get_exp_from_unfair_coin(double lambda, double prob1, double prob2) {
    
    return ((-1.0/lambda)*log(get_uniform_from_three_event_generator(prob1, prob2)));
}

double get_gaussian_from_unfair_coin(double prob1, double prob2) {
    return (sqrt(-2.0*log(get_uniform_from_three_event_generator(prob1, prob2)))*cos(2*3.141592654*get_uniform_from_three_event_generator(prob1, prob2)));
}



int main (int argc, char * argv[])
{
    
    double lambda;
    long no_of_trials;
    
    sscanf (argv[1], "%ld", &no_of_trials);
    sscanf (argv[2], "%lf", &lambda);
    ofstream outfile_1 (argv[3]);
    ofstream outfile_2 (argv[4]);
    
    assign_random_probabilities_to_three_events();
    
    
    cout << "Generating i.i.d. Unit-Normals and Exponentials from a 3-Event Generator" << endl;
    cout << "Probability of EVENT1                     = " << prob1 << endl;
    cout << "Probability of EVENT2                     = " << prob2 << endl;
    cout << "Probability of EVENT3                     = " << prob3 << endl;
    cout << "Number of trials                          = " << no_of_trials << endl;
    cout << "Output File Name for the Unit-Normal Data = " << argv[3] << endl;
    cout << "Output File Name for the Exponential Data = " << argv[4] << endl;
    cout << "Rate of Exponential RVs                   = " << lambda << endl;
    
    // i.i.d. exponential r.v. generation
    int count[100];
    for (int i = 0; i < 100; i++)
        count[i] = 0;
    
    for (int i = 0; i < no_of_trials; i++)
    {
        double y = get_exp_from_unfair_coin(lambda, prob1, prob2);
        for (int j = 0; j < 100; j++)
            if (y < ((double) j/10))
                count[j]++;
    }
    
    // Exponential CDF -- experimental vs. theory
    for (int j = 0; j < 100; j++)
        outfile_2 << ((double) j/10) << ", " << ((double) count[j])/((double) no_of_trials) << ", " << (1.0-exp(-1.0*lambda*((double) j/10))) << endl;
    
    // i.i.d. unit-normal generation
    double cdf_gaussian[100];
    for (int i = 0; i < 100; i++)
    {
        count[i] = 0;
        cdf_gaussian[i] = 0.0;
    }
    
    for (int i = 0; i < no_of_trials; i++)
    {
        double y = get_gaussian_from_unfair_coin(prob1, prob2);
        for (int j = 0; j < 100; j++)
            if (y < ((float) (j-50)/10))
                count[j]++;
    }
    // Unit-Normal CDF -- experimental vs. theory
    for (int j = 1; j < 100; j++)
        cdf_gaussian[j] = cdf_gaussian[j-1] +
        ((1/sqrt(2*3.141592654))*exp(-1.0*(((double) (j-50)/10)*
                                           ((double) (j-50)/10))/2)/10.0);
    
    for (int j = 0; j < 100; j++)
        outfile_1 << ((double) (j-50)/10) << ", " <<
        ((double) count[j])/((double) no_of_trials) << ", " << cdf_gaussian[j]
        << endl;
}
