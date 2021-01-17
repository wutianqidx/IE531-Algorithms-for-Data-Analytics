#include "generator.h"

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

