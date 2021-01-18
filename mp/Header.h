//
//  Header.h
//  IE531_MP1
//
//  Created by Tianqi Wu on 2/15/18.
//  Copyright © 2018 Tianqi Wu. All rights reserved.
//

#ifndef Header_h
#define Header_h


class Time
{
    private :
    
    public :
    double get_uniform();
    int three_event_generator(double prob1, double prob2, double prob3) ;
    int simulate_fair_coin_toss_from_three_event_generator(double prob1, double prob2, double prob3);
    double get_uniform_from_three_event_generator(double prob1, double prob2);
    double get_exp_from_unfair_coin(double lambda, double prob1, double prob2);
    
};

#endif
