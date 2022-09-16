//
//  main.cpp
//  Memo_Ame_option
//
//  Created by Yu Gao on 12/4/17.
//  Copyright © 2017 Yu Gao. All rights reserved.
//
#define _CRT_SECURE_NO_WARNINGS
#include <iostream>
#include <iomanip>
#include <cmath>
#include <fstream>
#include <cstdlib>

using namespace std;

double up_factor, uptick_prob, risk_free_rate, strike_price,downtick_prob,notick_prob;
double initial_stock_price, expiration_time, volatility, R;
int no_of_divisions;

double max(double a, double b) {
    return (b < a ) ? a:b;
}
double american_call_option(int k, int i,double**ame_call_mem) {
    if (k == no_of_divisions){
        return max(0.0, (initial_stock_price*pow(up_factor,i) - strike_price));}
    if(ame_call_mem[k][no_of_divisions+i]==-1)
    {
        ame_call_mem[k][no_of_divisions+i] =
        max((initial_stock_price * pow(up_factor,i) - strike_price),
            (uptick_prob*american_call_option(k+1, i+1, ame_call_mem) +
             notick_prob*american_call_option(k+1, i, ame_call_mem) +
             downtick_prob*american_call_option(k+1, i-1, ame_call_mem))/R);
        return ame_call_mem[k][no_of_divisions+i];
    }
    else
    {
        return ame_call_mem[k][no_of_divisions+i];
    }
}

double american_put_option(int k, int i, double**ame_put_mem) {
    if (k == no_of_divisions)
        return max(0.0, (strike_price - initial_stock_price*pow(up_factor,i)));
    if(ame_put_mem[k][no_of_divisions+i]==-1){
        ame_put_mem[k][no_of_divisions+i] =
        max((strike_price - initial_stock_price*pow(up_factor,i)),
            (uptick_prob*american_put_option(k+1, i+1, ame_put_mem)+
             notick_prob*american_put_option(k+1, i, ame_put_mem)+
             downtick_prob*american_put_option(k+1,i-1, ame_put_mem))/R);
        return ame_put_mem[k][no_of_divisions+i];
    }
    else
    {
        return ame_put_mem[k][no_of_divisions+i];
    }
}

int main (int argc, char* argv[])
{
    sscanf (argv[1], "%lf", &expiration_time);
    sscanf (argv[2], "%d", &no_of_divisions);
    sscanf (argv[3], "%lf", &risk_free_rate);
    sscanf (argv[4], "%lf", &volatility);
    sscanf (argv[5], "%lf", &initial_stock_price);
    sscanf (argv[6], "%lf", &strike_price);
    
    up_factor = exp(volatility*sqrt(2*expiration_time/((double) no_of_divisions)));
    R = exp(risk_free_rate*expiration_time/((double) no_of_divisions));
    uptick_prob = pow((sqrt(R)-1/sqrt(up_factor))/(sqrt(up_factor)-1/sqrt(up_factor)),2);
    downtick_prob = pow((sqrt(up_factor)-sqrt(R))/(sqrt(up_factor)-1/sqrt(up_factor)),2);
    notick_prob = 1 - downtick_prob - uptick_prob;
    cout << "Recursive Trinomial American Option Pricing" << endl;
    cout << "Expiration Time (Years) = " << expiration_time << endl;
    cout << "Number of Divisions = " << no_of_divisions << endl;
    cout << "Risk Free Interest Rate = " << risk_free_rate << endl;
    cout << "Volatility (%age of stock value) = " << volatility*100 << endl;
    cout << "Initial Stock Price = " << initial_stock_price << endl;
    cout << "Strike Price = " << strike_price << endl;
    cout << "--------------------------------------" << endl;
    cout << "Up Factor = " << up_factor << endl;
    cout << "Uptick Probability = " << uptick_prob << endl;
    cout << "Notick Probability = " << notick_prob << endl;
    cout << "Downtick Probability = " << downtick_prob << endl;
    cout << "--------------------------------------" << endl;
    
    double** ame_call_mem;
    ame_call_mem = new double * [no_of_divisions];
    for (int i = 0; i < no_of_divisions; i++)
        ame_call_mem[i] = new double[2 * no_of_divisions + 1];
    for (int i = 0; i < no_of_divisions; i++)
        for (int j = 0; j < 2 * no_of_divisions+1; j++)
            ame_call_mem[i][j] = -1;
    
    double** ame_put_mem;
    ame_put_mem = new double * [no_of_divisions];
    for (int i = 0; i < no_of_divisions; i++)
        ame_put_mem[i] = new double[2 * no_of_divisions + 1];
    for (int i = 0; i < no_of_divisions; i++)
        for (int j = 0; j < 2 * no_of_divisions + 1; j++)
            ame_put_mem[i][j] = -1;
    
    double call_price = american_call_option(0, 0,ame_call_mem);
    cout << "Trinomial Price of an American Call Option = " << call_price << endl;
    double put_price = american_put_option(0, 0, ame_put_mem);
    cout << "Trinomial Price of an American Put Option = " << put_price << endl;
    cout << "--------------------------------------" << endl;
    
    
    delete[] ame_put_mem; delete[] ame_call_mem;
}