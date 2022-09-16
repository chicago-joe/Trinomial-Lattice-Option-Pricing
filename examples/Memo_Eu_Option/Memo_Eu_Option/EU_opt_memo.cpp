//
//  main.cpp
//  Memo_Eu_Option
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
    return (b < a )? a:b;
}
double N(const double& z) {
    if (z > 6.0) { return 1.0; }; // this guards against overflow
    if (z < -6.0) { return 0.0; };
    double b1 = 0.31938153;
    double b2 = -0.356563782;
    double b3 = 1.781477937;
    double b4 = -1.821255978;
    double b5 = 1.330274429;
    double p = 0.2316419;
    double c2 = 0.3989423;
    double a=fabs(z);
    double t = 1.0/(1.0+a*p);
    double b = c2*exp((-z)*(z/2.0));
    double n = ((((b5*t+b4)*t+b3)*t+b2)*t+b1)*t;
    n = 1.0-b*n;
    if ( z < 0.0 ) n = 1.0 - n;
    return n;
};

double option_price_put_black_scholes(const double& S,      // spot price
                                      const double& K,      // Strike (exercise) price,
                                      const double& r,      // interest rate
                                      const double& sigma,  // volatility
                                      const double& time)
{
    double time_sqrt = sqrt(time);
    double d1 = (log(S/K)+ r * time)/(sigma * time_sqrt) + 0.5 * sigma * time_sqrt;
    double d2 = d1 - (sigma * time_sqrt);
    return K*exp(-r*time)*N(-d2) - S*N(-d1);
};

double option_price_call_black_scholes(const double& S,       // spot (underlying) price
                                       const double& K,       // strike (exercise) price,
                                       const double& r,       // interest rate
                                       const double& sigma,   // volatility
                                       const double& time)    // time to maturity
{
    double time_sqrt = sqrt(time);
    double d1 = (log(S/K)+ r * time)/(sigma * time_sqrt) + 0.5 * sigma*time_sqrt;
    double d2 = d1-(sigma*time_sqrt);
    return S*N(d1) - K*exp(-r*time)*N(d2);
};



double european_call_option(int k, int i,double**eu_call_mem) {
    if (k == no_of_divisions){
        return max(0.0, (initial_stock_price*pow(up_factor, i) - strike_price));}
    if (eu_call_mem[k][no_of_divisions+i]==-1){
        eu_call_mem[k][no_of_divisions+i]=
        (uptick_prob*european_call_option(k+1, i+1, eu_call_mem) +
         notick_prob*european_call_option(k+1, i, eu_call_mem) +
         downtick_prob*european_call_option(k+1, i-1, eu_call_mem))/R;
        return eu_call_mem[k][no_of_divisions+i];
    }
    else
    {
        return eu_call_mem[k][no_of_divisions + i];
    }
}

double european_put_option(int k, int i,double**eu_put_mem) {
    if (k == no_of_divisions){
        return max(0.0, strike_price - (initial_stock_price*pow(up_factor, i)));}
    if (eu_put_mem[k][no_of_divisions+i]==-1){
        eu_put_mem[k][no_of_divisions+i]=
        (uptick_prob*european_put_option(k+1, i+1, eu_put_mem) +
         notick_prob*european_put_option(k+1, i,eu_put_mem) +
         downtick_prob*european_put_option(k+1,i-1,eu_put_mem))/R;
        return eu_put_mem[k][no_of_divisions+i];
    }
    else
    {
        return eu_put_mem[k][no_of_divisions+i];
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
    downtick_prob=pow((sqrt(up_factor)-sqrt(R))/(sqrt(up_factor)-1/sqrt(up_factor)),2);
    notick_prob=1-downtick_prob-uptick_prob;
    cout << "Recursive Binomial European Option Pricing" << endl;
    cout << "Expiration Time (Years) = " << expiration_time << endl;
    cout << "Number of Divisions = " << no_of_divisions << endl;
    cout << "Risk Free Interest Rate = " << risk_free_rate << endl;
    cout << "Volatility (%age of stock value) = " << volatility*100 << endl;
    cout << "Initial Stock Price = " << initial_stock_price << endl;
    cout << "Strike Price = " << strike_price << endl;
    cout << "--------------------------------------" << endl;
    cout << "Up Factor = " << up_factor << endl;
    cout << "Uptick Probability = " << uptick_prob << endl;
    cout << "--------------------------------------" << endl;
    
    
    double** eu_call_mem;
    eu_call_mem = new double * [no_of_divisions];
    for (int i = 0; i < no_of_divisions; i++)
        eu_call_mem[i] = new double[2 * no_of_divisions + 1];
    for (int i = 0; i < no_of_divisions; i++)
        for (int j=0;j<2*no_of_divisions+1;j++)
            eu_call_mem[i][j] = -1;
    
    double** eu_put_mem;
    eu_put_mem = new double * [no_of_divisions];
    for (int i = 0; i < no_of_divisions; i++)
        eu_put_mem[i] = new double[2 * no_of_divisions + 1];
    for (int i = 0; i < no_of_divisions; i++)
        for (int j = 0; j < 2 * no_of_divisions+1; j++)
            eu_put_mem[i][j]=-1;
    
    double eucall,euput,call_bs,put_bs;
    eucall = european_call_option(0,0,eu_call_mem);
    euput = european_put_option(0,0,eu_put_mem);
    call_bs = option_price_call_black_scholes(initial_stock_price, strike_price, risk_free_rate,
                                              volatility, expiration_time);
    put_bs = option_price_put_black_scholes(initial_stock_price, strike_price, risk_free_rate,
                                            volatility, expiration_time);

    cout <<"Trinomial Price of an European Call Option = "<< eucall <<endl;
    cout << "Call Price according to Black-Scholes = " << call_bs << endl;
    cout << "--------------------------------------" << endl;
    cout <<"Trinomial Price of an European Put Option = "<< euput <<endl;
    cout << "Put Price according to Black-Scholes = " << put_bs << endl;
    cout << "--------------------------------------" << endl;
    cout << "--------------------------------------" << endl;
    cout << "Verifying Put-Call Parity: S+P-C = Kexp(-r*T)" << endl;
    
    float a=strike_price*exp(-risk_free_rate*expiration_time);
    float b=initial_stock_price+euput-eucall;
    if (a == b)
    {
        cout <<  initial_stock_price << " + " << euput << " - " << eucall;
        cout << " = " << strike_price << "exp(-" << risk_free_rate << " * " << expiration_time << ")" << endl;
        cout << initial_stock_price + euput - eucall << " = " << strike_price*exp(-risk_free_rate*expiration_time) << endl;
    }
    else
    {
        cout <<  initial_stock_price << " + " << euput << " - " << eucall;
        cout << " != " << strike_price << "exp(-" << risk_free_rate << " * " << expiration_time << ")" << endl;
        cout << initial_stock_price + euput - eucall << " != " << strike_price*exp(-risk_free_rate*expiration_time) << endl;
        cout<<"It seems put-call parity does not hold."<<endl;
    }

    delete[] eu_put_mem; delete[] eu_call_mem;
}