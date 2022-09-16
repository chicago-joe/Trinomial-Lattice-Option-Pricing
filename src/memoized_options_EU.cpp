// 
//	memoized_options_EU.cpp
//	Trinomial European Option Pricing using Recursion & Memoization
//	
//	Modified from Prof. Sreenivas' work on "Dynamic Trinomial Option Pricing"
//	
//	Created by Joseph Loss on 11/23/2018.
//  Copyright © 2018 Joseph Loss.
// 

#define _CRT_SECURE_NO_WARNINGS
#include <iostream>
#include <iomanip>
#include <cmath>
#include <fstream>
#include <cstdlib>
using namespace std;

double up_factor, uptick_prob, downtick_prob, notick_prob;
double risk_free_rate, strike_price, initial_stock_price;
double expiration_time, volatility, R;
int no_of_divisions;

double max(double a, double b)
{
	return (b < a) ? a : b;
}

// Overflow Protection // 
double N(const double& z) 
{ 
	if (z > 6.0) { return 1.0; };  
	if (z < -6.0) { return 0.0; }; 

	double b1 = 0.31938153; 
	double b2 = -0.356563782; 
	double b3 = 1.781477937; 
	double b4 = -1.821255978; 
	double b5 = 1.330274429; 
	double p = 0.2316419; 
	double c2 = 0.3989423; 
	double a = fabs(z); 
	double t = 1.0/(1.0 + a*p); 
	double b = c2 * exp((-z) * (z/2.0)); 
	double n = ((((b5 * t+b4) * t+b3) * t+b2) * t+b1) *t; 

	n = 1.0 - b*n; 
	if (z < 0.0) n = 1.0 - n;

	return n; 
}

double bsm_call_price(double& S, double& K, double& r, double& sigma, double& T)
{
	double d1 = (log(S/K) + r*T)/(sigma*sqrt(T)) + 0.5*sigma*sqrt(T);
	double d2 = d1 - sigma*sqrt(T);

	return S*N(d1) - K*exp(-r * T) * N(d2);
};

double european_call_option(int k, int i, double **memoized_call)
{
	if (k == no_of_divisions)
	{
		return max(0.0, (initial_stock_price * pow(up_factor, i) - strike_price));
	}
	if (memoized_call[k][no_of_divisions + i] == -1)
	{
		memoized_call[k][no_of_divisions + i] = 
		(uptick_prob * european_call_option(k+1, i+1, memoized_call) + 
			notick_prob * european_call_option(k+1, i, memoized_call) + 
			downtick_prob * european_call_option(k+1, i-1, memoized_call))/R;
		return memoized_call[k][no_of_divisions + i];	
	}
	else
	{
		return memoized_call[k][no_of_divisions + i];
	}
}

double bsm_put_price(double& S, double& K, double& r, double& sigma, double& T)
{
	double d1 = (log(S/K) + r*T)/(sigma*sqrt(T)) + 0.5*sigma*sqrt(T);
	double d2 = d1 - sigma*sqrt(T);

	return K*exp(-r * T) * N(-d2) - S*N(-d1);
};

double european_put_option(int k, int i, double **memoized_put)
{
	if (k == no_of_divisions)
		{
			return max(0.0, strike_price - initial_stock_price * pow(up_factor, i));
		}
	if (memoized_put[k][no_of_divisions + i] == -1)
	{
		memoized_put[k][no_of_divisions + i] = 
		(uptick_prob * european_put_option(k+1, i+1, memoized_put) + 
			notick_prob * european_put_option(k+1, i, memoized_put) + 
			downtick_prob * european_put_option(k+1, i-1, memoized_put))/R;
		return memoized_put[k][no_of_divisions + i];	
	}
	else
	{
		return memoized_put[k][no_of_divisions + i];
	}	
}



int main(int argc, char* argv[])
{
	sscanf(argv[1], "%lf", &expiration_time);
	sscanf(argv[2], "%d", &no_of_divisions);
	sscanf(argv[3], "%lf", &risk_free_rate);
	sscanf(argv[4], "%lf", &volatility);
	sscanf(argv[5], "%lf", &initial_stock_price);
	sscanf(argv[6], "%lf", &strike_price);

	up_factor = exp(volatility * sqrt(2.0 * expiration_time / ((double) no_of_divisions)));
	R = exp(risk_free_rate * expiration_time / ((double) no_of_divisions));
	
	uptick_prob = pow(((sqrt(R) - 1/sqrt(up_factor)) / 
		(sqrt(up_factor) - (1/sqrt(up_factor)))), 2.0);
	downtick_prob = pow(((sqrt(up_factor) - sqrt(R)) / 
		(sqrt(up_factor) - (1/sqrt(up_factor)))), 2.0);
	notick_prob = 1 - downtick_prob - uptick_prob;

/*	MEMOIZATION:
*	use a double pointer to store the result,
*	then re-use the result as needed,
*	rather than repeating previous calculations again.
*/

// EUROPEAN CALL OPTION (memoized)
	double **memoized_call;
	memoized_call = new double *[no_of_divisions];
	for (int i = 0; i < no_of_divisions; i++)
	{
		memoized_call[i] = new double[2 * no_of_divisions + 1];
	}
	for (int i = 0; i < no_of_divisions; i++)
	{
		for (int j = 0; j < 2 * no_of_divisions + 1; j++)
		{
			memoized_call[i][j] = -1;
		}
	}

// EUROPEAN PUT OPTION (memoized)
	double **memoized_put;
	memoized_put = new double *[no_of_divisions];
	for (int i = 0; i < no_of_divisions; i++)
	{
		memoized_put[i] = new double[2 * no_of_divisions + 1];
	}
	for (int i = 0; i < no_of_divisions; i++)
	{
		for (int j = 0; j < (2*no_of_divisions+1); j++)
		{
			memoized_put[i][j] = -1;
		}
	}

	cout << endl << "Recursive Trinomial European Option Pricing (with memoization) " << endl;
	cout << "Expiration Time (Years) = " << expiration_time << endl;
	cout << "Number of Divisions = " << no_of_divisions << endl;
	cout << "Risk Free Interest Rate = " << risk_free_rate << endl;
	cout << "Volatility (%age of stock value) = " << volatility * 100 << endl;
	cout << "Initial Stock Price = " << initial_stock_price << endl;
	cout << "Strike Price = " << strike_price << endl;	
	cout << "---------------------------------------" << endl;
	//	cout << "R = " << R << endl;
	cout << "Up Factor = " << up_factor << endl;
	cout << "Uptick Probability = " << uptick_prob << endl;
	cout << "Downtick Probability = " << downtick_prob << endl;
	cout << "Notick Probability = " << notick_prob << endl;
	cout << "---------------------------------------" << endl;
	
	double call_price, put_price, blackscholes_call_prc, blackscholes_put_prc; 
	call_price = european_call_option(0, 0, memoized_call);
	put_price = european_put_option(0, 0, memoized_put);

	blackscholes_call_prc = bsm_call_price(initial_stock_price, strike_price, 
		risk_free_rate, volatility, expiration_time); 
	blackscholes_put_prc = bsm_put_price(initial_stock_price, strike_price, 
		risk_free_rate, volatility, expiration_time);

	cout << "Trinomial Price of an European Call Option (memoized) = " << call_price << endl;
	cout << "Call Price according to Black-Scholes = " << blackscholes_call_prc << endl;
	cout << "---------------------------------------" << endl;
	cout << "Trinomial Price of an European Put Option (memoized) = " << put_price << endl;
	cout << "Put Price according to Black-Scholes = " << blackscholes_put_prc << endl;
	cout << "---------------------------------------" << endl << endl;
	
	delete[] memoized_call;
	delete[] memoized_put;

	system("pause");
}	