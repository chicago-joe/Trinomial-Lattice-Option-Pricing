// 
//	memoized_options_AM.cpp
//	Trinomial American Option Pricing using Recursion & Memoization
//	
//	Modified from Prof. Sreenivas' work on "Dynamic Trinomial Option Pricing"
//	
//	Created by Joseph Loss on 11/22/2018.
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

double american_call_option(int k, int i, double **memoized_call)
{
	if (k == no_of_divisions)
	{
		return max(0.0, (initial_stock_price * pow(up_factor, i) - strike_price));
	}
	if (memoized_call[k][no_of_divisions + i] == -1)
	{
		memoized_call[k][no_of_divisions + i] = 
		max((initial_stock_price * pow(up_factor, i) - strike_price),
			(uptick_prob * american_call_option(k+1, i+1, memoized_call) + 
				notick_prob * american_call_option(k+1, i, memoized_call) + 
				downtick_prob * american_call_option(k+1, i-1, memoized_call))/R);
		return memoized_call[k][no_of_divisions + i];	
	}
	else
	{
		return memoized_call[k][no_of_divisions + i];
	}
}

double american_put_option(int k, int i, double **memoized_put)
{
	if (k == no_of_divisions)
		{
			return max(0.0, (strike_price - initial_stock_price * pow(up_factor, i)));
		}
	if (memoized_put[k][no_of_divisions + i] == -1)
	{
		memoized_put[k][no_of_divisions + i] = 
		max((strike_price - initial_stock_price * pow(up_factor, i)),
			(uptick_prob * american_put_option(k+1, i+1, memoized_put) + 
				notick_prob * american_put_option(k+1, i, memoized_put) + 
				downtick_prob * american_put_option(k+1, i-1, memoized_put))/R);
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

// AMERICAN CALL OPTION (memoized)
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

// AMERICAN PUT OPTION (memoized)
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

	cout << endl << "Recursive Trinomial American Option Pricing (with memoization) " << endl;
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
	
	double call_price = american_call_option(0, 0, memoized_call);
	cout << "Trinomial Price of an American Call Option (memoized) = " << call_price << endl;
	double put_price = american_put_option(0, 0, memoized_put);
	cout << "Trinomial Price of an American Put Option (memoized) = " << put_price << endl;
	cout << "---------------------------------------" << endl << endl;
	
	delete[] memoized_call;
	delete[] memoized_put;

	system("pause");
}	