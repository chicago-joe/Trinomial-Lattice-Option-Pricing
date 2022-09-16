// trinomial.cpp
//
// Calculating the price of an American Option using Recursive
// Programming where the Probability Matrix is derived from the
// Trinomial Lattice.  
// Written by Prof. Sreenivas
//
// Template modified for pricing American Options using Trinomial Methodology
// Modified by Joseph Loss
//

#define _CRT_SECURE_NO_WARNINGS

#include <iostream>
#include <iomanip>
#include <cmath>
#include <fstream>
#include <cstdlib>

using namespace std;

double up_factor, uptick_prob, risk_free_rate, strike_price;
double downtick_prob, notick_prob;
double initial_stock_price, expiration_time, volatility, R;
int no_of_divisions;

double max(double a, double b) {
	return (b < a) ? a : b;
}

double american_call_option(int k, int i, double current_stock_price)
{
	if (k == no_of_divisions)
		return max(0.0, (current_stock_price - strike_price));
	else
		return (max((current_stock_price - strike_price),
		(uptick_prob*american_call_option(k + 1, i + 1, current_stock_price*up_factor) +
			(notick_prob*american_call_option(k + 1, i, current_stock_price)) +
			(downtick_prob*american_call_option(k + 1, i - 1, current_stock_price / up_factor))) / R));
}

double american_put_option(int k, int i, double current_stock_price) {
	if (k == no_of_divisions)
		return max(0.0, (strike_price - current_stock_price));
	else
		return (max((strike_price - current_stock_price),
		(uptick_prob*american_put_option(k + 1, i + 1, current_stock_price*up_factor) +
			(notick_prob*american_put_option(k + 1, i, current_stock_price)) +
			(downtick_prob*american_put_option(k + 1, i - 1, current_stock_price / up_factor))) / R));
}

int main(int argc, char* argv[])
{
	sscanf(argv[1], "%lf", &expiration_time);
	sscanf(argv[2], "%d", &no_of_divisions);
	sscanf(argv[3], "%lf", &risk_free_rate);
	sscanf(argv[4], "%lf", &volatility);
	sscanf(argv[5], "%lf", &initial_stock_price);
	sscanf(argv[6], "%lf", &strike_price);

	up_factor = exp(volatility*sqrt(2.0*(expiration_time / ((double)no_of_divisions))));
	R = exp(risk_free_rate*expiration_time / ((double)no_of_divisions));
	uptick_prob = pow(((sqrt(R) - (1 / sqrt(up_factor))) / (sqrt(up_factor) - (1 / sqrt(up_factor)))), 2.0);
	downtick_prob = pow(((sqrt(up_factor) - sqrt(R)) / (sqrt(up_factor) - (1 / sqrt(up_factor)))), 2.0);
	notick_prob = 1 - uptick_prob - downtick_prob;

	cout << "Recursive Trinomial American Option Pricing" << endl;
	cout << "Expiration Time (Years) = " << expiration_time << endl;
	cout << "Number of Divisions = " << no_of_divisions << endl;
	cout << "Risk Free Interest Rate = " << risk_free_rate << endl;
	cout << "Volatility (%age of stock value) = " << volatility * 100 << endl;
	cout << "Initial Stock Price = " << initial_stock_price << endl;
	cout << "Strike Price = " << strike_price << endl;
	cout << "--------------------------------------" << endl;
	cout << "R = " << R << endl;
	cout << "Up Factor = " << up_factor << endl;
	cout << "Up-tick Probability = " << uptick_prob << endl;
	cout << "Down-tick Probability = " << downtick_prob << endl;
	cout << "No-tick Probability = " << notick_prob << endl;
	cout << "--------------------------------------" << endl;
	double call_price = american_call_option(0, 0, initial_stock_price);
	cout << "Trinomial Price of an American Call Option (recursive) = " << call_price << endl;
	double put_price = american_put_option(0, 0, initial_stock_price);
	cout << "Trinomial Price of an American Put Option (recursive) = " << put_price << endl;
	cout << "--------------------------------------" << endl << endl;
	
	system("pause");
	
	/* 
	cout << "--------------------------------------" << endl;
	cout << "Let us verify the Put-Call Parity: S+P-C = Kexp(-r*T) for American Options" << endl;
	cout <<  initial_stock_price << " + " << put_price << " - " << call_price;
	cout << " = " << strike_price << "exp(-" << risk_free_rate << " * " << expiration_time << ")" << endl;
	cout << initial_stock_price + put_price - call_price << " ?=? " << strike_price*exp(-risk_free_rate*expiration_time) << endl;
	if (abs(initial_stock_price + put_price - call_price - strike_price*exp(-risk_free_rate*expiration_time)) <= 1e-3)
		cout << "Looks like Put-Call Parity holds within three decimal places" << endl;
	else
		cout << "Looks like Put-Call Parity does NOT hold" << endl;
	cout << "--------------------------------------" << endl;
	*/
}