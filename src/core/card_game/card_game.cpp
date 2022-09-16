#include <iostream>
#include <cstdio>
#include <ctime>
#include "card_game.h"

using namespace std;
	
int main (int argc, char * const argv[]) 
{
	clock_t start;
	double duration;

	int total_number_of_cards;
	double **table;
	cout << endl;
	
	start = clock();
//	cin >> total_number_of_cards;
	sscanf_s (argv[1], "%d", &total_number_of_cards);
	
	table = new double *[total_number_of_cards/2 + 1];
	for (int i = 0; i < total_number_of_cards/2 + 1; i++)
	{
		table[i] = new double[total_number_of_cards/2 + 1];
		for (int j = 0; j < total_number_of_cards/2 + 1; j++)
		{
			table[i][j] = -1;
		}
	}

	cout << "Total Number of Cards = " << total_number_of_cards << endl;
	cout << "Value of the game = " << value(total_number_of_cards/2, total_number_of_cards/2, table) << endl;

	duration = (clock() - start) / (double) CLOCKS_PER_SEC;
	cout << "Time: " << duration << endl << endl;
	system("Pause");
    return 0;
}