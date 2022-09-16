//
//  card_game.h
//  cardgame
//
//  Created by Yu Gao on 9/17/17.
//  Copyright © 2017 Yu Gao. All rights reserved.
//

/*
 *  card_game.h
 *  Card Game
 *
 *  Created by Ramavarapu Sreenivas
 */

#ifndef	CARDGAME_H
#define CARDGAME_H
#include <algorithm>

using namespace std;

double value(int r, int b, double **mem)
{
    
    if (0 == r)
        return ((double) b);
    if (0 == b)
        return (0);
    
    if (mem[r][b] == -1)
        
    {
    double term1 = ((double) r/(r+b)) * value(r-1, b, mem);
    
    double term2 = ((double) b/(r+b)) * value(r, b-1, mem);
    
    mem[r][b] = max((term1 + term2), (double) (b - r));
    }
    
    //return std::max((term1 + term2), (double) (b - r));
    return mem[r][b];
}
#endif
