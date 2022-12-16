//
//  Champ.h
//  ChampElectrique
//
//  Created by Zhifan Song on 17/12/2021.
//

#ifndef Champ_h
#define Champ_h

#include "systeme.h"

typedef struct{
    int Nx,Ny;
    double **Ex;
    double **Ey;
}ChampE;

ChampE creerChamp( int dimX,int dimY );
void freeChamp( ChampE *e );
ChampE gradientV( int dimX, int dimY, System v3, double h );
double** amplitudeE( ChampE champ );


#endif /* Champ_h */
