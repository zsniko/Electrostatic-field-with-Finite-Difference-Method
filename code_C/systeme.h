//
//  systeme.h
//  ChampElectrique
//
//  Created by Zhifan Song on 22/11/2021.
//

#ifndef systeme_h
#define systeme_h

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "listeSC.h"

typedef struct{
    
    int Nx,Ny;
    double** tension;
    double** charge;
    
}System;


System creerSystem(int dimX,int dimY);
void freeSystem(System *s);
void writefileMatrix(char *fName,double** v,int dimX,int dimY);
double** creerMatrice(int dimX,int dimY);
void freeMatrice(double ***mat,int dimX);
void afficherMatrice(double** mat,int dimX,int dimY);
double norm(double** v,int dimX, int dimY);
double** substractMatrix(double** mat1, double**mat2,int dimX, int dimY);
double distancePoint(int x1, int y1, int x2, int y2);

#endif /* systeme_h */
