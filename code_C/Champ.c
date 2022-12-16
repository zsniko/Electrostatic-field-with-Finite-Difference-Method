//
//  Champ.c
//  ChampElectrique
//
//  Created by Zhifan Song on 17/12/2021.
//

#include "Champ.h"

ChampE creerChamp(int dimX,int dimY){
    ChampE res = {dimX,dimY};
    res.Ex = creerMatrice(dimX, dimY);
    res.Ey = creerMatrice(dimX, dimY);
    return res;
}

void freeChamp(ChampE *e){
    for(int i=0;i<e->Nx;i++){
        free(e->Ex[i]);
    }
    for(int i=0;i<e->Nx;i++){
        free(e->Ey[i]);
    }
    free(e->Ex);
    free(e->Ey);
    e->Nx = 0;
    e->Ny = 0;
}

ChampE gradientV(int dimX, int dimY, System v,double h){
    
    ChampE E = creerChamp(dimX, dimY);
    // Ey
    // calcul points interieurs
    for(int i=0;i<dimX;i++){
        for(int j=1;j<dimY-1;j++){
            E.Ey[i][j] = -(v.tension[i][j+1] - v.tension[i][j-1])/(2*h);
        }
    }
    // points bords gauche et droite
    for(int i=0;i<dimX;i++){
        E.Ey[i][0] = -(v.tension[i][1] - v.tension[i][0])/h;
        E.Ey[i][dimY-1] = -(v.tension[i][dimY-2] - v.tension[i][dimY-1])/h;
    }
    //Ex
    // calcul points interieurs
    for(int j=0;j<dimY;j++){
        for(int i=1;i<dimX-1;i++){
            E.Ex[i][j] = -(v.tension[i+1][j] - v.tension[i-1][j])/(2*h);
        }
    }
    // points bords gauche et droite
    for(int j=0;j<dimY;j++){
        E.Ex[0][j] = -(v.tension[1][j] - v.tension[0][j])/h;
        E.Ex[dimX-1][j] = -(v.tension[dimX-2][j] - v.tension[dimX-1][j])/h;
    }
    return E;
}

double** amplitudeE(ChampE champ){
    double** ampE = creerMatrice(champ.Nx, champ.Ny);
    for(int i=0;i<champ.Nx;i++){
        for(int j=0;j<champ.Ny;j++){
            ampE[i][j] = sqrt(champ.Ex[i][j]*champ.Ex[i][j] + champ.Ey[i][j]*champ.Ey[i][j]);
        }
    }
    
    return ampE;
}
