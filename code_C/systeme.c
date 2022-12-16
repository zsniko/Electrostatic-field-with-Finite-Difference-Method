//
//  systeme.c
//  ChampElectrique
//
//  Created by Zhifan Song on 22/11/2021.
//

#include "systeme.h"

System creerSystem(int dimX,int dimY){
    System res = {dimX,dimY}; // Number of X-grids and number of Y-grids
    
    // v = zeros(Nx,Ny)
    res.tension = calloc(res.Nx,sizeof(double*));
    for(int i=0;i<res.Nx;i++){
        res.tension[i] = calloc(res.Ny,sizeof(double));
    }
    // Initialize edge potentials
    // Top-Wall potential and Bottom Wall potential = 0
    for(int i=0;i<res.Ny;i++){
        res.tension[0][i] = 0;
        res.tension[res.Nx-1][i] = 0;
    }
    // Left-Wall potential and Right Wall potential = 0
    for(int i=0;i<res.Nx;i++){
        res.tension[i][0] = 0;
        res.tension[res.Nx-1][0] = 0;
    }
    
    // Initializing corner potentials
    res.tension[0][0] = 0.5*(res.tension[0][1] + res.tension[1][0]);
    res.tension[res.Nx-1][0] = 0.5*(res.tension[res.Nx-2][0] + res.tension[res.Nx-1][1]);
    res.tension[0][res.Ny-1] = 0.5*(res.tension[0][res.Ny-2] + res.tension[1][res.Ny-1]);
    res.tension[res.Nx-1][res.Ny-1] = 0.5*(res.tension[res.Nx-1][res.Ny-2] + res.tension[res.Nx-2][res.Ny-1]);
    
    // matrice des charges
    res.charge = calloc(res.Nx,sizeof(double*));
    for(int i=0;i<res.Nx;i++){
        res.charge[i] = calloc(res.Ny,sizeof(double));
    }
        
    return res;
}

void freeSystem(System *s){
    
    for(int i=0;i<s->Nx;i++){
        free(s->tension[i]);
    }
    free(s->tension);
    
    for(int i=0;i<s->Nx;i++){
        free(s->charge[i]);
    }
    free(s->charge);
    
    s->Nx = 0;
    s->Ny = 0;
}

void writefileMatrix(char *fName,double** v,int dimX,int dimY){
    FILE* fichier = fopen(fName,"w");
    
    if(fichier != NULL){
        
        for(int i=0;i<dimX;i++){
            for(int j=0;j<dimY;j++){
                fprintf(fichier,"%lf ",v[i][j]);
            }
            fprintf(fichier,"\n");
        }
        fclose(fichier);
    }
    else printf("Erreur d'ouverture du fichier!\n");
    
}

double** creerMatrice(int dimX,int dimY){
    
    double** res = calloc(dimX,sizeof(double*));
    for(int i=0;i<dimX;i++){
        res[i] = calloc(dimY,sizeof(double));
    }
    return res;
    
}

void afficherMatrice(double** mat,int dimX,int dimY){
    for(int i=0;i<dimX;i++){
        for(int j=0;j<dimY;j++){
            printf("%lf ",mat[i][j]);
        }
        printf("\n");
    }
}

void freeMatrice(double ***mat,int dimX){
    for(int i=0;i<dimX;i++){
        free((*mat)[i]);
    }
    free(*mat);
    *mat = NULL;
}

double** substractMatrix(double** mat1, double**mat2,int dimX,int dimY){
    double **mat = creerMatrice(dimX, dimY);
    for(int i=0;i<dimX;i++){
        for(int j=0;j<dimY;j++){
            mat[i][j] = mat1[i][j] - mat2[i][j];
        }
    }
    return mat;
}

double norm(double** v,int dimX, int dimY){
    double max;
    double *norm = calloc(dimX,sizeof(double));
    for(int i=0;i<dimX;i++){
        for(int j=0;j<dimY;j++){
            norm[i] = norm[i] + v[i][j];
        }
    }
    for(int a=0;a<dimX;a++){
        norm[a] = fabs(norm[a]);
    }
    max = norm[0];
    for(int k=0;k<dimX;k++){
        if(max<norm[k]) max = norm[k];
    }
    return max; 
}

double distancePoint(int x1, int y1, int x2, int y2){
    double distance = 0;
    distance = sqrt(pow((x1-x2),2)+pow((y1-y2),2));
    return distance;
}
