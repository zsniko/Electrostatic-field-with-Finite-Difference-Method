//
//  methodes.c
//  ChampElectrique
//
//  Created by Zhifan Song on 16/12/2021.
//

#include "methodes.h"

System MethodeJacobi(System v, int dimX, int dimY, int length_plate, int position_plate, double tension_plaque_up, double tension_plaque_down, double h, int *k, int plaque){
    
    double **Vn = creerMatrice(dimX, dimY); // Matrice memorisation pour la methode Jacobi
    *k = 0; // initialisation compteur nombre d'iteration
    double eps0 = 8.85419e-12; //8.85419e-12 permittivite du vide
    double epsr = 1;
    double eps = eps0*epsr;
    
    int lp = floor(length_plate/2);
    // mid points
    int mpx = ceil(dimX/2);
    int mpy = ceil(dimY/2);
    int pp1 = mpx + position_plate;
    int pp2 = mpx - position_plate;
    double ep = 1; // critere d‘arret 0.0001
    
    // Methode de jacobi
    while(ep>0.0001){
        
        Vn = v.tension;
        v.tension = creerMatrice(dimX, dimY);
        (*k)++; // MAJ compteur nb iteration
        
        // Calcul (sans les 4 points sur les bords)
        for(int i=1;i<=v.Nx-2;i++){
            for(int j=1;j<=v.Ny-2;j++){
                
                // plaque: + - 100V imposees pour chaque iteration
                if (plaque == 1){
                  for(int a=mpy-lp-1;a<=mpy+lp-1;a++){
                      v.tension[pp1-1][a] = tension_plaque_up;
                      v.tension[pp2-1][a] = tension_plaque_down;
                  }
                }
                
                v.tension[i][j] = 0.25*(Vn[i+1][j]+Vn[i-1][j]+Vn[i][j+1]+Vn[i][j-1]+h*h*v.charge[i][j]/eps);
            }
        }
        
        ep = norm(substractMatrix(v.tension, Vn, dimX, dimY),dimX,dimY)/norm(v.tension,dimX,dimY);
        
    }
    
    return v;
    
}

System MethodeGaussSeidel(System v, int dimX, int dimY, int length_plate, int position_plate, double tension_plaque_up, double tension_plaque_down, double h, int *k, int plaque){
    
    
    double **VnGs = creerMatrice(dimX, dimY); // Matrice memorisation pour la methode GS
    *k = 0; // initialisation compteur nombre d'iteration
    double eps0 = 8.85419e-12; //8.85419e-12 permittivite du vide
    double epsr = 1;
    double eps = eps0*epsr;
    
    int lp = floor(length_plate/2);
    // mid points
    int mpx = ceil(dimX/2);
    int mpy = ceil(dimY/2);
    int pp1 = mpx + position_plate;
    int pp2 = mpx - position_plate;
    double ep = 1; // critere d‘arret 0.0001
    
    // Methode de Gauss-Seidel
    while(ep>0.0001){
        
        VnGs = v.tension;
        v.tension = creerMatrice(dimX, dimY); // ESSENTIEL
        (*k)++; // MAJ compteur nb iteration
        
        // Calcul (sans les 4 points sur les bords)
        for(int i=1;i<=v.Nx-2;i++){  // boucle pour x
            for(int j=1;j<=v.Ny-2;j++){ // boucle pour y
                
                // plaque: + - 100V imposees pour chaque iteration
                if (plaque == 1){
                  for(int a=mpy-lp-1;a<=mpy+lp-1;a++){
                      v.tension[pp1-1][a] = tension_plaque_up;
                      v.tension[pp2-1][a] = tension_plaque_down;
                  }
                }
                
                // methode Gs - mise a jour
                 v.tension[i][j] = 0.25*(VnGs[i+1][j]+v.tension[i-1][j]+VnGs[i][j+1]+v.tension[i][j-1]+h*h*v.charge[i][j]/eps);
                
            }
        }
        ep = norm(substractMatrix(v.tension, VnGs, dimX, dimY),dimX,dimY)/norm(v.tension,dimX,dimY);
    }
    return v;
}

System MethodeRelaxation(System v, int dimX, int dimY, int length_plate, int position_plate, double tension_plaque_up, double tension_plaque_down, double h, int *k, int plaque){
    
    double **VnR = creerMatrice(dimX, dimY); // Memorisation matrice pour la methode relaxation
    *k = 0; // initialisation compteur nombre d'iteration
    double eps0 = 8.85419e-12; //8.85419e-12 permittivite du vide
    double epsr = 1;
    double eps = eps0*epsr;
    
    int lp = floor(length_plate/2); // la position de depart par rapport a l'axe X
    // mid points
    int mpx = ceil(dimX/2);
    int mpy = ceil(dimY/2);
    int pp1 = mpx + position_plate;
    int pp2 = mpx - position_plate;
    double ep = 1; // critere d‘arret 0.0001
    
    double omega = (8-sqrt(64-16*(cos(M_PI/dimX)+cos(M_PI/dimY))*(cos(M_PI/dimX)+cos(M_PI/dimY))))/((cos(M_PI/dimX)+cos(M_PI/dimY))*(cos(M_PI/dimX)+cos(M_PI/dimY))); // 1.939676
    
    // Methode de Sur-relaxation sucessive
    while(ep>0.0001){
        
        VnR = v.tension;
        v.tension = creerMatrice(dimX, dimY); // ESSENTIEL
        (*k)++; // MAJ compteur nb iteration
        
        // Calcul (sans les 4 points sur les bords)
        for(int i=1;i<=v.Nx-2;i++){
            for(int j=1;j<=v.Ny-2;j++){
                
                // plaque: + - 100V imposees pour chaque iteration
                if (plaque == 1){
                  for(int a=mpy-lp-1;a<=mpy+lp-1;a++){
                      v.tension[pp1-1][a] = tension_plaque_up;
                      v.tension[pp2-1][a] = tension_plaque_down;
                  }
                }
                
                 v.tension[i][j] = VnR[i][j]+omega*(0.25*(VnR[i+1][j]+v.tension[i-1][j]+VnR[i][j+1]+v.tension[i][j-1]+h*h*v.charge[i][j]/eps)-VnR[i][j]);
                
            }
        }
        ep = norm(substractMatrix(v.tension, VnR, dimX, dimY),dimX,dimY)/norm(v.tension,dimX,dimY);
        //printf("%lf\n",ep);
    }
    return v;
}
