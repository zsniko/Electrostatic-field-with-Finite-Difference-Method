//
//  main.c
//  ChampElectrique
//
//  Created by Zhifan Song on 22/11/2021.
//
//  Les variables sont donnees dans le main, on pourrait ecrire un nouveau programme qui demande a l'utilisateur d'entrer les dimensions 
//  ainsi que la longueur, la position de la plaque (boucle do while pour redemander si chiffres impossibles dans les dimensions choisies... )

#include "systeme.h"
#include "methodes.h"
#include "Champ.h"

#define Jacobi "/Users/zhifansong/Desktop/Informatique/Projet/ChampElectrique/ChampElectrique/resultats/Vjacobi.txt"
#define GaussSeidel "/Users/zhifansong/Desktop/Informatique/Projet/ChampElectrique/ChampElectrique/resultats/Vgauss.txt"
#define Relaxation "/Users/zhifansong/Desktop/Informatique/Projet/ChampElectrique/ChampElectrique/resultats/Vrelax.txt"
#define ChampEx "/Users/zhifansong/Desktop/Informatique/Projet/ChampElectrique/ChampElectrique/resultats/Ex.txt"
#define ChampEy "/Users/zhifansong/Desktop/Informatique/Projet/ChampElectrique/ChampElectrique/resultats/Ey.txt"
#define ampE "/Users/zhifansong/Desktop/Informatique/Projet/ChampElectrique/ChampElectrique/resultats/E.txt"
#define Charge_GS "/Users/zhifansong/Desktop/Informatique/Projet/ChampElectrique/ChampElectrique/resultats/chargeGS.txt"
#define Charge_R "/Users/zhifansong/Desktop/Informatique/Projet/ChampElectrique/ChampElectrique/resultats/chargeR.txt"


int main(int argc, const char * argv[]) {
    
    int dimX = 101, dimY = 101; // dimensions X Y, par exemple 101 donc 102x102 avec les bords et 100x100 a l'interieur sans bords
    double tension_plaque_up = 100; // tension de la plaque en haut
    double tension_plaque_down = -100; // tension de la plaque en bas
    int longueur_plaque = 51; // la longeur des 2 plaques, nombres de points
    int position_plaque = 15;// l'endroit ou on met le condensateur, par rapport au point milieu de l'axe X
    double h = 1; // Pas
    // reconditionner la dimension en fonction de h
    /*if (h < 1) {
    dimX = floor(dimX/h) - 1;
    dimY = floor(dimY/h) - 1;
    }*/
    
    int k = 0, kgs = 0, ksr = 0; // compteurs d'iterations des 3 methodes utilisees respectivement

    double eps0 = 8.85419e-12; //permittivite du vide, multiplie par e-12 pour la valeur reelle!
    double epsr = 3; // permittivite relative du dielectrique dans le condensateur
    double eps = eps0*epsr;
    
    int lp = floor(longueur_plaque/2);
    // mid points
    int mpx = ceil(dimX/2);
    int mpy = ceil(dimY/2);
    int pp1 = mpx + position_plaque;
    int pp2 = mpx - position_plaque;
    
    // PARTIE 1 - CHAMP ELECTRIQUE DANS LE VIDE
    System vv1 = creerSystem(dimX, dimY);
    int kvv1 = 0,kvv2 = 0;// compteurs iterations
    vv1.charge[10][10] = 2e-9; // mettre une charge, par exemple a [10,10]
    double realV = (1/(4*M_PI*eps0))*vv1.charge[10][10]*h*h/distancePoint(12,12,10,10);
    vv1 = MethodeGaussSeidel(vv1, dimX, dimY, 0, 0, 0, 0, h, &kvv1, 0);
    /*printf("La valeur reelle: %lf\n",realV);
    printf("La valeur MN: %lf\n\n",vv1.tension[12][12]);*/
    System vv2 = creerSystem(dimX, dimY);
    vv2.charge[10][10] = 2e-9; // mettre une charge, par exemple a [10,10]
    vv2 = MethodeRelaxation(vv2, dimX, dimY, 0, 0, 0, 0, h, &kvv2, 0);
    printf("Nombre d'iterations methode Gauss-Seidel: %d\n",kvv1);
    printf("Nombre d'iterations methode Sur-Relaxation Successive : %d\n\n",kvv2);
    writefileMatrix(Charge_GS, vv1.tension,vv1.Nx,vv1.Ny);
    writefileMatrix(Charge_R, vv2.tension,vv2.Nx,vv2.Ny);
    freeSystem(&vv1);
    freeSystem(&vv2);
    
    // PARTIE 2 - CHAMP ELECTRIQUE _ CONDENSATEUR PLAN
    // Le system de potentiels cree par la methode Jacobi
    System v1 = creerSystem(dimX, dimY);
    v1 = MethodeJacobi(v1, dimX, dimY, longueur_plaque, position_plaque, tension_plaque_up, tension_plaque_down,h, &k, 1);
    // Le system de potentiels cree par la methode Gauss-Seidel
    System v2 = creerSystem(dimX, dimY);
    v2 = MethodeGaussSeidel(v2, dimX, dimY, longueur_plaque, position_plaque,tension_plaque_up, tension_plaque_down, h, &kgs, 1);
    // Le systeme de potentiels cree par la methode Sur-Relaxation Sucessive
    System v3 = creerSystem(dimX, dimY);
    v3 = MethodeRelaxation(v3, dimX, dimY, longueur_plaque, position_plaque, tension_plaque_up, tension_plaque_down, h, &ksr, 1);
    
    // Affichage nombres d'iterations
    printf("Nombre d'iterations methode Jacobi: %d\n",k);
    printf("Nombre d'iterations methode Gauss-Seidel: %d\n",kgs);
    printf("Nombre d'iterations methode Sucessive over-relexation : %d\n\n",ksr);
    
    // Calcul champ electrique en utilisant systems de potentiels calcules par les methodes numeriques
    ChampE champ_electrique_jacobi = gradientV(dimX,dimY,v1,h);
    ChampE champ_electrique_GS = gradientV(dimX,dimY,v2,h);
    ChampE champ_electrique_relax = gradientV(dimX,dimY,v3,h);
    
    // Calcul amplitude champ electrique
    double** E1 = amplitudeE(champ_electrique_jacobi);
    double** E2 = amplitudeE(champ_electrique_GS);
    double** E3 = amplitudeE(champ_electrique_relax);
    
    // Calcul de l'energie electrostatique d'un condensateur
    // on utilise les resultats calcules par la methodde de sur-relaxation sucessive
    double Ee = 0;// l'energie : variable - somme
    for(int i=pp2-1;i<=pp1-1;i++){ // parcourir l'endroit ou on a mis le condensateur
        for(int j=mpy-lp-1;j<=mpy+lp-1;j++){
            Ee = Ee + 0.5*eps*(E3[i][j]*E3[i][j])*h*h;
        }
    }
    // la valeur de la capacite calculee par la methode numerique:
    double Cmn = 1e12* 2*Ee/((tension_plaque_up-tension_plaque_down)*(tension_plaque_up-tension_plaque_down));
    // la valeur reelle de la capacite:
    double C = 1e12*eps*longueur_plaque/(pp1-pp2);  // multiplie par 1e12 pour afficahge
    
    printf("La valeur de la capacite calculee par la methode numerique :\n");
    printf("Cmn = %lf pF\n",Cmn);
    printf("La valeur reelle de la capacite:\n");
    printf("C = %lf pF\n",C);

    // Ecrire les donnees dans les fichiers txt pour ensuite importer dans Matalab pour visualisation
    writefileMatrix(Jacobi, v1.tension,v1.Nx,v1.Ny);
    writefileMatrix(GaussSeidel, v2.tension,v2.Nx,v2.Ny);
    writefileMatrix(Relaxation, v3.tension,v3.Nx,v3.Ny);
    writefileMatrix(ChampEx, champ_electrique_relax.Ex,champ_electrique_relax.Nx,champ_electrique_relax.Ny);
    writefileMatrix(ChampEy, champ_electrique_relax.Ey,champ_electrique_relax.Nx,champ_electrique_relax.Ny);
    writefileMatrix(ampE, E3,champ_electrique_relax.Nx,champ_electrique_relax.Ny);
    
    // Vider les matrices creees
    freeSystem(&v1);
    freeSystem(&v2);
    freeSystem(&v3);
    freeChamp(&champ_electrique_jacobi);
    freeChamp(&champ_electrique_GS);
    freeChamp(&champ_electrique_relax);

    return 0;
}
