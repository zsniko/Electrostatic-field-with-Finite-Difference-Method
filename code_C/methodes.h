//
//  methodes.h
//  ChampElectrique
//
//  Created by Zhifan Song on 16/12/2021.
//

#ifndef methodes_h
#define methodes_h

#include "systeme.h"

System MethodeJacobi(System v, int dimX, int dimY, int length_plate, int position_plate, double tension_plaque_up, double tension_plaque_down, double h, int *k, int plaque);
System MethodeGaussSeidel(System v, int dimX, int dimY, int length_plate, int position_plate, double tension_plaque_up, double tension_plaque_down, double h, int *k, int plaque);
System MethodeRelaxation(System v, int dimX, int dimY, int length_plate, int position_plate, double tension_plaque_up, double tension_plaque_down, double h, int *k, int plaque);

#endif /* methodes_h */
