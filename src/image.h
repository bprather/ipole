/*
 * image.h
 *
 *  Created on: Sep 9, 2019
 *      Author: bprather
 */

#ifndef IMAGE_H
#define IMAGE_H

#include "decs.h"

/* imaging */
void make_ppm(REAL p[], int nx, int ny, REAL freq, char filename[]);
void rainbow_palette(REAL data, REAL min, REAL max, int *pRed,
                     int *pGreen, int *pBlue);

#endif // IMAGE_H
