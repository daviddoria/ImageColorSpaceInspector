#ifndef Conversions_H
#define Conversions_H


void RGBtoCIELab(unsigned char rgb[3], float cieLab[3]);

void RGBToShapedHSV(float rgb[3], float hsv[3]);

#endif
