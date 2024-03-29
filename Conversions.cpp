#include "Conversions.h"

#include <algorithm> // min,max
#include <cmath>
#include <iostream>

#include <vtkMath.h>

void RGBtoCIELab(unsigned char rgb[3], float cieLab[3])
{
  // The input RGB must be between [0,255]
  // The output should be between [0, 100]

  // Normalize RGB values.
  double r = static_cast<double>(rgb[0])/static_cast<double>(255);
  double g = static_cast<double>(rgb[1])/static_cast<double>(255);
  double b = static_cast<double>(rgb[2])/static_cast<double>(255);

  if(r > 0.04045)
  {
    r = std::pow(((r + 0.055) / 1.055), 2.4);
  }
  else
  {
    r = r / 12.92;
  }

  if(g > 0.04045)
  {
    g = std::pow(((g + 0.055) / 1.055), 2.4);
  }
  else
  {
    g = g / 12.92;
  }

  if(b > 0.04045)
  {
    b = std::pow(((b + 0.055) / 1.055), 2.4);
  }
  else
  {
    b = b / 12.92;
  }

  r = r * 100.0;
  g = g * 100.0;
  b = b * 100.0;

  double X = r * 0.4124 + g * 0.3576 + b * 0.1805;
  double Y = r * 0.2126 + g * 0.7152 + b * 0.0722;
  double Z = r * 0.0193 + g * 0.1192 + b * 0.9505;

  double X2, Y2, Z2;
  X = X / 95.047;
  Y = Y / 100.0;
  Z = Z / 108.883;

  if(X > 0.008856)
  {
    X2 = std::exp(std::log(X)/3.0);
  }
  else
  {
    X2 = (7.787 * X) + (16.0 / 116.0);
  }
  if(Y > 0.008856)
  {
    Y2 = exp(log(Y)/3.0);
  }
  else
  {
    Y2 = (7.787 * Y) + (16.0 / 116.0);
  }
  if(Z > 0.008856)
  {
    Z2 = exp(log(Z)/3.0);
  }
  else
  {
    Z2 = (7.787 * Z) + (16.0 / 116.0);
  }

  cieLab[0] = (116.0 * Y2) - 16.0; // L
  cieLab[1] = 500.0 * (X2 - Y2); // a
  cieLab[2] = 200.0 * (Y2 - Z2); // b
  
//   if(cieLab[0] < 0 || cieLab[0] > 100 || cieLab[1] < 0 || cieLab[1] > 100 || cieLab[2] < 0 || cieLab[2] > 100)
//     {
//     std::cout << "Something is wrong!" << std::endl;
//     std::cout << "RGB " << rgb[0] << " " << rgb[1] << " " << rgb[2] << std::endl;
//     std::cout << "LAB " << cieLab[0] << " " << cieLab[1] << " " << cieLab[2] << std::endl;
//     }
}

void RGBToShapedHSV(float rgb[3], float hsv[3])
{
  if(rgb[0] < 0 || rgb[0] > 1 || rgb[0] < 0 || rgb[0] > 1 || rgb[0] < 0 || rgb[0] > 1)
    {
    std::cerr << "Error: RGBToShapedHSV input must be [0,1]!" << std::endl;
    exit(-1);
    }
  float original_hsv[3];
  vtkMath::RGBToHSV(rgb, original_hsv);

  float h = original_hsv[0];
  float s = original_hsv[1];
  float v = original_hsv[2];

  float r = s; // Radius of cylinder
  float theta = h; // Angle

  float lengthScale = 2.0f;
  float z = v * lengthScale; // The spacing of the Z/V slices of the cylinder are way too far apart without this scaling
  float x = r*cos(theta * 2.0f * vtkMath::Pi());
  float y = r*sin(theta * 2.0f * vtkMath::Pi());
  hsv[0] = x;
  hsv[1] = y;
  hsv[2] = z;
}
