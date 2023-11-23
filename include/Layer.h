#ifndef LAYER_H
#define LAYER_H

#include <cmath>

#include "include/SRIM.h"

// The atomic class to handle the layers calculations
class Layer
{
private:

SRIM srim;

int Element;
double ThicknessStep;

double Xi(double E);
double Emax(double E);
double MakeBeta(double E);

public:

Layer(int LayerElement);
~Layer( );

void setThicknessStep(double Step){ this->ThicknessStep = Step; };

double GetK(double E){return this->Xi(E)/this->Emax(E);};
double GetXi(double E){return this->Xi(E);};
double GetBeta(double E){return this->MakeBeta(E);};

double EvaluateZiegler(double AtEnergy);

double GetGVL(double E);
double GetVVL(double E);
double GetDEML(double E);
};

#endif