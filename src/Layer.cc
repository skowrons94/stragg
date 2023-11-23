#include "include/Layer.h"

Layer::Layer(int LayerElement)
{
 Element = LayerElement;
 ThicknessStep = 10; 
}

Layer::~Layer( )
{
}

//Evaluate the k and x factors
double Layer::MakeBeta(double E)
{
 double mp = 938272; //proton's mass in keV
 return std::sqrt(2 * mp * E + E * E) / (mp + E);
}

double Layer::Xi(double E)
{
 //Xi = 2.5507e-7 * Z * M * n / (A * b^2) [keV]
 double XiFactor = 2.5507e-7;
 double TotalCharge = Element;
 double TotalMass = srim.getElementMass(Element);
 double MolarMass = srim.getElementMass(Element);
 double Beta = this->MakeBeta(E);
 double f = (XiFactor * TotalCharge * MolarMass * this->ThicknessStep)/(TotalMass * Beta * Beta);
 return f;
}

double Layer::Emax(double E)
{
 //E = (2 * m * c^2 * b^2 * g^2)/(1 + 2*g*m/M)
 double me = 511; //electron's mass in keV
 double mp = 938272; //proton's mass in keV
 double Beta = this->MakeBeta(E);
 double Gamma = 1/std::sqrt(1-Beta*Beta);
 double f = (2 * me * Beta * Beta * Gamma * Gamma) / (1 + (2*Gamma*me)/mp);
 return f;
}

//Get the Bohr variance of the current layer
double Layer::GetGVL(double E)
{
 double BohrFactor = (8 * this->GetXi(E)) / (3);
 double Beta = this->MakeBeta(E);
 double me = 511; //electron mass in keV
 double SumIonization = srim.getElementIonization(Element) * std::log(2 * me * Beta * Beta / srim.getElementIonization(Element));
 return std::sqrt(BohrFactor * SumIonization);
}

//Get the Vavilov variance of the current layer
double Layer::GetVVL(double E)
{
 double Xi = this->GetXi(E);
 double Beta = this->GetBeta(E);
 double K = this->GetK(E);
 return std::sqrt(Xi*Xi*(1-Beta*Beta/2)/K);
}

double Layer::GetDEML(double E)
{
 return this->EvaluateZiegler(E) * this->ThicknessStep * 0.001;
}

// Evaluate the Ziegler function
double Layer::EvaluateZiegler(double AtEnergy)
{
 if ( AtEnergy == 0)
 {
    return 0;
 }
 else if ( AtEnergy >= 0 && AtEnergy < 25)  // Simplification to get Ziegler(0)=0
 {
    return srim.getElementValue1( Element )*std::pow(AtEnergy,srim.getElementValue2( Element )) + srim.getElementValue3( Element )*std::pow(AtEnergy,srim.getElementValue4( Element ));
 }
 else
 {
    double StoppingLow, StoppingHigh, Stopping;
    StoppingLow = srim.getElementValue1( Element )*std::pow(AtEnergy,srim.getElementValue2( Element )) + srim.getElementValue3( Element )*std::pow(AtEnergy,srim.getElementValue4( Element ));
    StoppingHigh = (srim.getElementValue5( Element )/std::pow(AtEnergy,srim.getElementValue6( Element )))*std::log( (srim.getElementValue7( Element )/AtEnergy) + (srim.getElementValue8( Element )*AtEnergy));
    Stopping = (StoppingLow*StoppingHigh)/(StoppingHigh+StoppingLow);
    return Stopping;
 }
}