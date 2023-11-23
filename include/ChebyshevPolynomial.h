// Chebyshev polynomial
class ChebyshevPolynomial
{
private:
double T(unsigned int Grade, double Argument);
public:
ChebyshevPolynomial(){};
double GetValue(unsigned int Grade, double Argument){return this->T(Grade,Argument);};
};