#include <iostream>
#include "armadillo"
using namespace std;
using namespace arma;


double legendre(int n, double x){

        // This is eq 5.11 from chapter 5.3.1 recursive relation
    double L_plus = 1 , L = 0, L_minus;


    for (int i = 0; i< n; i++)
    {
        L_minus = L;
        L = L_plus;
        L_plus = (2*i +1)*x*L - i*L_minus ;
        L_plus /= i+1;

    }


    return L_plus;
}
int main()
{

double limit = 5;
double alpha = 2;
double weight = 1;
double x = 1;

double L0 = legendre(2,2);
cout<<L0<<endl;




    return 0;
}

void GaussianLegendre(double upperlimit, double lowerlimit, double x, int n){

    double a = ((lowerlimit-upperlimit)/2) ;
    double b = (lowerlimit + upperlimit)/2;
    double t = a*x + b;



}
