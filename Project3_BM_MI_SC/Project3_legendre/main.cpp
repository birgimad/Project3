#include <iostream>
#include "include\armadillo"
using namespace std;
using namespace arma;

//   The following function computes
//   the weights and function values for the use of the gauss-legendre method
double GaussianLegendre(double upperlimit, double lowerlimit, double x[], double w[], int n){
// using info from 5.3.5 to change limits to go from -1 to 1
    double a = ((lowerlimit-upperlimit)/2) ;
    double b = (lowerlimit + upperlimit)/2;
    //double t = a*x + b;
    double pi = 3.14159265359;
    int m  = (n+1)/2;
    double root;
    double dev_L;
    double z;

    double *x_temp1;
    double *w_temp1;
    double *x_temp2;
    double *w_temp2;

    x_temp1 = x;
    w_temp1 = w;
    x_temp2 = x + n - 1;
    w_temp2 = w + n - 1;

    for (int i = 1; i<=m;i++){
        root = cos(pi * (4*i-1) / (4*n + 2 ));  //approximation of the root of the n'th polynomial
    do{
        // This is eq 5.11 from chapter 5.3.1 recursive relation to compute the Legendre polynomials
        double L_plus = 1.0 , L = 0.0, L_minus;
        for (int j = 0; j < n; j++)
        {
        L_minus = L;
        L = L_plus;
        L_plus = (2.0*j +1)*root*L - j*L_minus ;
        L_plus /= j+1;
        }

        dev_L = (-n*root*L_plus + n*L)/(1-root*root);   //derivative of the Legendre polynomial (L_plus)

        z = root;
        root  = z - L_plus/dev_L;                   // Newton's method
       } while(fabs(root - z) > pow(10,-6));

        //compute values of x and w
        * w_temp1 = a*2/((1-root*root)*(dev_L*dev_L));
        *(x_temp1++) = a*root+b;
        *(x_temp2--) = -a*root+b;
        *(w_temp2--) = *(w_temp1++);
    }
}

double int_function(double x1, double x2, double x3, double y1, double y2, double y3)
{
  int alpha = 2.0;
  double denominator = pow((x1-y1),2)+pow((x2-y2),2)+pow((x3-y3),2);
  double exponent = exp(-2.0*alpha*(sqrt(x1*x1+x2*x2+x3*x3)+sqrt(y1*y1+y2*y2+y3*y3)));
  if (denominator < pow(10,-6)){return 0;}
  else return exponent/sqrt(denominator);
}

int main()
{
     int n;
     double upper_limit;
     double lower_limit;
     cout << "Please enter value of n:\n>";
     cin >> n;
     cout << "Please enter value of lower integral limit:\n>";
     cin >> lower_limit;
     cout << "Please enter value of upper integral limit:\n>";
     cin >> upper_limit;
     cout << "n = " << n << endl;
     cout << "lower limit = " << lower_limit << endl;
     cout << "upper limit = " << upper_limit << endl;

     double x[n] = {};
     double w[n]={};

//   Computing the mesh points and weights
GaussianLegendre(lower_limit,upper_limit,x,w,n);

//  Computing the result from the computed mesh points and weights
double int_GaussianLegendre = 0;

for (int f = 0; f<n; f++){
for (int g = 0; g<n; g++){
for (int i = 0; i<n; i++){
for (int j = 0; j<n; j++){
for (int k = 0; k<n; k++){
for (int l = 0; l<n; l++){
    int_GaussianLegendre += w[f]*w[g]*w[i]*w[j]*w[k]*w[l]*int_function(x[f],x[g],x[i],x[j],x[k],x[l]);
}}}}}}

cout << "Integral value = " << int_GaussianLegendre << endl;


    return 0;
}
