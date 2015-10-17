#include <iostream>
#include "include\armadillo"
using namespace std;
using namespace arma;

/*
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
*/

double GaussianLegendre(double upperlimit, double lowerlimit, double x[], double w[], int n){
// using info from 5.3.5
    double a = ((lowerlimit-upperlimit)/2) ;
    double b = (lowerlimit + upperlimit)/2;
    //double t = a*x + b;
    double pi = 3.14159265359;
    int m  = (n + 2)/2;
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


    for (int i = 1; i<m;i++){
        root = cos(pi * (4*i-1) / (4*n + 2 ));


    do{
            // This is eq 5.11 from chapter 5.3.1 recursive relation
            double L_plus = 1.0 , L = 0.0, L_minus;
        for (int j = 0; j < n; j++)
        {
        L_minus = L;
        L = L_plus;
        L_plus = (2.0*j +1)*root*L - j*L_minus ;
        L_plus /= j+1;
        }

        dev_L = (-n*root*L_plus + n*L)/(1-root*root);

        z = root;
        root  = z - L_plus/dev_L;                   // Newton's method
       } while(fabs(root - z) > pow(10,-6));


        * w_temp1 = a*2/((1-root*root)*(dev_L*dev_L));
        *(x_temp1++) = a*root+b;
        *(x_temp2--) = -a*root+b;
        *(w_temp2--) = *(w_temp1++);

    }

}

double int_function(double x1, double x2, double x3, double y1, double y2, double y3)
{
  int alpha = 2.0;
  double denominator = sqrt(pow((x1-y1),2)+pow((x2-y2),2)+pow((x3-y3),2));
  double exponent = exp(-2.0*alpha*(sqrt(x1*x1+x2*x2+x3*x3)+sqrt(y1*y1+y2*y2+y3*y3)));
  if (denominator < pow(10,-6)){return 0;}
  else return exponent/denominator;
}

double arybla(double x[])
{
    double *x_2;
    x_2 = x;
    for (int i = 0; i<10; i++)
    {
        *(x_2++) = x[i]+3;
    }
}

int main()
{

//double limit = 5;
//double alpha = 2;
//   weights and function values for the use of the gauss-legendre
//   method
     int n = 30;
     double x[n] = {};
     double w[n]={};
     double upper_limit = 5;
     double lower_limit = -5;


//   set up the mesh points and weights
GaussianLegendre(lower_limit,upper_limit,x,w,n);
/*
cout << "x" << endl;
for (int j = 0; j<n; j++)
{
    cout << x[j] << endl;
}
cout << "w=" << endl;
for (int j = 0; j<n; j++)
{
    cout << w[j] << endl;
}
*/

double int_GaussianLegendre = 0;

for (int f = 0; f<n; f++){
for (int g = 0; g<n; g++){
for (int i = 0; i<n; i++){
for (int j = 0; j<n; j++){
for (int k = 0; k<n; k++){
for (int l = 0; l<n; l++){
    int_GaussianLegendre += w[f]*w[g]*w[i]*w[j]*w[k]*w[l]*int_function(x[f],x[g],x[i],x[j],x[k],x[l]);
}}}}}}

cout << "result=" << int_GaussianLegendre << endl;

/*
double int_gauss = 0.;
for ( int i = 0;  i < n; i++){
   int_gauss+=w[i]*int_function(x[i]);
}

    for (int j = 0; j<10; j++)
    {
        cout << x[j] << endl;
    }

    arybla(x);

    for (int j = 0; j<10; j++)
    {
        cout << x[j] << endl;
    }
*/
    return 0;
}
