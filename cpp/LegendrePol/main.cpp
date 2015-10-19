#include <iostream>
#include <fstream>
#include <cmath>
using namespace std;
#include <iomanip>
using std::setw;
#include "include\armadillo"
using namespace arma;

int main()
{
    int n;
    cout << "Please enter value of n:\n>";
    cin >> n;

    double x[n] = {};
    double w[n]={};

    double pi = 3.14159265359;
    double root, dev_L, z;
    int m  = (n + 2)/2;

    double *x_temp1;
    double *w_temp1;
    double *x_temp2;
    double *w_temp2;

    x_temp1 = x;
    w_temp1 = w;
    x_temp2 = x + n - 1;
    w_temp2 = w + n - 1;

    for (int i = 1; i<m;i++){
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
        * w_temp1 = 2/((1-root*root)*(dev_L*dev_L));
        *(x_temp1++) = root;
        *(x_temp2--) = -root;
        *(w_temp2--) = *(w_temp1++);
    }
    cout << "roots:" << setw(20) << "weight:" << endl;
    for (int i=0; i<n; i++ )
    {
       cout << x[i] << setw(15) << w[i] << endl;
    }

    ofstream myfile ("LegendrePol_n_5.txt");
        if (myfile.is_open())
        {
            myfile << "Computed roots and their weights of the n'th Legendre polynomial:" <<endl;
            myfile << "n = " << n << endl;
            myfile << "roots:" << setw(20) << "weight:" << endl;
            for(int i = 0; i<n ; i++){
            myfile << x[i] << setw(15) << w[i] << endl;
            }

        }

}



