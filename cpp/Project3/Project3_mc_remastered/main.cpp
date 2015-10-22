#include <cmath>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
using namespace std;


double func(double x1,double y1,double z1,double x2,double y2, double z3);

double probvalu(double x1,double x2, double x3, double x4, double x5, double x6, double N);

int main()
    {
        int N ;

        //X = new double[n];


        cout << "Read in the number of Monte-Carlo samples: " << endl;
        cin >> N;
        double *x = new double [6];
        double Montec = 0;
        double Montesqr = 0;
        double fx;
        double length = 5;
        double jacobidet = pow(exp(1),6);
        double px;

        //double x[N]={};

        srand(time(NULL));

        for (int i=0;i<N;i++){
            for (int j=0;j<6;j++){
                x[j] = ((double) rand() / (RAND_MAX)); //random numbers generated in the interval(0,1)

                //cout<<x[j]<<endl;

            }

            fx = func(x[0],x[1],x[2],x[3],x[4],x[5]);
            px = probvalu(x[0],x[1],x[2],x[3],x[4],x[5], ((double) N));
            Montec += fx/px/((double) N);
            Montesqr += (fx/px)*(fx/px);

       }

       double MC = Montec;  // /((double) N )
       double Montesqr1= Montesqr/((double) N);
       double variance= Montesqr1-MC*MC;
       double volume = pow(length,6);
       cout << " Integral = " << jacobidet*MC << endl;
       cout << "sigma = " << volume*sqrt(variance/ ((double) N)) << endl;


      return 0;
    }
        double func(double r1, double r2, double theta1, double theta2, double phi1, double phi2)
        {
           double alpha = 2.;

           double drdr = r1*r1*r2*r2*(-sin(theta1))*(-sin(theta2));
           double denominator = r1*r1+r2*r2-2*r1*r2*(cos(theta1)*cos(theta2)+sin(theta1)*sin(theta2)*cos(phi1-phi2));
           double exponential = exp(-2.0*alpha*(r1+r2));
           if (fabs(denominator) < pow(10,-6)){return 0;}
           else return drdr*exponential/sqrt(denominator);
        }

double probvalu(double x1,double x2, double x3, double x4, double x5, double x6, double N){
    double y1 = (pow((1-x1),(-1/(N -1)))-1);
    double y2 = (pow((1-x2),(-1/(N -1)))-1);
    double y3 = (pow((1-x3),(-1/(N -1)))-1);
    double y4 = (pow((1-x4),(-1/(N -1)))-1);
    double y5 = (pow((1-x5),(-1/(N -1)))-1);
    double y6 = (pow((1-x6),(-1/(N -1)))-1);

    //cout << x1 << setw(20) <<exp(-(y1))<<endl;
    return exp(-(y1+y2+y3+y4+y5+y6));
}


