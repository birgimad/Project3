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
        clock_t start,finish;

        //double x[N]={};



        srand(time(NULL));


        for (int i=0;i<N;i++){
            for (int j=0;j<6;j++){
                x[j] = ((double) rand() / (RAND_MAX)); //random numbers generated in the interval(0,1)

                //cout<<x[j]<<endl;

            }

            double pi = 3.1415;
            double y1 = -log(1-x[0]);
            double y2 = -log(1-x[1]);
            double y3 = x[2]*pi;
            double y4 = x[3]*pi;
            double y5 = x[4]*2*pi;
            double y6 = x[5]*2*pi;

            start=clock();
            fx = func(y1,y2,y3,y4,y5,y6);
            px = 1/(4*pi*pi*pi*pi);
            Montec += fx/px;
            Montesqr += (fx/px)*(fx/px);
            finish = clock();
            double t = ((finish-start)/CLOCKS_PER_SEC);

       }

       double MC = Montec/((double) N );
       double Montesqr1= Montesqr/((double) N);
       double variance= Montesqr1-MC*MC;
       double volume = pow(length,6);
       cout << " Integral = " << MC << endl;
       cout << "sigma = " << sqrt(variance/ ((double) N)) << endl;

        ofstream myfile ("Monte_Carlo10.txt");
        if (myfile.is_open())
             {
                 myfile << "Computed integral value using Improved Monte carlo:" <<endl;
                 myfile << "N = " << N << endl;
                 myfile << "Integral value:" << setw(20) << int_GaussianLaguerre << endl;


         //ofstream myfile
                // if (myfile.is_open())

                     myfile << "Computed CPU time for: " <<endl;
                    // myfile << "n = " << n << endl;
                     myfile << "TIME:"<< setw(5) << t <<"sec" << endl;
                 }



      return 0;
    }
        double func(double r1, double r2, double theta1, double theta2, double phi1, double phi2)
        {
           double alpha = 2.;

           double drdr = r1*r1*r2*r2*(-sin(theta1))*(-sin(theta2));
           double denominator = sqrt(r1*r1+r2*r2-2*r1*r2*(cos(theta1)*cos(theta2)+sin(theta1)*sin(theta2)*cos(phi1-phi2)));
           //double exponential = exp(-2.0*alpha*(r1+r2));
           if (fabs(denominator) < pow(10,-6)){return 0;}
           else return drdr/(pow(2*alpha,5)*denominator); // exponential removed
        }




