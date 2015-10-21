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
        double jacobidet = pow((2*length),6);

        //double x[N]={};



        srand(time(NULL));


        for (int i=0;i<N;i++){
            for (int j=0;j<6;j++){
                x[j] = - length +2*length*((double) rand() / (RAND_MAX)); //random numbers generated in the interval(0,1)

                //cout<<x[j]<<endl;

            }

            fx = func(x[0],x[1],x[2],x[3],x[4],x[5]);
            Montec += fx;
            Montesqr += fx*fx;

       }

       double MC = Montec/((double) N );
       double Montesqr1= Montesqr/((double) N);
       double variance= Montesqr1-MC*MC;
       double volume = pow(length,6);
       cout << " Integral = " << jacobidet*MC << endl;
       cout << "sigma = " << volume*sqrt(variance/ ((double) N)) << endl;



      return 0;
    }
        double func(double x1, double y1, double z1, double x2, double y2, double z2)
        {
           double alpha = 2.;

           double r1=sqrt(x1*x1+y1*y1+z1*z1);
           double r2=sqrt(x2*x2+y2*y2+z2*z2);
           double deno=sqrt(pow((x1-x2),2)+pow((y1-y2),2)+pow((z1-z2),2));

           if(deno <pow(10.,-6.)) { return 0;}
           else return exp(-2*alpha*(r1+r2))/deno;
        }
