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
        double *x = new double [N];
        double Montec = 0;
        double Montesqr = 0;
        double fx;

        //double x[N]={};



        srand(time(NULL));


        for (int i=0;i<N;i++){
            for (int j=0;j<N;j++){
                for (int k=0;k<N;k++){
                    for (int l=0;l<N;l++){
                        for (int m=0;m<N;m++){
                            for (int n=0;n<N;n++){
                                x[i] = ((double) rand() / (RAND_MAX)); //random numbers generated in the interval(0,1)

                                fx = func(x[i],x[j],x[k],x[l],x[m],x[n]);
                                           Montec += fx;
                                           Montesqr += fx*fx;
                            }
                        }
                    }
                }
            }


       }

                  double MC = Montec/((double) N );
                        double Montesqr1= Montesqr/((double) N);
                        //double variance=Montesqr1-MC*MC;
                        cout << " Integral = " << MC << endl;


      return 0;
    }
        double func(double x1, double y1, double z1, double x2, double y2, double z2)
        {
           double alpha = 2.;

           double exp1=-2*alpha*sqrt(x1*x1+y1*y1+z1*z1);
           double exp2=-2*alpha*sqrt(x2*x2+y2*y2+z2*z2);
           double deno=sqrt(pow((x1-x2),2)+pow((y1-y2),2)+pow((z1-z2),2));

           if(deno <pow(10.,-6.)) { return 0;}
          else return exp(exp1+exp2)/deno;
        }


