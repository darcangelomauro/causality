//g++ -std=c++11 -o main main.cpp -I/home/pmxmd10/gsl/include -L/home/pmxmd10/gsl/lib -lgsl -lopenblas -lm 

#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <gsl/gsl_rng.h>
#include <ctime>
#include <cmath>

using namespace std;

#define N 20

string cc_to_name(const double& g2)
{
    ostringstream osg2;
    osg2 << g2;
    string sg2 = osg2.str();
    replace(sg2.begin(), sg2.end(), '.', 'd');
    replace(sg2.begin(), sg2.end(), '-', 'n');
    
    return sg2;
}

int idx(const int& i, const int& j)
{
    return N*i + j;
}

double H(int* lattice, const double& beta)
{
    // quadratic part
    double H2=0;

    for(int i=0; i<N; ++i)
    {
        for(int j=i+1; j<N; ++j)
            H2 += lattice[idx(i,j)]*lattice[idx(j,i)];
    }


    // quartic part
    double H4=0;

    for(int i=0; i<N; ++i)
    {
        for(int j=0; j<N; ++j)
        {
            if(i!=j)
            {
                for(int k=0; k<N; ++k)
                {
                    if( (i!=k) && (j!=k) )
                        H4 += lattice[idx(i,j)]*lattice[idx(k,j)]*lattice[idx(k,i)]*(1+lattice[idx(k,i)]);
                }
            }
        }
    }


    return (H2 - beta*H4)/(N*N);
}



int main()
{
    // initialize random number generator
    gsl_rng* engine = gsl_rng_alloc(gsl_rng_ranlxd2);
    gsl_rng_set(engine, time(NULL));
    
    
    // build lattice
    int* lattice = new int [N*N];


    double beta = 3;

    while(beta < 15)
    {
        // initialize lattice
        for(int i=0; i<N*N; ++i)
            lattice[i] = 0;
        
        double prev_H;

        // out file
        string name_h = "H_b" + cc_to_name(beta) + ".txt";
        string name_l = "L_b" + cc_to_name(beta) + ".txt";
        ofstream out_h(name_h);
        ofstream out_l(name_l);
        
        // simulation
        for(int i=0; i<10000; ++i)
        {
            // sweep the lattice
            for(int j=0; j<N*N; ++j)
            {
                // store previous value of spin (j)
                int prev_spin = lattice[j];

                // find initial energy
                double en_i;
                if(i==0 && j==0)
                    en_i = H(lattice, beta);
                else
                    en_i = prev_H;

                // randomly change value of spin (j)
                if(gsl_rng_uniform(engine) < 0.5)
                    lattice[j] += 1;
                else
                    lattice[j] += 2;

                if(lattice[j] == 2)
                    lattice[j] = -1;
                if(lattice[j] == 3)
                    lattice[j] = 0;

                // find final energy
                double en_f = H(lattice, beta);

                // accept with probability min[1, exp(en_i-en_f)]
                if(en_f > en_i)
                {
                    double e = exp(en_i-en_f);
                    if(gsl_rng_uniform(engine) > e)
                    {
                        // restore old configuration
                        prev_H = en_i;
                        lattice[j] = prev_spin;
                    }
                    else
                        prev_H = en_f;
                }
                else
                    prev_H = en_f;
            }

            // print H and lattice
            out_h << prev_H << endl;
            for(int k=0; k<N*N; ++k)
                out_l << lattice[k] << " ";
            out_l << endl;
        }

        out_l.close();
        out_h.close();
        
        beta += 0.1;
    }
    
        
        
    delete [] lattice;
}
