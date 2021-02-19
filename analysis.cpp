//g++ -std=c++11 -o analysis analysis.cpp  

#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm>

using namespace std;

#define N 10
#define path "N10/"
#define bmin 7
#define bmax 30
#define step 1
#define ntherm 5000
#define nsim 10000

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

double trans_violation(int* lattice)
{
    double res = 0;

    for(int i=0; i<N; ++i)
    {
        for(int j=0; j<N; ++j)
        {
            if(i != j)
            {
                if(lattice[idx(i,j)] == 1)
                {
                    for(int k=0; k<N; ++k)
                    {
                        if(k != j)
                        {
                            if(lattice[idx(j,k)] == 1)
                            {
                                if(lattice[idx(i,k)] != 1) ++res;
                            }
                        }
                    }
                }
            }
        }
    }

    return res;
}

double antisymm_violation(int* lattice)
{
    double res = 0;

    for(int i=0; i<N; ++i)
    {
        for(int j=i+1; j<N; ++j)
        {
            if(lattice[idx(i,j)]*lattice[idx(j,i)] == 1) ++res;
        }
    }

    return res;
}
            




int main()
{
    cout << "Starting analysis with following parameters" << endl;
    cout << "N = " << N << endl;
    cout << "path = " << path << endl;
    cout << "beta min = " << bmin << endl;
    cout << "beta max = " << bmax << endl;
    cout << "beta step = " << step << endl;
    cout << "thermalization sweeps = " << ntherm << endl;
    cout << "total sweeps = " << nsim << endl;

    
    // out file
    ofstream out_h(path + string("H.txt"));
    ofstream out_antisymm(path + string("antisymm_violation.txt"));
    ofstream out_trans(path + string("trans_violation.txt"));


    double beta = bmin;

    while(beta < bmax)
    {
        // in file
        string name_h = path + string("H_b") + cc_to_name(beta) + ".txt";
        string name_l = path + string("L_b") + cc_to_name(beta) + ".txt";
        ifstream in_h(name_h);
        ifstream in_l(name_l);
        

        double h = 0;
        double trans = 0;
        double antisymm = 0;

        // read data
        for(int i=0; i<nsim; ++i)
        {
            // build lattice
            int* lattice = new int [N*N];

            for(int j=0; j<N*N; ++j)
                in_l >> lattice[j];

            double h_temp;
            in_h >> h_temp;

            if(i > ntherm)
            {
                trans += trans_violation(lattice);
                antisymm += antisymm_violation(lattice);
                h += h_temp;
            }

            delete [] lattice;
        }
                
        out_h << beta << " " << h/(nsim-ntherm) << endl;
        out_antisymm << beta << " " << antisymm/(nsim-ntherm) << endl;
        out_trans << beta << " " << trans/(nsim-ntherm) << endl;

        beta += step;
    
        in_h.close();
        in_l.close();
    }

    out_antisymm.close();
    out_trans.close();
}

