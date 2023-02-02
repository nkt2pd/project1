#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <vector>
#include <ctime>
#include <fstream>
#include <iostream>
#include <chrono>
#include <random>

#define chain_size 50
#define trials_per_temp 100000

using namespace std;

int update_site(int k, double beta);
double get_energy();
int get_magnetization();
void init_chain();
void set_coordinates();
void set_nn();
int mod(int x, int m);
double sweep(double beta);
void Monte_Carlo_Sim(double beta);
void randomize();
void init_rand();
double rand1();

mt19937 rng(chrono::steady_clock::now().time_since_epoch().count());

const int L = 50;       // chain size
const int Ns = L;       // number of spins
const int N_nn1 = 2;

const double Jnn = -1;  // nearest-neighbor ferromagnetic coupling

// a simple data structure for square-lattice spin
class Site {
public:
    int idx;    // index 0, 1, 2, ... Ns-1
    int x;   // coordinates
    int Sz;     // Ising spins
    
    // idx of its four neighbors
    Site *nn1[N_nn1];
};

Site spin[Ns];

int main() {

    //initialize chain of Ising Variables
    init_chain();

    set_coordinates();
    set_nn();

    randomize();

    for(double T = 5; T>=0; T -= 0.1) {
    
        cout << "T = " << T << endl;
        
        Monte_Carlo_Sim(1./T);
    
    }

    return 0;
}

void init_chain() {
    for (int i=0; i < Ns; i++) {
        spin[i].Sz = 1;
    }
}

void set_coordinates() {
    
    for(int x=0; x<L; x++) {
            spin[x].idx = x;
            spin[x].x = x;
        }
}

void set_nn() {
    int j;
    for(int i=0; i<Ns; i++) {
        
        j = mod(i+1, L);
        spin[i].nn1[0] = &spin[j];
        
        j = mod(i-1, L);
        spin[i].nn1[1] = &spin[j];
        
    }
}

void randomize() {
    for(int i=0; i<Ns; i++) spin[i].Sz = (rand1() > 0.5) ? +1 : -1;
}

int mod(int x, int m) { //implements periodic boundary condition
    if (x>=0 && x<m)
        return x;
    else if (x<0)
        return m-1-mod(-1-x,m);
    else
        return x%m;
}

//test a new state, accept or reject the state based on parameters. 
int update_site(int k, double beta) {

    int Sz_nn = 0;
    double delE = 0;

    for(int l=0; l<N_nn1; l++) {
        Sz_nn += spin[k].nn1[l]->Sz;
    }

    double Hamiltonian = Jnn * Sz_nn;

    //change in energy when spin is flipped
    delE = 2. * -spin[k].Sz * Hamiltonian;

    double r = rand1();

    if (delE == 0) {
        //if there is no change in the energy, coin flip to see if accepted
        if (r < 0.5) {
            spin[k].Sz *= -1;
            return 1;
        } else return 0;
    } else {
        if (r < exp(-delE * beta)) {
            spin[k].Sz *= -1;
            return 1;
        } else {
        return 0;
        }
    }
}

double sweep(double beta) {
    int hits = 0;
    for(int i=0; i<Ns; i++) {
        hits += update_site(i, beta);
    }
    return ((double) hits)/((double) Ns);      // success rate
}

int get_magnetization() {
    int sum;
    for (int i = 0; i < Ns; i++) {
        sum += spin[i].Sz;
    }
    return sum;
}

double get_energy() {
    int sum = 0;
    for(int i=0; i<Ns; i++) {
        for(int k=0; k<N_nn1; k++) {
            sum += spin[i].Sz * spin[i].nn1[k]->Sz;
        }
    }
    return 0.5 * Jnn * sum;
}

void Monte_Carlo_Sim(double beta) {

    int thermalize = 5000;
    int nsweep = 20;
    int ndata = 600000;

    //Run 5000 sweeps of the system to achieve equilibrium
    double accepted = 0;
    for (int i = 0; i < thermalize; i++) {
        accepted += sweep(beta);
    }
    cout << "spin update rate = " << accepted/((double) thermalize) << endl;

    //now start to collect data
    double E1 = 0, E2 = 0;
    double M1 = 0, M2 = 0, M4 = 0;
    double avg_accept = 0;

    int system_check = 100000;

    for (int n = 0; n < ndata; n++) {

        if(n % system_check == 0) cout << "n = " << n << endl;

        accepted = 0;
        for (int r = 0; r < nsweep; r++) {

            accepted += sweep(beta);

        }
        accepted /= ((double) nsweep);
        avg_accept = (n * avg_accept + accepted) / (n + 1.);

        double e = get_energy();
        E1 = (n * E1 + e) / (n + 1.);
        E2 = (n * E2 + e*e) / (n + 1.);

        double mag = get_magnetization();

        M1 = (n * M1 + mag) / (n + 1.);
        M2 = (n * M2 + pow(mag, 2)) / (n + 1.);
        M4 = (n * M4 + pow(mag, 4)) / (n + 1.);
    }

    ofstream ising_energy;
    ising_energy.open("ising_energy.dat", fstream::app);
    ofstream ising_heat;
    ising_heat.open("ising_heat.dat", fstream::app);

    ising_energy << 1./beta << ", " << E1/((double) Ns) << "\n";
    ising_heat << 1./beta << ", " << pow(beta, 2) * (E2 - E1*E1) / ((double) Ns) << '\n';

    ising_energy.close();
    ising_heat.close();
}

void init_rand() {
    
}

double rand1() {
    
    std::uniform_real_distribution<double> dis(0.0, 1.0);
    double rand1 = dis(rng);
    
    return rand1;
}