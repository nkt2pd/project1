#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <vector>
#include <ctime>
#include <fstream>
#include <iostream>
#include <chrono>
#include <random>
#include <set>
#include <cstring>
#include <algorithm>

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
double MC_sweep(double beta);
void Wolff_Cluster_Sim(double beta);
void Monte_Carlo_Sim(double beta);
void randomize();
double rand1();
void clear_Wolff_files();
void clear_Metrop_files();
double Wolff_Sweep(double beta);
inline int index_site(int x, int y);

mt19937 rng(chrono::steady_clock::now().time_since_epoch().count());

const int L = 30;       // chain size
const int D = 2;        // D-dimensional lattice
const int Ns = L * L;      // number of spins
const int N_nn1 = 4;

const double Jnn = -1;  // nearest-neighbor ferromagnetic coupling

// a simple data structure for square-lattice spin
class Site {
public:
    int idx;    // index 0, 1, 2, ... Ns-1
    int x, y;   // coordinates
    int Sz;     // Ising spins
    
    // idx of its four neighbors
    int idx_nn1[N_nn1];

};

class Cluster {
public:
    vector<int> sites;
    double size;
};

Site spin[Ns];

int main(int argc, char *argv[]) {
    if (argc != 2) {
        cout << "Please specify if you would like to use the Metropolis Algorithm" << endl;
        cout << "or the Wolff Cluster Algorithm" << endl;
        cout << "----------------------------------------------------" << endl;
        cout << "Usage: ./ising_2d_test <Metropolis/Wolff>";

        return 1;
    }
    
    const char Wolff[] = "Wolff";
    const char Metrop[] = "Metropolis";
    int mode = 0;

    if (strcmp(argv[1], Wolff) == 0) {
        mode = 1;
    } else if (strcmp(argv[1], Metrop) == 0) {
        mode = 2;
    } else {
        cout << "Please specify if you would like to use the Metropolis Algorithm" << endl;
        cout << "or the Wolff Cluster Algorithm" << endl;
        cout << "----------------------------------------------------" << endl;
        cout << "Usage: ./ising_2d <Metropolis/Wolff>";

        return 1;
    }
    
    //initialize chain of Ising Variables
    init_chain();
    
    set_coordinates();
    
    set_nn();
    randomize();
    
    if (mode == 1) {
        clear_Wolff_files();
        cout << "Simulating System w/ Wolff Cluster Algorithm" << endl;
        for(double T = 5; T>=0; T -= 0.1) {
        
            cout << "T = " << T << endl;
        
            Wolff_Cluster_Sim(1./T);
        }
    } else {
        clear_Metrop_files();
        cout << "Simulating System w/ Metropolis ALgorithm" << endl;
        for(double T = 5; T>=0; T -= 0.1) {
    
            cout << "T = " << T << endl;
        
            Monte_Carlo_Sim(1./T);
        }
    }    

    return 0;
}

void init_chain() {
    for (int i=0; i < Ns; i++) {
        spin[i].Sz = 1;
    }
}

inline int index_site(int x, int y) {
    // xm, ym take care of the periodic boundary condition
    int xm = mod(x, L);
    int ym = mod(y, L);
    return xm + L * ym;
}

void set_coordinates() {
    for(int x=0; x<L; x++)
    for(int y=0; y<L; y++) {
            
        int idx = index_site(x, y);
        spin[idx].idx = idx;
        spin[idx].x = x;
        spin[idx].y = y;
    }
}

void randomize() {
    for(int i=0; i<Ns; i++) {
        spin[i].Sz = (rand1() > 0.5) ? +1 : -1;
    } 
}

void set_nn() {
    int j;
    for(int i=0; i<Ns; i++) {
        
        int x = spin[i].x;
        int y = spin[i].y;
        
        j = index_site(x+1, y);
        spin[i].idx_nn1[0] = j;
        
        j = index_site(x-1, y);
        spin[i].idx_nn1[1] = j;
        
        j = index_site(x, y+1);
        spin[i].idx_nn1[2] = j;
        
        j = index_site(x, y-1);
        spin[i].idx_nn1[3] = j;
    }
}

int mod(int x, int m) {
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
        int m = spin[k].idx_nn1[l];
        Sz_nn += spin[m].Sz;
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

double MC_sweep(double beta) {
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
            int j = spin[i].idx_nn1[k];
            sum += spin[i].Sz * spin[j].Sz;
        }
    }
    return 0.5 * Jnn * sum;
}

double Wolff_Sweep(double beta) {
    
    double prob = 1. - exp(2.0*beta*Jnn);
    double r = rand1();
    double R = 0;
    r *= L*L;

    if (r == 900.) r--;

    Cluster C;
    Cluster F_old;
    Cluster F_new;

    F_old.sites.clear();
    C.sites.clear();
    C.size = 0;
    C.sites.push_back((int)r);
    F_old.sites.push_back((int)r);

    int in_cluster[Ns] = {0};
    in_cluster[(int)r] = 1;

    while(!F_old.sites.empty()) {
        F_new.sites.clear();

        for (int i = 0; i < F_old.sites.size(); i++) {

            for (int j = 0; j < N_nn1; j++) {

                int k = spin[F_old.sites[i]].idx_nn1[j];

                if (spin[k].Sz == spin[F_old.sites[i]].Sz && in_cluster[k] == 0) {

                    if (rand1() < prob) {

                        F_new.sites.push_back(spin[k].idx);
                        C.sites.push_back(spin[k].idx);
                        in_cluster[k] = 1;
                        C.size += 1.;

                        //cout << C.size;
                    }
                }
            }
        }
            F_old.sites = F_new.sites;
    }
    
    //Once cluster is made, flip all the spins
    for (int k = 0; k < C.sites.size(); k++) {
        spin[C.sites[k]].Sz *= -1;
    }

    return C.size;
}

void Wolff_Cluster_Sim(double beta) {

    int thermalize = 5000;
    int nsweep = (int) (0.5 * pow(beta, -3.5));
    int ndata = 600000;

    //Run 5000 sweeps of the system to achieve equilibrium
    double total_members = 0;
    for (int i = 0; i < thermalize; i++) {
        total_members += Wolff_Sweep(beta);
    }
    cout << "Average Cluster Size = " << total_members/((double) thermalize) << endl;

    //now start to collect data
    double E1 = 0, E2 = 0;
    double M1 = 0, M2 = 0, M4 = 0;
    double avg_size = 0;

    int system_check = 100000;

    for (int n = 0; n < ndata; n++) {

        if(n % system_check == 0) cout << "n = " << n << endl;

        total_members = 0;
        for (int r = 0; r < nsweep; r++) {
            total_members += Wolff_Sweep(beta);
        }
        total_members /= ((double) nsweep);
        avg_size = (n * avg_size + total_members) / (n + 1.);

        double e = get_energy();
        E1 = (n * E1 + e) / (n + 1.);
        E2 = (n * E2 + e*e) / (n + 1.);

        double mag = get_magnetization();

        M1 = (n * M1 + mag) / (n + 1.);
        M2 = (n * M2 + pow(mag, 2)) / (n + 1.);
        M4 = (n * M4 + pow(mag, 4)) / (n + 1.);
    }

    ofstream Wolff_energy;
    Wolff_energy.open("Wolff_energy.dat", fstream::app);
    ofstream Wolff_heat;
    Wolff_heat.open("Wolff_heat.dat", fstream::app);

    Wolff_energy << 1./beta << ", " << E1/((double) Ns) << "\n";
    Wolff_heat << 1./beta << ", " << pow(beta, 2) * (E2 - E1*E1) / ((double) Ns) << '\n';

    Wolff_energy.close();
    Wolff_heat.close();
}

void Monte_Carlo_Sim(double beta) {

    int thermalize = 5000;
    int nsweep = 20;
    int ndata = 600000;

    //Run 5000 sweeps of the system to achieve equilibrium
    double accepted = 0;
    for (int i = 0; i < thermalize; i++) {
        accepted += MC_sweep(beta);
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

            accepted += MC_sweep(beta);

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

    ofstream ising_energy_2d;
    ising_energy_2d.open("ising_energy_2d.dat", fstream::app);
    ofstream ising_heat_2d;
    ising_heat_2d.open("ising_heat_2d.dat", fstream::app);

    ising_energy_2d << 1./beta << ", " << E1/((double) Ns) << "\n";
    ising_heat_2d << 1./beta << ", " << pow(beta, 2) * (E2 - E1*E1) / ((double) Ns) << '\n';

    ising_energy_2d.close();
    ising_heat_2d.close();
}

double rand1() {
    
    std::uniform_real_distribution<double> dis(0.0, 1.0);
    double rand1 = dis(rng);
    
    return rand1;
}

void clear_Metrop_files() {
    ofstream clear1;
    clear1.open("ising_energy_2d.dat", fstream::trunc);
    ofstream clear2;
    clear2.open("ising_heat_2d.dat", fstream::trunc);

    clear1.close();
    clear2.close();
}

void clear_Wolff_files() {
    ofstream clear1;
    clear1.open("Wolff_energy.dat", fstream::trunc);
    ofstream clear2;
    clear2.open("Wolff_heat.dat", fstream::trunc);

    clear1.close();
    clear2.close();
}