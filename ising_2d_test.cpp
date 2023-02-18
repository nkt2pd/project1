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

const int N_nn1 = 4;

class Site {
public:
    int idx;    // index 0, 1, 2, ... Ns-1
    int x, y;   // coordinates
    int Sz;     // Ising spins
    int Ns;
    int L;
    // idx of its four neighbors
    Site *nn1[N_nn1];

};

int update_site(int k, double beta, Site *spin, int *Ns, int *L);
double get_energy(Site *spin, int *Ns, int *L);
int get_magnetization(Site *spin, int *Ns, int *L);
void init_chain(Site *spin, int *Ns, int *L);
void set_coordinates(Site *spin, int *Ns, int *L);
void set_nn(Site *spin, int *Ns, int *L);
int mod_LR(int x, int m);
int mod_UD(int x, int m);
double MC_sweep(double beta, Site *spin, int *Ns, int *L);
void Wolff_Cluster_Sim(double beta, int mode, Site *spin, int *Ns, int *L, int pick_file);
void Monte_Carlo_Sim(double beta, Site *spin, int *Ns, int *L);
void randomize(Site *spin, int *Ns, int *L);
double rand1();
void clear_Wolff_files();
void clear_Metrop_files();
double Wolff_Sweep(double beta, Site *spin, int *Ns, int *L);
void print_Wolff(int i, int mode, double E1, double E2, double M1, double M2, double M4, double beta, int *Ns);

mt19937 rng(chrono::steady_clock::now().time_since_epoch().count());

const int L1 = 30;       // chain size
const int L2 = 60;
const int L3 = 120;
const int L4 = 240;
const int L5 = 360;
const int D = 2;        // D-dimensional lattice
const int Ns1 = L1 * L1;      // number of spins
const int Ns2 = L2 * L2;
const int Ns3 = L3 * L3;
const int Ns4 = L4 * L4;
const int Ns5 = L5 * L5;

const double Jnn = -1;  // nearest-neighbor ferromagnetic coupling

// a simple data structure for square-lattice spin

class Cluster {
public:
    vector<int> sites;
    double size;
};

Site spin1[Ns1];
Site spin2[Ns2];
Site spin3[Ns3];
Site spin4[Ns4];
Site spin5[Ns5];

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
        cout << "Usage: ./ising_2d_test <Metropolis/Wolff>";

        return 1;
    }

    spin1[0].Ns = Ns1;
    spin1[0].L = L1;
    spin2[0].Ns = Ns2;
    spin2[0].L = L2;
    spin3[0].Ns = Ns3;
    spin3[0].L = L3;
    spin4[0].Ns = Ns4;
    spin4[0].L = L4;
    spin5[0].Ns = Ns5;
    spin5[0].L = L5;

    if (mode == 1) {
        clear_Wolff_files();

        Site *spin;
        int *Ns;
        int *L;
        for (int i = 0; i < 5; i++) {
            if (i == 0) {
                spin = spin1;
                Ns = &spin1[0].Ns;
                L = &spin1[0].L;
                cout << "Lattice Length: 30" << endl;
            } else if (i == 1) {
                spin = spin2;
                Ns = &spin2[0].Ns;
                L = &spin2[0].L;
                cout << "Lattice Length: 60" << endl;
            } else if (i == 2) {
                spin = spin3;
                Ns = &spin3[0].Ns;
                L = &spin3[0].L;
                cout << "Lattice Length: 120" << endl;
            } else if (i == 3) {
                spin = spin4;
                Ns = &spin4[0].Ns;
                L = &spin4[0].L;
                cout << "Lattice Length: 240" << endl;
            } else {
                spin = spin5;
                Ns = &spin5[0].Ns;
                L = &spin5[0].L;
                cout << "Lattice Length: 360" << endl;
            }
            
            init_chain(spin, Ns, L);
            
            set_coordinates(spin, Ns, L);
            set_nn(spin, Ns, L);

            randomize(spin, Ns, L);

            double del_T = 0.1;
            int mode_dat = 0;

            cout << "Simulating System w/ Wolff Cluster Algorithm" << endl;
            
            for(double T = 2.8; T>=2; T -= del_T) {
                
                //collect finite size scaling data for relevant temp range
                if (T <= 2.81 && T > 2) {
                    del_T = .01;
                    mode_dat = 1;
                } else {
                    del_T = .1;
                    mode_dat = 0;
                }

                cout << "T = " << T << endl;
            
                Wolff_Cluster_Sim(1./T, mode_dat, spin, Ns, L, i);
            }
        }
    } else {
        Site *spin = spin1;
        int *Ns = &spin1[0].Ns;
        int *L = &spin1[0].L;

        init_chain(spin, Ns, L);

        set_coordinates(spin1, Ns, L);
        set_nn(spin1, Ns, L);

        randomize(spin1, Ns, L);

        clear_Metrop_files();
        cout << "Simulating System w/ Metropolis ALgorithm" << endl;
        for(double T = 5; T>=0; T -= 0.1) {
    
            cout << "T = " << T << endl;
        
            Monte_Carlo_Sim(1./T, spin, Ns, L);
        }
    }    

    return 0;
}
    
void init_chain(Site *spin, int *Ns, int *L) {
    for (int i=0; i < *Ns; i++) {
        spin[i].Sz = 1;
    }
}

void set_coordinates(Site *spin, int *Ns, int *L) {
    int y = 0;
    int x = 0;
    for(int i=0; i<spin[0].Ns; i++) {
            spin[i].idx = i;
            spin[i].x = x;
            spin[i].y = y;
            x++;
            if (x%*L == 0) {
                y++;
                x = 0;
            }
        }

    /* for (int j = 0; j < L; j++) {
        for (int k = 0; k < L; k++) {
            cout << "(" << spin[k + (j*L)].x << ", " << spin[k + (j*L)].y << ") ";
        }
        cout << endl;
    } */
}

void randomize(Site *spin, int *Ns, int *L) {
    for(int i=0; i<*Ns; i++) {
        spin[i].Sz = (rand1() > 0.5) ? +1 : -1;
    } 
}

void set_nn(Site *spin, int *Ns, int *L) {
    int k;
    for(int i=0; i<spin[0].Ns; i++) {
        
        k = mod_LR(spin[i].x + 1, *L);
        k += (*L * spin[i].y);
        spin[i].nn1[0] = &spin[k];
        
        k = mod_LR(spin[i].x - 1, *L);
        k += (*L * spin[i].y);
        spin[i].nn1[1] = &spin[k];
        
        k = mod_UD(i-*L, *Ns);
        spin[i].nn1[2] = &spin[k];

        k = mod_UD(i+*L, *Ns);
        spin[i].nn1[3] = &spin[k];

    }

    /*for(int m = 0; m < L; m++) {
        for(int n = 0; n < L; n++) {
            cout << "(" << spin[n + m*L].nn1[0]->idx << ", " << spin[n + m*L].nn1[1]->idx << ", " <<
                    spin[n + m*L].nn1[2]->idx << ", " << spin[n + m*L].nn1[3]->idx << ") ";
        }
        cout << endl;
    }*/
    
}

int mod_LR(int x, int m) { //implements periodic boundary condition in x-direction
    if (x>=0 && x<m) {
        return x;
    } else if (x<0) {
        return m-1-mod_LR(-1-x,m);
    } else {
        return x%m;
    }
}

int mod_UD(int x, int m) { //implements periodic boundary condition in y-direction
    if (x>=0 && x<m)
        return x;
    else if (x<0)
        return m-1-mod_UD(-1-x,m);
    else
        return x%m;
}

//test a new state, accept or reject the state based on parameters. 
int update_site(int k, double beta, Site *spin, int *Ns, int *L) {

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

double MC_sweep(double beta, Site *spin, int *Ns, int *L) {
    int hits = 0;
    for(int i=0; i<*Ns; i++) {
        hits += update_site(i, beta, spin, Ns, L);
    }
    return ((double) hits)/((double) *Ns);      // success rate
}

int get_magnetization(Site *spin, int *Ns, int *L) {
    int sum = 0;
    for (int i = 0; i < *Ns; i++) {
        sum += spin[i].Sz;
    }

    return sum;
}

double get_energy(Site *spin, int *Ns, int *L) {
    int sum = 0;
    for(int i=0; i<*Ns; i++) {
        for(int k=0; k<N_nn1; k++) {
            sum += spin[i].Sz * spin[i].nn1[k]->Sz;
        }
    }
    return 0.5 * Jnn * sum;
}

double Wolff_Sweep(double beta, Site *spin, int *Ns, int *L) {
    double prob = 1 - exp(2*beta*Jnn);
    double r = rand1();
    double R = 0;
    r *= *L* (*L);

    if (r == 900.) r--;

    Cluster C;
    Cluster F_old;
    Cluster F_new;

    F_old.sites.clear();
    C.sites.clear();
    C.size = 0;
    C.sites.push_back((int)r);
    F_old.sites.push_back((int)r);

    int in_cluster[*Ns] = {0};
    in_cluster[(int)r] = 1;

    while(!F_old.sites.empty()) {
        F_new.sites.clear();

        for (int i = 0; i < F_old.sites.size(); i++) {

            for (int j = 0; j < N_nn1; j++) {

                int k = spin[F_old.sites[i]].nn1[j]->idx;

                if (spin[F_old.sites[i]].nn1[j]->Sz == spin[F_old.sites[i]].Sz && in_cluster[k] == 0) {

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

void Wolff_Cluster_Sim(double beta, int mode, Site *spin, int *Ns, int *L, int pick_file) {

    int thermalize = 5000;
    int nsweep = (int) (0.5 * pow(beta, -3.5));;
    int ndata = 600000;

    //Run 5000 sweeps of the system to achieve equilibrium
    double total_members = 0;
    for (int i = 0; i < thermalize; i++) {
        total_members += Wolff_Sweep(beta, spin, Ns, L);
    }
    cout << "Average Cluster Size = " << total_members/((double) thermalize) << endl;

    //now start to collect data
    double E1 = 0, E2 = 0;
    double M1 = 0, M2 = 0, M4 = 0, abs_M = 0;
    double avg_size = 0;

    int system_check = 100000;

    for (int n = 0; n < ndata; n++) {

        if(n % system_check == 0) cout << "n = " << n << endl;

        total_members = 0;

        for (int r = 0; r < nsweep; r++) {
            total_members += Wolff_Sweep(beta, spin, Ns, L);
        }

        total_members /= ((double) nsweep);
        avg_size = (n * avg_size + total_members) / (n + 1.);

        double e = get_energy(spin, Ns, L);
        E1 = (n * E1 + e) / (n + 1.);
        E2 = (n * E2 + e*e) / (n + 1.);

        double mag = abs(get_magnetization(spin, Ns, L));
        M1 = (n * M1 + mag) / (n + 1.);
        M2 = (n * M2 + pow(mag, 2)) / (n + 1.);
        M4 = (n * M4 + pow(mag, 4)) / (n + 1.);

    }

    print_Wolff(pick_file, mode, E1, E2, M1, M2, M4, beta, Ns);

}

void Monte_Carlo_Sim(double beta, Site *spin, int *Ns, int *L) {

    int thermalize = 5000;
    int nsweep = 20;
    int ndata = 600000;

    //Run 5000 sweeps of the system to achieve equilibrium
    double accepted = 0;
    for (int i = 0; i < thermalize; i++) {
        accepted += MC_sweep(beta, spin, Ns, L);
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

            accepted += MC_sweep(beta, spin, Ns, L);

        }
        accepted /= ((double) nsweep);
        avg_accept = (n * avg_accept + accepted) / (n + 1.);

        double e = get_energy(spin, Ns, L);

        E1 = (n * E1 + e) / (n + 1.);
        E2 = (n * E2 + e*e) / (n + 1.);

        double mag = get_magnetization(spin, Ns, L);

        M1 = (n * M1 + mag) / (n + 1.);
        M2 = (n * M2 + pow(mag, 2)) / (n + 1.);
        M4 = (n * M4 + pow(mag, 4)) / (n + 1.);
    }

    ofstream ising_energy_2d;
    ising_energy_2d.open("ising_energy_2d.dat", fstream::app);
    ofstream ising_heat_2d;
    ising_heat_2d.open("ising_heat_2d.dat", fstream::app);

    ising_energy_2d << 1./beta << ", " << E1/((double) *Ns) << "\n";
    ising_heat_2d << 1./beta << ", " << pow(beta, 2) * (E2 - E1*E1) / ((double) *Ns) << '\n';

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
    ofstream clear1, clear2, clear3, clear4, clear5,
             clear6, clear7, clear8, clear9, clear10,
             clear11, clear12, clear13, clear14, clear15,
             clear16, clear17, clear18, clear19, clear20,
             clear21, clear22, clear23, clear24, clear25;
    
    clear1.open("Wolff_energy_30.dat", fstream::trunc);
    clear2.open("Wolff_heat_30.dat", fstream::trunc);
    clear3.open("Wolff_binder4_30.dat", fstream::trunc);
    clear4.open("Wolff_m_30.dat", fstream::trunc);
    clear5.open("Wolff_susceptibility_30.dat", fstream::trunc);
    clear6.open("Wolff_energy_60.dat", fstream::trunc);
    clear7.open("Wolff_heat_60.dat", fstream::trunc);
    clear8.open("Wolff_binder4_60.dat", fstream::trunc);
    clear9.open("Wolff_m_60.dat", fstream::trunc);
    clear10.open("Wolff_susceptibility_60.dat", fstream::trunc);
    clear11.open("Wolff_energy_120.dat", fstream::trunc);
    clear12.open("Wolff_heat_120.dat", fstream::trunc);
    clear13.open("Wolff_binder4_120.dat", fstream::trunc);
    clear14.open("Wolff_m_120.dat", fstream::trunc);
    clear15.open("Wolff_susceptibility_120.dat", fstream::trunc);
    clear16.open("Wolff_energy_240.dat", fstream::trunc);
    clear17.open("Wolff_heat_240.dat", fstream::trunc);
    clear18.open("Wolff_binder4_240.dat", fstream::trunc);
    clear19.open("Wolff_m_240.dat", fstream::trunc);
    clear20.open("Wolff_susceptibility_240.dat", fstream::trunc);
    clear21.open("Wolff_energy_360.dat", fstream::trunc);
    clear22.open("Wolff_heat_360.dat", fstream::trunc);
    clear23.open("Wolff_binder4_360.dat", fstream::trunc);
    clear24.open("Wolff_m_360.dat", fstream::trunc);
    clear25.open("Wolff_susceptibility_360.dat", fstream::trunc);

    clear1.close();
    clear2.close();
    clear3.close();
    clear4.close();
    clear5.close();
    clear6.close();
    clear7.close();
    clear8.close();
    clear9.close();
    clear10.close();
    clear11.close();
    clear12.close();
    clear13.close();
    clear14.close();
    clear15.close();
    clear16.close();
    clear17.close();
    clear18.close();
    clear19.close();
    clear20.close();
    clear21.close();
    clear22.close();
    clear23.close();
    clear24.close();
    clear25.close();
}

void print_Wolff(int i, int mode, double E1, double E2, double M1, double M2, double M4, double beta, int *Ns) {
    if (i == 0) {
        ofstream Wolff_energy_30;
        Wolff_energy_30.open("Wolff_energy_30.dat", fstream::app);
        ofstream Wolff_heat_30;
        Wolff_heat_30.open("Wolff_heat_30.dat", fstream::app);

        Wolff_energy_30 << 1./beta << ", " << E1/((double) *Ns) << endl;
        Wolff_heat_30 << 1./beta << ", " << pow(beta, 2) * (E2 - E1*E1) / ((double) *Ns) << endl;

        Wolff_energy_30.close();
        Wolff_heat_30.close();

        if (mode == 1) {
            ofstream Wolff_m_30;
            Wolff_m_30.open("Wolff_m_30.dat", fstream::app);
            ofstream Wolff_susceptibility_30;
            Wolff_susceptibility_30.open("Wolff_susceptibility_30.dat", fstream::app);
            ofstream Wolff_binder4_30;
            Wolff_binder4_30.open("Wolff_binder4_30.dat", fstream::app);

            Wolff_m_30 << 1./beta << ", " << M1/((double) *Ns) << endl;
            Wolff_susceptibility_30 << 1./beta << ", " << beta*(M2 - M1*M1)/((double) *Ns) << endl;
            Wolff_binder4_30 << 1./beta << ", " << 1. - (M4/(3 * pow(M2, 2))) << endl;

            Wolff_m_30.close();
            Wolff_susceptibility_30.close();
            Wolff_binder4_30.close();
        }
    } else if (i == 1) {
        ofstream Wolff_energy_60;
        Wolff_energy_60.open("Wolff_energy_60.dat", fstream::app);
        ofstream Wolff_heat_60;
        Wolff_heat_60.open("Wolff_heat_60.dat", fstream::app);

        Wolff_energy_60 << 1./beta << ", " << E1/((double) *Ns) << endl;
        Wolff_heat_60 << 1./beta << ", " << pow(beta, 2) * (E2 - E1*E1) / ((double) *Ns) << endl;

        Wolff_energy_60.close();
        Wolff_heat_60.close();

        if (mode == 1) {
            ofstream Wolff_m_60;
            Wolff_m_60.open("Wolff_m_60.dat", fstream::app);
            ofstream Wolff_susceptibility_60;
            Wolff_susceptibility_60.open("Wolff_susceptibility_60.dat", fstream::app);
            ofstream Wolff_binder4_60;
            Wolff_binder4_60.open("Wolff_binder4_60.dat", fstream::app);

            Wolff_m_60 << 1./beta << ", " << M1/((double) *Ns) << endl;
            Wolff_susceptibility_60 << 1./beta << ", " << beta*(M2 - M1*M1)/((double) *Ns) << endl;
            Wolff_binder4_60 << 1./beta << ", " << 1. - (M4/(3 * pow(M2, 2))) << endl;

            Wolff_m_60.close();
            Wolff_susceptibility_60.close();
            Wolff_binder4_60.close();
        }
    } else if (i == 2) {
        ofstream Wolff_energy_120;
        Wolff_energy_120.open("Wolff_energy_120.dat", fstream::app);
        ofstream Wolff_heat_120;
        Wolff_heat_120.open("Wolff_heat_120.dat", fstream::app);

        Wolff_energy_120 << 1./beta << ", " << E1/((double) *Ns) << endl;
        Wolff_heat_120 << 1./beta << ", " << pow(beta, 2) * (E2 - E1*E1) / ((double) *Ns) << endl;

        Wolff_energy_120.close();
        Wolff_heat_120.close();

        if (mode == 1) {
            ofstream Wolff_m_120;
            Wolff_m_120.open("Wolff_m_120.dat", fstream::app);
            ofstream Wolff_susceptibility_120;
            Wolff_susceptibility_120.open("Wolff_susceptibility_120.dat", fstream::app);
            ofstream Wolff_binder4_120;
            Wolff_binder4_120.open("Wolff_binder4_120.dat", fstream::app);

            Wolff_m_120 << 1./beta << ", " << M1/((double) *Ns) << endl;
            Wolff_susceptibility_120 << 1./beta << ", " << beta*(M2 - M1*M1)/((double) *Ns) << endl;
            Wolff_binder4_120 << 1./beta << ", " << 1. - (M4/(3 * pow(M2, 2))) << endl;

            Wolff_m_120.close();
            Wolff_susceptibility_120.close();
            Wolff_binder4_120.close();
        }
    } else if (i == 3) {
        ofstream Wolff_energy_240;
        Wolff_energy_240.open("Wolff_energy_240.dat", fstream::app);
        ofstream Wolff_heat_240;
        Wolff_heat_240.open("Wolff_heat_240.dat", fstream::app);

        Wolff_energy_240 << 1./beta << ", " << E1/((double) *Ns) << endl;
        Wolff_heat_240 << 1./beta << ", " << pow(beta, 2) * (E2 - E1*E1) / ((double) *Ns) << endl;

        Wolff_energy_240.close();
        Wolff_heat_240.close();

        if (mode == 1) {
            ofstream Wolff_m_240;
            Wolff_m_240.open("Wolff_m_240.dat", fstream::app);
            ofstream Wolff_susceptibility_240;
            Wolff_susceptibility_240.open("Wolff_susceptibility_240.dat", fstream::app);
            ofstream Wolff_binder4_240;
            Wolff_binder4_240.open("Wolff_binder4_240.dat", fstream::app);

            Wolff_m_240 << 1./beta << ", " << M1/((double) *Ns) << endl;
            Wolff_susceptibility_240 << 1./beta << ", " << beta*(M2 - M1*M1)/((double) *Ns) << endl;
            Wolff_binder4_240 << 1./beta << ", " << 1. - (M4/(3 * pow(M2, 2))) << endl;

            Wolff_m_240.close();
            Wolff_susceptibility_240.close();
            Wolff_binder4_240.close();
        }
    } else {
        ofstream Wolff_energy_360;
        Wolff_energy_360.open("Wolff_energy_360.dat", fstream::app);
        ofstream Wolff_heat_360;
        Wolff_heat_360.open("Wolff_heat_360.dat", fstream::app);

        Wolff_energy_360 << 1./beta << ", " << E1/((double) *Ns) << endl;
        Wolff_heat_360 << 1./beta << ", " << pow(beta, 2) * (E2 - E1*E1) / ((double) *Ns) << endl;

        Wolff_energy_360.close();
        Wolff_heat_360.close();

        if (mode == 1) {
            ofstream Wolff_m_360;
            Wolff_m_360.open("Wolff_m_360.dat", fstream::app);
            ofstream Wolff_susceptibility_360;
            Wolff_susceptibility_360.open("Wolff_susceptibility_360.dat", fstream::app);
            ofstream Wolff_binder4_360;
            Wolff_binder4_360.open("Wolff_binder4_360.dat", fstream::app);

            Wolff_m_360 << 1./beta << ", " << M1/((double) *Ns) << endl;
            Wolff_susceptibility_360 << 1./beta << ", " << beta*(M2 - M1*M1)/((double) *Ns) << endl;
            Wolff_binder4_360 << 1./beta << ", " << 1. - (M4/(3 * pow(M2, 2))) << endl;

            Wolff_m_360.close();
            Wolff_susceptibility_360.close();
            Wolff_binder4_360.close();
        }
    }
}