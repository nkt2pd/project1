#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <vector>
#include <iostream>
#include <fstream>
#include <string>

void scale_binder4(std::vector<double> temps_30, std::vector<double> temps_60, std::vector<double> temps_120);
void scale_m(std::vector<double> temps_30, std::vector<double> temps_60, std::vector<double> temps_120);
void scale_suscept(std::vector<double> temps_30, std::vector<double> temps_60, std::vector<double> temps_120);
void scale_heat(std::vector<double> temps_30, std::vector<double> temps_60, std::vector<double> temps_120);
void get_temps(std::vector<double>& temps_30, std::vector<double>& temps_60, std::vector<double>& temps_120);

int main() {

    std::vector<double> temps_30;
    std::vector<double> temps_60;
    std::vector<double> temps_120;

    get_temps(temps_30, temps_60, temps_120);

    scale_binder4(temps_30, temps_60, temps_120);
    scale_m(temps_30, temps_60, temps_120);
    scale_suscept(temps_30, temps_60, temps_120);
    scale_heat(temps_30, temps_60, temps_120);

    return 0;
}

void get_temps(std::vector<double>& temps_30, std::vector<double>& temps_60, std::vector<double>& temps_120){
    for (double T = 2.8; T>= 2; T-=.01) {
        temps_30.push_back((T-2.269)*30);
        temps_60.push_back((T-2.269)*60);
        temps_120.push_back((T-2.269)*120);
    }
}

void scale_binder4(std::vector<double> temps_30, std::vector<double> temps_60, std::vector<double> temps_120) {

    std::string read;
    std::ofstream Wolff_binder4_fss;
    std::ifstream Wolff_binder4_30;
    std::ifstream Wolff_binder4_60;
    std::ifstream Wolff_binder4_120;

    Wolff_binder4_fss.open("Wolff_binder4_fss.dat");
    Wolff_binder4_30.open("Wolff_binder4_30.dat");
    Wolff_binder4_60.open("Wolff_binder4_60.dat");
    Wolff_binder4_120.open("Wolff_binder4_120.dat");

    std::vector<double> binder4_30;
    std::vector<double> binder4_60;
    std::vector<double> binder4_120;

    int get = 0;
    int counter = 0;

    while (Wolff_binder4_30 >> read) {
        if (get == 1) {
            binder4_30[counter] = std::stod(read);
            get = 0;
            counter++;
        }
        get++;
    }

    counter = 0;
    get = 0;

    while (Wolff_binder4_60 >> read) {
        if (get == 1) {
            binder4_60[counter] = std::stod(read);
            get = 0;
            counter++;
        }
        get++;
    }

    counter = 0;
    get = 0;

    while (Wolff_binder4_120 >> read) {
        if (get == 1) {
            binder4_120[counter] = std::stod(read);
            get = 0;
            counter++;
        }
        get++;
    }

    for (int i = 0; i < binder4_30.size(); i++) {
        Wolff_binder4_fss << temps_30[i] << ", " << binder4_30[i] << ", " << temps_60[i] << ", " << binder4_60[i] << ", " 
                          << temps_120[i] << ", " << binder4_120[i] << std::endl;
    }

    Wolff_binder4_fss.close();
    Wolff_binder4_30.close();
    Wolff_binder4_60.close();
    Wolff_binder4_120.close();

}

void scale_m(std::vector<double> temps_30, std::vector<double> temps_60, std::vector<double> temps_120) {

}

void scale_suscept(std::vector<double> temps_30, std::vector<double> temps_60, std::vector<double> temps_120) {

}

void scale_heat(std::vector<double> temps_30, std::vector<double> temps_60, std::vector<double> temps_120) {

}