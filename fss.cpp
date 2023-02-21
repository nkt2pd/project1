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

    int get = 1;
    int counter = 0;
    double val = 0;

    while (Wolff_binder4_30 >> read) {
        if (get%2 == 0) {
            val = std::stod(read);
            binder4_30.push_back(val);
            counter++;
        }
        get++;
    }

    counter = 0;
    get = 1;

    while (Wolff_binder4_60 >> read) {
        if (get%2 == 0) {
            val = std::stod(read);
            binder4_60.push_back(val);
            counter++;
        }
        get++;
    }

    counter = 0;
    get = 1;

    while (Wolff_binder4_120 >> read) {
        if (get%2 == 0) {
            val = std::stod(read);
            binder4_120.push_back(val);
            counter++;
        }
        get++;
    }

    for (int i = 0; i <= 80; i++) {
        Wolff_binder4_fss << temps_30[i] << ", " << binder4_30[i] << ", " << temps_60[i] << ", " << binder4_60[i] << ", " 
                          << temps_120[i] << ", " << binder4_120[i] << std::endl;
    }

    Wolff_binder4_fss.close();
    Wolff_binder4_30.close();
    Wolff_binder4_60.close();
    Wolff_binder4_120.close();

}

void scale_m(std::vector<double> temps_30, std::vector<double> temps_60, std::vector<double> temps_120) {

    std::string read;
    std::ofstream Wolff_m_fss;
    std::ifstream Wolff_m_30;
    std::ifstream Wolff_m_60;
    std::ifstream Wolff_m_120;

    Wolff_m_fss.open("Wolff_m_fss.dat");
    Wolff_m_30.open("Wolff_m_30.dat");
    Wolff_m_60.open("Wolff_m_60.dat");
    Wolff_m_120.open("Wolff_m_120.dat");

    std::vector<double> m_30;
    std::vector<double> m_60;
    std::vector<double> m_120;

    int get = 1;
    int counter = 0;
    double val = 0;

    while (Wolff_m_30 >> read) {
        if (get%2 == 0) {
            val = std::stod(read);
            m_30.push_back(val);
            counter++;
        }
        get++;
    }

    counter = 0;
    get = 1;

    while (Wolff_m_60 >> read) {
        if (get%2 == 0) {
            val = std::stod(read);
            m_60.push_back(val);
            counter++;
        }
        get++;
    }

    counter = 0;
    get = 1;

    while (Wolff_m_120 >> read) {
        if (get%2 == 0) {
            val = std::stod(read);
            m_120.push_back(val);
            counter++;
        }
        get++;
    }

    for (int i = 0; i <= 80; i++) {
        Wolff_m_fss << temps_30[i] << ", " << m_30[i] * pow(30, .125) << ", " << temps_60[i] << ", " << m_60[i] * pow(60, .125) << ", " 
                          << temps_120[i] << ", " << m_120[i] * pow(120, .125) << std::endl;
    }

    Wolff_m_fss.close();
    Wolff_m_30.close();
    Wolff_m_60.close();
    Wolff_m_120.close();

}

void scale_suscept(std::vector<double> temps_30, std::vector<double> temps_60, std::vector<double> temps_120) {

    std::string read;
    std::ofstream Wolff_suscept_fss;
    std::ifstream Wolff_suscept_30;
    std::ifstream Wolff_suscept_60;
    std::ifstream Wolff_suscept_120;

    Wolff_suscept_fss.open("Wolff_susceptibility_fss.dat");
    Wolff_suscept_30.open("Wolff_susceptibility_30.dat");
    Wolff_suscept_60.open("Wolff_susceptibility_60.dat");
    Wolff_suscept_120.open("Wolff_susceptibility_120.dat");

    std::vector<double> suscept_30;
    std::vector<double> suscept_60;
    std::vector<double> suscept_120;

    int get = 1;
    int counter = 0;
    double val = 0;

    while (Wolff_suscept_30 >> read) {
        if (get%2 == 0) {
            val = std::stod(read);
            suscept_30.push_back(val);
            counter++;
        }
        get++;
    }

    counter = 0;
    get = 1;

    while (Wolff_suscept_60 >> read) {
        if (get%2 == 0) {
            val = std::stod(read);
            suscept_60.push_back(val);
            counter++;
        }
        get++;
    }

    counter = 0;
    get = 1;

    while (Wolff_suscept_120 >> read) {
        if (get%2 == 0) {
            val = std::stod(read);
            suscept_120.push_back(val);
            counter++;
        }
        get++;
    }

    for (int i = 0; i <= 80; i++) {
        Wolff_suscept_fss << temps_30[i] << ", " << suscept_30[i] * pow(30, -1.75) << ", " << temps_60[i] << ", " << suscept_60[i] * pow(60, -1.75) << ", " 
                          << temps_120[i] << ", " << suscept_120[i] * pow(120, -1.75) << std::endl;
    }

    Wolff_suscept_fss.close();
    Wolff_suscept_30.close();
    Wolff_suscept_60.close();
    Wolff_suscept_120.close();

}

void scale_heat(std::vector<double> temps_30, std::vector<double> temps_60, std::vector<double> temps_120) {

    std::string read;
    std::ofstream Wolff_heat_fss;
    std::ifstream Wolff_heat_30;
    std::ifstream Wolff_heat_60;
    std::ifstream Wolff_heat_120;

    Wolff_heat_fss.open("Wolff_heat_fss.dat");
    Wolff_heat_30.open("Wolff_heat_30.dat");
    Wolff_heat_60.open("Wolff_heat_60.dat");
    Wolff_heat_120.open("Wolff_heat_120.dat");

    std::vector<double> heat_30;
    std::vector<double> heat_60;
    std::vector<double> heat_120;

    int get = 1;
    int counter = 0;
    double val = 0;

    while (Wolff_heat_30 >> read) {
        if (get%2 == 0) {
            val = std::stod(read);
            heat_30.push_back(val);
            counter++;
        }
        get++;
    }

    counter = 0;
    get = 1;

    while (Wolff_heat_60 >> read) {
        if (get%2 == 0) {
            val = std::stod(read);
            heat_60.push_back(val);
            counter++;
        }
        get++;
    }

    counter = 0;
    get = 1;

    while (Wolff_heat_120 >> read) {
        if (get%2 == 0) {
            val = std::stod(read);
            heat_120.push_back(val);
            counter++;
        }
        get++;
    }

    for (int i = 0; i <= 80; i++) {
        Wolff_heat_fss << temps_30[i] << ", " << heat_30[i] / log(30) << ", " << temps_60[i] << ", " << heat_60[i] / log(60) << ", " 
                          << temps_120[i] << ", " << heat_120[i] / log(120) << std::endl;
    }

    Wolff_heat_fss.close();
    Wolff_heat_30.close();
    Wolff_heat_60.close();
    Wolff_heat_120.close();

}