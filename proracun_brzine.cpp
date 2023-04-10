#include <iostream>
#include <istream>
#include <fstream>
#include <vector>
#include <cmath>

class Motor
{
    double U1 = 0;  // napon
    double f = 0;   // frekvencija
    double R1 = 0;  // aktivni otpor statora
    double R2 = 0;  // aktivni otpor rotora
    double X1 = 0;  // induktivni otpor statora
    double X2 = 0;  // induktivni otpor rotora
    double J = 0;   // moment inercije rotora
    double n_n = 0; // nominalna brzina motora
    double n_s = 0; // sinhrona brzina motora

public:
    friend std::istream &operator>>(std::istream &in, Motor &motor);
    friend double proracun_momenta(const Motor &motor, double s);
    friend void proracun_brzine(const double &moment, std::vector<std::pair<double, double>> &brzine, const Motor &motor, const double &dt);
    friend void proracun_brzine_u_funkciji_vremena(const std::vector<std::pair<double, double>> &momenti, std::vector<std::pair<double, double>> &brzine, const Motor &motor);
};

/* UNOS PODATAKA */

std::istream &operator>>(std::istream &in, Motor &motor)
{
    std::cout << "Napon: ";
    in >> motor.U1;
    std::cout << "Frekvencija: ";
    in >> motor.f;
    std::cout << "Sinhrona brzina motora: ";
    in >> motor.n_s;
    std::cout << "Nominalna brzina motora: ";
    in >> motor.n_n;
    std::cout << "Aktivni otpor statora: ";
    in >> motor.R1;
    std::cout << "Aktivni otpor rotora: ";
    in >> motor.R2;
    std::cout << "Induktivni otpor statora: ";
    in >> motor.X1;
    std::cout << "Induktivni otpor rotora: ";
    in >> motor.X2;
    std::cout << "Moment inercije rotora: ";
    in >> motor.J;
    return in;
}

void unos_specifikacija(Motor &motor)
{
    std::cout << std::endl;
    std::cout << "Unesite specifikacije motora:\n"
              << std::endl;
    std::cin >> motor;
}

/* PRORACUNI */

double proracun_momenta(const Motor &motor, double s)
{
    return ((3 * motor.R2 * pow(motor.U1, 2)) / (2 * 3.14159 * motor.f * s * (pow((motor.R1 + (motor.R2 / s)), 2) + pow((motor.X1 + motor.X2), 2))));
}

void proracun_momenata_u_funkciji_klizanja(const Motor &motor, std::vector<std::pair<double, double>> &momenti)
{
    for (double i = 1; i >= 0; i -= 0.01)
    {
        momenti.push_back(std::make_pair(proracun_momenta(motor, i), i));
    }
}

void proracun_brzine(const double &moment, std::vector<std::pair<double, double>> &brzine, const Motor &motor, const double &dt)
{
    if (brzine.empty())
    { 
        brzine.push_back(std::make_pair(((moment / motor.J) * dt), dt));
    }
    else
    {
        brzine.push_back(std::make_pair((((moment / motor.J) * dt) + brzine[brzine.size() - 1].first), dt));
    }
}

void proracun_brzine_u_funkciji_vremena(const std::vector<std::pair<double, double>> &momenti, std::vector<std::pair<double, double>> &brzine, const Motor &motor)
{
    double s_n = (motor.n_s - motor.n_n) / motor.n_s;
    double omega = 0;
    double n = 0;
    double s = 0;
    double dt = 0;
    double moment = 0;
    int brojac = 0;
    while (brojac < momenti.size())
    {
        n = 0.1047198 * omega;
        s = (motor.n_s - n) / motor.n_s;
        int i = 0;
        while (momenti[i].second > s && momenti[i].second >= s_n)
        {
            ++i;
        }
        if (i != 0)
        {
            --i;
        }
        moment = momenti[i].first;
        proracun_brzine(moment, brzine, motor, dt);
        omega = brzine[brzine.size() - 1].first;
        ++brojac;
        dt += 0.1;
    }
}

/* ISPISI U FILE */

void ispis_momenata_u_file(const std::vector<std::pair<double, double>> &momenti)
{
    std::ofstream out("momenti.txt");
    for (int i = 0; i < momenti.size(); ++i)
    {
        out << momenti[i].second << " " << momenti[i].first << std::endl;
    }
}

void ispis_brzina_u_file(const std::vector<std::pair<double, double>> &brzine)
{
    std::ofstream out("brzine.txt");
    for (int i = 0; i < brzine.size(); ++i)
    {
        out << brzine[i].second << " " << brzine[i].first << std::endl;
    }
}

int main(int argc, char const *argv[])
{
    Motor motor;
    std::vector<std::pair<double, double>> momenti;
    std::vector<std::pair<double, double>> brzine;

    unos_specifikacija(motor);
    proracun_momenata_u_funkciji_klizanja(motor, momenti);
    proracun_brzine_u_funkciji_vremena(momenti, brzine, motor);

    ispis_momenata_u_file(momenti);
    ispis_brzina_u_file(brzine);

    return 0;
}
