#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>

#include "ODEsolvers.hpp"
#include "LaneEmdenLibrary.hpp"

using namespace std;

// costanti fisiche in unita' SI
const double G = 6.67430e-11;
const double hbar = 1.054571817e-34;
const double c = 2.99792458e8;
const double m_p = 1.67262192369e-27;
const double M_sun = 1.98847e30;
const double R_sun = 6.957e8;

int main()
{
    const double dxi = 1e-3; // passo di integrazione
    const double eps = 1e-6; // punto iniziale
    const double xi_max = 50.0; // limite massimo per ricerca xi_s

    std::cout << setprecision(12);
    std::cout << "\n\n----------------------------\n";
    std::cout << "PARAMETRI PER L'INTEGRAZIONE\n";
    std::cout << "----------------------------\n\n";

    std::cout << "passo di integrazione: dxi = " << dxi << ",\n";
    std::cout << "traslazione dall'origine: eps = " << eps << ",\n";
    std::cout << "raggio massimo: xi_max = " << xi_max << ".\n\n";

    std::cout << "---------\n";
    std::cout << "SOLUZIONI\n";
    std::cout << "---------\n\n";

    // caso n = 0 (laneemden_n0.dat)
    {
        double n = 0.0;
        g_n = n;

        ofstream out("laneemden_n0.dat");
        out << setprecision(12);

        double eps2 = eps * eps;
        double eps4 = eps2 * eps2;
        double xi = eps;
        
        double Y[2];
        Y[0] = 1.0 - eps2 / 6.0 + (n * eps4) / 120.0;
        Y[1] = -eps / 3.0 + (n * eps2 * eps) / 30.0;

        // valore iniziali
        out << 0.0 << "\t\t\t" << 1.0 << "\t\t\t" << 1.0 << "\n";
        out << xi  << "\t\t\t" << Y[0] << "\t\t\t" << pow(Y[0], n) << "\n";

        // integrazione
        double xi_prev = xi, th_prev = Y[0];
        const int neq = 2;
        const int N = (int)((xi_max - eps) / dxi);

        for (int i = 0; i < N; i++) {
            
            RK4(xi, Y, LaneEmdenRHS, dxi, neq);
            xi = eps + (i + 1) * dxi;

            if (th_prev > 0.0 && Y[0] <= 0.0) {
                double t = th_prev / (th_prev - Y[0]);
                double xi_s = xi_prev + t * (xi - xi_prev);

                // densita' = 0 per convenzione, all'interno vale rho_c
                out << xi_s << "\t\t" << 0.0 << "\t\t" << 0.0 << "\n";
                break;
            }

            out << xi << "\t\t" << Y[0] << "\t\t" << ((Y[0] > 0.0) ? pow(Y[0], n) : 0.0) << "\n";

            xi_prev = xi;
            th_prev = Y[0];

        }

        out.close();

        // soluzione analitica
        {
            ofstream outa("laneemden_n0_anal.dat");
            outa << setprecision(12);

            int N = (int)(xi_max / dxi);

            for (int i = 0; i <= N; i++) {
                double xi = i * dxi;
                double theta_a = 1.0 - (xi * xi) / 6.0;

                // la soluzione analitica e' considerata fino al primo zero xi_s = sqrt(6)
                if (theta_a <= 0.0) {
                    outa << sqrt(6.0) << "\t\t" << 0.0 << "\n";
                    break;
                }

                outa << xi << "\t\t" << theta_a << "\n";
            }

            outa.close();
        }

        double xi_s = SolveLaneEmden(n, dxi, xi_max, eps);

        std::cout << "--- CASO n = 0 ---\n";
        if (g_found_surface) {
            std::cout << "zero (soluzione numerica): xi_s = " << xi_s << ",\n";
            std::cout << "zero (soluzione analitica): sqrt(6) = " << sqrt(6.0) << ",\n";
            std::cout << "derivata calcolata allo zero: theta'(xi_s) = " << g_thetaPrime_s << ",\n";
            std::cout << "errore massimo = " << g_max_err_n0 << ".\n\n";
        }

    }

    // caso n = 1 (laneemden_n1.dat)
    {
        double n = 1.0;
        g_n = n;

        ofstream out("laneemden_n1.dat");
        out << setprecision(12);

        double eps2 = eps * eps;
        double eps4 = eps2 * eps2;
        double xi = eps;

        double Y[2];
        Y[0] = 1.0 - eps2 / 6.0 + (n * eps4) / 120.0;
        Y[1] = -eps / 3.0 + (n * eps2 * eps) / 30.0;

        // valori iniziali 
        out << 0.0 << "\t\t\t" << 1.0 << "\t\t\t" << 1.0 << "\n";
        out << xi  << "\t\t\t" << Y[0] << "\t\t\t" << pow(Y[0], n) << "\n";

        // integrazione
        double xi_prev = xi, th_prev = Y[0];
        const int neq = 2;
        const int N = (int)((xi_max - eps) / dxi);

        for (int i = 0; i < N; i++) {
            
            RK4(xi, Y, LaneEmdenRHS, dxi, neq);
            xi = eps + (i + 1) * dxi;

            if (th_prev > 0.0 && Y[0] <= 0.0) {
                double t = th_prev / (th_prev - Y[0]);
                double xi_s = xi_prev + t * (xi - xi_prev);
                
                out << xi_s << "\t\t" << 0.0 << "\t\t" << 0.0 << "\n";
                break;
            }

            out << xi << "\t\t" << Y[0] << "\t\t" << ((Y[0] > 0.0) ? pow(Y[0], n) : 0.0) << "\n";

            xi_prev = xi;
            th_prev = Y[0];

        }

        out.close();

        // soluzione analitica 
        {
            ofstream outa("laneemden_n1_anal.dat");
            outa << setprecision(12);

            int N = (int)(xi_max / dxi);

            for (int i = 0; i <= N; i++) {
                double xi = i * dxi;
                double theta_a;

                // nell'origine si usa sin(xi)/xi -> 1
                if (xi == 0.0)
                    theta_a = 1.0;
                else
                    theta_a = sin(xi) / xi;

                // la soluzione e' considerata fino al primo zero xi_s = pi
                if (theta_a <= 0.0) {
                    outa << M_PI << "\t\t" << 0.0 << "\n";
                    break;
                }

                outa << xi << "\t\t" << theta_a << "\n";
            }

            outa.close();
        }

        double xi_s = SolveLaneEmden(n, dxi, xi_max, eps);

        std::cout << "--- CASO n = 1 ---\n";
        if (g_found_surface) {
            std::cout << "zero (soluzione numerica): xi_s = " << xi_s << ",\n";
            std::cout << " zero (soluzione analitica): pi = " << M_PI << ",\n";
            std::cout << "derivata calcolata allo zero: theta'(xi_s) = " << g_thetaPrime_s << ",\n";
            std::cout << "errore massimo = " << g_max_err_n1 << ".\n\n";
        }

    }

    // caso n = 3/2 (laneemden_n1p5.dat)
    {
        double n = 1.5;
        g_n = n;

        ofstream out("laneemden_n1p5.dat");
        out << setprecision(12);

        double eps2 = eps * eps;
        double eps4 = eps2 * eps2;
        double xi = eps;

        double Y[2];
        Y[0] = 1.0 - eps2 / 6.0 + (n * eps4) / 120.0;
        Y[1] = -eps / 3.0 + (n * eps2 * eps) / 30.0;

        // valori iniziali 
        out << 0.0 << "\t\t\t" << 1.0 << "\t\t\t" << 1.0 << "\n";
        out << xi  << "\t\t\t" << Y[0] << "\t\t\t" << pow(Y[0], n) << "\n";

        // integrazione 
        double xi_prev = xi, th_prev = Y[0];
        const int neq = 2;
        const int N = (int)((xi_max - eps) / dxi);

        for (int i = 0; i < N; i++) {
            
            RK4(xi, Y, LaneEmdenRHS, dxi, neq);
            xi = eps + (i + 1) * dxi;

            if (th_prev > 0.0 && Y[0] <= 0.0) {
                double t = th_prev / (th_prev - Y[0]);
                double xi_s = xi_prev + t * (xi - xi_prev);
                
                out << xi_s << "\t\t" << 0.0 << "\t\t" << 0.0 << "\n";
                break;
            }

            out << xi << "\t\t" << Y[0] << "\t\t" << ((Y[0] > 0.0) ? pow(Y[0], n) : 0.0) << "\n";

            xi_prev = xi;
            th_prev = Y[0];

        }

        out.close();

        double xi_s = SolveLaneEmden(n, dxi, xi_max, eps);

        std::cout << "--- CASO n = 3/2 ---\n";
        if (g_found_surface) {
            std::cout << "zero (soluzione numerica): xi_s = " << xi_s << ",\n";
            std::cout << "derivata calcolata allo zero: theta'(xi_s) = " << g_thetaPrime_s << ".\n\n";
        }

    }

    // caso n = 3 (laneemden_n3.dat)
    {
        double n = 3.0;
        g_n = n;

        ofstream out("laneemden_n3.dat");
        out << setprecision(12);

        double eps2 = eps * eps;
        double eps4 = eps2 * eps2;
        double xi = eps;

        double Y[2];
        Y[0] = 1.0 - eps2 / 6.0 + (n * eps4) / 120.0;
        Y[1] = -eps / 3.0 + (n * eps2 * eps) / 30.0;

        // valori iniziali
        out << 0.0 << "\t\t\t" << 1.0 << "\t\t\t" << 1.0 << "\n";
        out << xi  << "\t\t\t" << Y[0] << "\t\t\t" << pow(Y[0], n) << "\n";

        // integrazione
        double xi_prev = xi, th_prev = Y[0];
        const int neq = 2;
        const int N = (int)((xi_max - eps) / dxi);

        for (int i = 0; i < N; i++) {
            
            RK4(xi, Y, LaneEmdenRHS, dxi, neq);
            xi = eps + (i + 1) * dxi;

            if (th_prev > 0.0 && Y[0] <= 0.0) {
                double t = th_prev / (th_prev - Y[0]);
                double xi_s = xi_prev + t * (xi - xi_prev);
                
                out << xi_s << "\t\t" << 0.0 << "\t\t" << 0.0 << "\n";
                break;
            }

            out << xi << "\t\t" << Y[0] << "\t\t" << ((Y[0] > 0.0) ? pow(Y[0], n) : 0.0) << "\n";

            xi_prev = xi;
            th_prev = Y[0];

        }

        out.close();

        double xi_s = SolveLaneEmden(n, dxi, xi_max, eps);

        std::cout << "--- CASO n = 3 ---\n";
        if (g_found_surface) {
            std::cout << "zero (soluzione numerica): xi_s = " << xi_s << ",\n";
            std::cout << "derivata calcolata allo zero: theta'(xi_s) = " << g_thetaPrime_s << ".\n\n";
        }

        if (g_found_surface) {
            const double K = (pow(3.0, 1.0/3.0) * pow(M_PI, 2.0/3.0)) / (pow(2.0, 4.0/3.0) * 4.0) * (hbar * c) / pow(m_p, 4.0/3.0); // costante eq. di stato politropica
            const double omega = - (xi_s * xi_s) * g_thetaPrime_s; // costante adimensionale di massa

            // densita' centrali [kg/m^3] per verifica indipendenza di M da rho_c
            const int nrho = 5;
            const double rho_c[nrho] = {1e8, 1e9, 1e10, 1e11, 1e12};

            double M[nrho];
            double mass_solar[nrho];

            double R[nrho];
            double radius_solar[nrho];

            // calcolo di masse e raggi in SI e conversione in masse e raggi solari
            for (int i = 0; i < nrho; i++) {
                M[i] = PolytropeMass(3.0, K, rho_c[i], xi_s, g_thetaPrime_s);
                R[i] = PolytropeRadius(3.0, K, rho_c[i], xi_s);

                mass_solar[i] = M[i] / M_sun;
                radius_solar[i] = R[i] / R_sun;
            }

            // ricerca della massa minima e massima
            double mass_min = mass_solar[0];
            double mass_max = mass_solar[0];

            for (int i = 1; i < nrho; i++) {
                if (mass_solar[i] < mass_min)
                    mass_min = mass_solar[i];

                if (mass_solar[i] > mass_max)
                    mass_max = mass_solar[i];
            }

            // massima variazione percentuale rispetto alla massa minima
            double mass_variation = 0.0;

            if (mass_min != 0.0)
                mass_variation = 100.0 * (mass_max - mass_min) / fabs(mass_min);

            std::cout << "- MASSA DI CHANDRASEKHAR -\n";
            std::cout << setprecision(12);
            std::cout << "costante eq. L-E per n = 3: omega = -xi_s^2*theta'(xi_s) = " << omega << ",\n";
            std::cout << "costante eq. di stato politropica: K = " << K << ".\n\n";

            std::cout << "rho_c [kg/m^3] \t\t M [M_sun] \t\t R [R_sun]\n";
            std::cout << "---------------------------------------------------------------\n";

            for (int i = 0; i < nrho; i++) {
                std::cout << scientific << setprecision(6) << rho_c[i] << " \t\t " << fixed << setprecision(12) << mass_solar[i] << " \t " << radius_solar[i] << "\n";
            }

            std::cout << "\nMassima variazione relativa della massa = " << scientific << setprecision(6) << mass_variation << " %.\n";
            std::cout << defaultfloat << "\n";
        }

    }

    // caso n = 4 (laneemden_n4.dat)
    {
        double n = 4.0;
        g_n = n;

        ofstream out("laneemden_n4.dat");
        out << setprecision(12);

        double eps2 = eps * eps;
        double eps4 = eps2 * eps2;
        double xi = eps;

        double Y[2];
        Y[0] = 1.0 - eps2 / 6.0 + (n * eps4) / 120.0;
        Y[1] = -eps / 3.0 + (n * eps2 * eps) / 30.0;

        // valori iniziali
        out << 0.0 << "\t\t\t" << 1.0 << "\t\t\t" << 1.0 << "\n";
        out << xi  << "\t\t\t" << Y[0] << "\t\t\t" << pow(Y[0], n) << "\n";

        // integrazione
        double xi_prev = xi, th_prev = Y[0];
        const int neq = 2;
        const int N = (int)((xi_max - eps) / dxi);

        for (int i = 0; i < N; i++) {
            
            RK4(xi, Y, LaneEmdenRHS, dxi, neq);
            xi = eps + (i + 1) * dxi;

            if (th_prev > 0.0 && Y[0] <= 0.0) {
                double t = th_prev / (th_prev - Y[0]);
                double xi_s = xi_prev + t * (xi - xi_prev);

                out << xi_s << "\t\t" << 0.0 << "\t\t" << 0.0 << "\n";
                break;
            }

            out << xi << "\t\t" << Y[0] << "\t\t" << ((Y[0] > 0.0) ? pow(Y[0], n) : 0.0) << "\n";

            xi_prev = xi;
            th_prev = Y[0];

        }

        out.close();

        double xi_s = SolveLaneEmden(n, dxi, xi_max, eps);

        std::cout << "--- CASO n = 4 ---\n";
        if (g_found_surface) {
            std::cout << "zero (soluzione numerica): xi_s = " << xi_s << ",\n";
            std::cout << "derivata calcolata allo zero: theta'(xi_s) = " << g_thetaPrime_s << ".\n\n";
        }

    }

    // caso n = 5 (laneemden_n5.dat)
    {
        double n = 5.0;
        g_n = n;

        ofstream out("laneemden_n5.dat");
        out << setprecision(12);

        double eps2 = eps * eps;
        double eps4 = eps2 * eps2;
        double xi = eps;

        double Y[2];
        Y[0] = 1.0 - eps2 / 6.0 + (n * eps4) / 120.0;
        Y[1] = -eps / 3.0 + (n * eps2 * eps) / 30.0;

        // valori iniziali
        out << 0.0 << "\t\t\t" << 1.0 << "\t\t\t" << 1.0 << "\n";
        out << xi  << "\t\t\t" << Y[0] << "\t\t\t" << pow(Y[0], n) << "\n";

        // integrazione
        double xi_prev = xi, th_prev = Y[0];
        const int neq = 2;
        const int N = (int)((xi_max - eps) / dxi);

        for (int i = 0; i < N; i++) {
            
            RK4(xi, Y, LaneEmdenRHS, dxi, neq);
            xi = eps + (i + 1) * dxi;

            if (th_prev > 0.0 && Y[0] <= 0.0) {
                double t = th_prev / (th_prev - Y[0]);
                double xi_s = xi_prev + t * (xi - xi_prev);
                
                out << xi_s << "\t\t" << 0.0 << "\t\t" << 0.0 << "\n";
                break;
            }

            out << xi << "\t\t" << Y[0] << "\t\t" << ((Y[0] > 0.0) ? pow(Y[0], n) : 0.0) << "\n";

            xi_prev = xi;
            th_prev = Y[0];

        }

        out.close();

        // soluzione analitica su una griglia fino a xi_max
        {
            ofstream outa("laneemden_n5_anal.dat");
            outa << setprecision(12);

            int N = (int)(xi_max / dxi);

            for (int i = 0; i <= N; i++) {
                
                double xi = i * dxi;
                double theta_a = 1.0 / sqrt(1.0 + (xi * xi) / 3.0);
                
                outa << xi << "\t\t" << theta_a << "\n";
            }

            outa.close();
        }

        SolveLaneEmden(n, dxi, xi_max, eps);

        std::cout << "--- CASO n = 5 ---\n";
        std::cout << "soluzione analitica: theta_5(xi) = 1/sqrt(1 + xi^2/3),\n";
        std::cout << "errore massimo (fino a xi_max) = " << g_max_err_n5 << ",\n";
        if (!g_found_surface) {
            std::cout << "superficie non trovata entro xi_max = " << xi_max << " (raggio infinito).\n\n";
        }

    }

    // calcolo coefficiente -xi_s^2*theta'(xi_s) per n tra 0 e 4.99
    ofstream coeff("laneemden_mass_coefficient.dat");
        coeff << setprecision(12);
        coeff << "# n\t\txi_s\t\t\t-xi_s^2 theta'(xi_s)\n";

    {
        // intervallo esteso per indivduare superficie per n -> 5
        const double xi_max = 1800.0;

        // n da 0 a 4 con passo di 0.5
        for (int i = 0; i <= 8; i++)
        {
            double n = 0.5*i;
            double xi_s = SolveLaneEmden(n, dxi, xi_max, eps);

            if (g_found_surface)
                coeff << n << "\t\t" << xi_s << "\t\t" << -xi_s*xi_s*g_thetaPrime_s << "\n";
        }

        // n da 4.1 a 4.9 con passo di 0.1
        for (int i = 0; i <= 8; i++)
        {
            double n = 4.1 + 0.1*i;
            double xi_s = SolveLaneEmden(n, dxi, xi_max, eps);

            if (g_found_surface)
                coeff << n << "\t\t" << xi_s << "\t\t" << -xi_s*xi_s*g_thetaPrime_s << "\n";
        }

        // n da 4.1 a 4.9 con passo di 0.01
        for (int i = 0; i <= 9; i++)
        {
            double n = 4.9 + 0.01*i;
            double xi_s = SolveLaneEmden(n, dxi, xi_max, eps);

            if (g_found_surface)
                coeff << n << "\t\t" << xi_s << "\t\t" << -xi_s*xi_s*g_thetaPrime_s << "\n";
        }
    }

    coeff.close();
    
    // test convergenza RK4 per n = 1
    ofstream conv("RK4_convergence.dat");

    conv << setprecision(12);

    // estremi dell'intervallo usato nel test
    const double eps_conv = 1e-6;
    const double xi_max_conv = 3.0;

    // numero di intervalli (ad ogni test N raddoppia, quindi h si dimezza)
    int N_values[] = {6, 12, 24, 48, 96, 192, 384, 768, 1536};
    const int ntests = sizeof(N_values) / sizeof(N_values[0]);

    conv << fixed << setprecision(6);
    conv << "# h" << "\t\t\t\t\t\tE_max" << "\t\t\t\t\tratio (R)" << "\t\t\t\torder (p)" << "\n";

    double old_error = 0.0;
    double old_h = 0.0;

    for (int i = 0; i < ntests; i++)
    {
        int N = N_values[i];

        // passo di integrazione
        double h = (xi_max_conv - eps_conv) / N;

        // errore massimo della soluzione numerica
        double error = ErrorRK4(N);

        if (i == 0) 
        {
            // per il primo valore non si ha ancora un errore precedente con cui confrontarlo
            conv << scientific << setprecision(12);
            conv << h << "\t\t" << error << "\t\t" << "NaN" << "\t\t\t\t" << "NaN" << "\n";
        }

        else 
        {
            // rapporto tra errori consecutivi (atteso 16)
            double ratio = old_error / error;

            // ordine di convergenza osservato
            double order = log(old_error / error) / log(old_h / h);

            conv << h << "\t\t" << error << "\t\t" << ratio << "\t\t" << order << "\n";
        }

        old_error = error;
        old_h = h;
    }

    conv.close();

    return 0;
}