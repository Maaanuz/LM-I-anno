#include <cmath>

#include "ODEsolvers.hpp"
#include "LaneEmdenLibrary.hpp"

// costanti fisiche
const double G = 6.67430e-11;

// variabili globali per LaneEmdenRHS e SolveLaneEmden
double g_n = 1.0;
int g_found_surface = 0;
double g_thetaPrime_s = 0.0; // theta'(xi_s)
double g_max_err_n0 = 0.0; // max errore per n=0
double g_max_err_n1 = 0.0; // max errore per n=1
double g_max_err_n5 = 0.0; // max errore per n=5

// funzione per il calcolo del RHS
// sistema del primo ordine con Y = {theta, theta'} e R = {theta', theta''}
void LaneEmdenRHS(double xi, double *Y, double *R)
{
    const double theta  = Y[0];
    const double dtheta = Y[1];
    R[0] = dtheta;
    double thetan = 0.0;

    // per theta <= 0 si annulla theta^n per evitare potenze con n non intero
    if (theta > 0.0) 
        thetan = pow(theta, g_n);

    // l'integrazione parte da eps > 0
    if (xi > 0.0) 
        R[1] = -(2.0 / xi) * dtheta - thetan;
    else 
        R[1] = 0.0;
}

// funzione per l'integrazione ODE
double SolveLaneEmden(double n, double dxi, double xi_max, double eps)
{
    // impostazione indice politropico 
    // azzeramento dell'integrazione precedente
    g_n = n;
    g_found_surface = 0;
    g_thetaPrime_s = 0.0;
    g_max_err_n0 = 0.0;
    g_max_err_n1 = 0.0;
    g_max_err_n5 = 0.0;

    const double eps2 = eps * eps;
    const double eps4 = eps2 * eps2;

    double xi = eps;
    double Y[2];

    // sviluppo in serie intorno a xi = 0
    Y[0] = 1.0 - eps2 / 6.0 + (n * eps4) / 120.0;
    Y[1] = -eps / 3.0 + (n * eps2 * eps) / 30.0;

    double xi_prev = xi;
    double th_prev = Y[0];
    double dth_prev = Y[1];

    const int neq = 2;
    const int N = (int)((xi_max - eps) / dxi); // numero di passi senza superare xi_max

    for (int i = 0; i < N; i++) {

        RK4(xi, Y, LaneEmdenRHS, dxi, neq);

        xi = eps + (i + 1) * dxi;

        // -------------------------------
        // aggiornamenti massimo errore su theta

        // n = 0, theta_0 = 1 - xi^2/6
        if (fabs(n - 0.0) < 1e-15) {

            double ana0 = 1.0 - (xi * xi) / 6.0;
            double err0 = fabs(Y[0] - ana0);

            if (err0 > g_max_err_n0)
                g_max_err_n0 = err0;
        }

        // n = 1, theta_1 = sin(xi)/xi
        if (fabs(n - 1.0) < 1e-15) {

            double ana1 = sin(xi) / xi;
            double err1 = fabs(Y[0] - ana1);

            if (err1 > g_max_err_n1) 
                g_max_err_n1 = err1;
        }

        // n = 5, theta_5 = 1/sqrt(1 + xi^2/3)
        if (fabs(n - 5.0) < 1e-15) {

            double ana5 = 1.0 / sqrt(1.0 + (xi * xi) / 3.0);
            double err5 = fabs(Y[0] - ana5);

            if (err5 > g_max_err_n5)
                g_max_err_n5 = err5;
        }
        // -------------------------------

        // ricerca del primo zero di theta
        if (th_prev > 0.0 && Y[0] <= 0.0) {
            
            // interpolazione lineare per xi_s
            double xi_s = xi_prev - th_prev * (xi - xi_prev) / (Y[0] - th_prev);

            // interpolazione lineare per theta'(xi_s)
            double thetaPrime_s = dth_prev + (Y[1] - dth_prev) * (xi_s - xi_prev) / (xi - xi_prev);

            g_thetaPrime_s = thetaPrime_s;
            g_found_surface = 1;

            return xi_s;
        }

        // salvataggio punto corrente per il prossimo passo
        xi_prev  = xi;
        th_prev  = Y[0];
        dth_prev = Y[1];
    }

    // superficie non trovata entro xi_max
    return -1.0;
}

// funzione per il calcolo della massa 
double PolytropeMass(double n, double K, double rho_c, double xi_s, double thetaPrime_s)
{
    const double factor = (n + 1.0) * K / (4.0 * M_PI * G);
    const double a0 = sqrt(factor);
    const double expn = (1.0 - n) / (2.0 * n);
    const double a = a0 * pow(rho_c, expn);
    const double term = - (xi_s * xi_s) * thetaPrime_s; // adimensionale

    return 4.0 * M_PI * a * a * a * rho_c * term; // [kg]
}

// funzione per il calcolo del raggio
double PolytropeRadius(double n, double K, double rho_c, double xi_s)
{
    const double factor = (n + 1.0) * K / (4.0 * M_PI * G);
    const double a0 = sqrt(factor);
    const double expn = (1.0 - n)/(2.0*n);
    const double a = a0 * pow(rho_c, expn);

    return a * xi_s; // [m]
}

// funzione per il test convergenza RK4 per n = 1
double ErrorRK4(int N)
{
    const double eps = 1e-6;
    const double xi_max = 3.0;
    const double n = 1.0;

    g_n = n;

    double eps2 = eps * eps;
    double eps4 = eps2 * eps2;
    double xi = eps;

    double Y[2];

    // condizioni iniziali da sviluppo intorno a 0
    Y[0] = 1.0 - eps2 / 6.0 + n * eps4 / 120.0;
    Y[1] = -eps / 3.0 + n * eps2 * eps / 30.0;

    // N passi sull'intervallo [eps, xi_max]
    double h = (xi_max - eps) / N;

    const int neq = 2;

    double max_error = 0.0;

    for (int i = 0; i < N; i++)
    {
        RK4(xi, Y, LaneEmdenRHS, h, neq);

        xi += h;

        // soluzione analitica per n = 1
        double theta_exact = sin(xi) / xi;

        // errore assoluto nel punto corrente
        double error = fabs(Y[0] - theta_exact);

        // si conserva l'errore massimo
        if (error > max_error)
            max_error = error;
    }

    return max_error;
}