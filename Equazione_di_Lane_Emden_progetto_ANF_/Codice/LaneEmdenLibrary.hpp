#ifndef LANE_EMDEN_LIBRARY_HPP
#define LANE_EMDEN_LIBRARY_HPP

// funzione per il calcolo del RHS
// restituisce il primo zero di theta(xi), -1 altrimenti
void LaneEmdenRHS(double xi, double *Y, double *R);

// funzione per l'integrazione ODE (RK4)
double SolveLaneEmden(double n, double dxi, double xi_max, double eps);

// funzione per il calcolo della massa della nana bianca
double PolytropeMass(double n, double K, double rho_c, double xi_s, double thetaPrime_s);

// funzione per il calcolo del raggio della nana bianca
double PolytropeRadius(double n, double K, double rho_c, double xi_s);

// funzione per il test convergenza RK4 per n = 1
// massimo errore assoluto per n = 1 con N passi
// confronta con sin(xi)/xi
double ErrorRK4(int N);

// indice politropico per RHS
extern double g_n;

extern int g_found_surface; // 1 se il primo zero e' stato trovato, 0 altrimenti
extern double g_thetaPrime_s; // theta'(xi_s), valida solo se g_found_surface = 1
 
// massimi errori assoluti per i casi con soluzione analitica.
extern double g_max_err_n0;
extern double g_max_err_n1;
extern double g_max_err_n5;

#endif