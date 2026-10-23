#include <iostream>
#include <cstdio>
#include <chrono>
#include <cmath>

#include <iomanip>



using namespace std;


double multiplier_vecteurs(double* x, int* f, int n){
    double result= 0.0;

    #pragma omp parallel for reduction(+:result)
    for (int i = 0; i < n; i++)
    {
        if (f[i] == 0) {
            result += x[i];
        }
    }
    return result;
}

double calcul_G(Matrice* matrice, double* vecteur, double* result, double prem_term, double sec_term, double alpha){
    int N = matrice->n;

    double xf = multiplier_vecteurs(vecteur, matrice->degree_Sort_Tot, N);
    double partie_constante = (prem_term*xf) + sec_term;

    double norme_totale = 0.0;

    #pragma omp parallel for reduction(+:norme_totale) schedule(dynamic, 1024)
    for (int j = 0; j<N; j++){
        double somme =0.0;
        
        for(long long i = matrice->col_ptr[j]; i < matrice->col_ptr[j+1]; i++){
           int row = matrice->row_idx[i];
            
            somme += vecteur[row] * matrice->values[i];
        }
        double nouveau_score = (alpha*somme) + partie_constante;
        result[j] = nouveau_score;
        double diff = nouveau_score - vecteur[j];
        norme_totale += (diff < 0) ? -diff : diff;
    }

    return norme_totale;

}


Result* calcul_pagerank(Matrice* P, double epsilon, double alpha, int max_iter, int max_temps){


    int N= P->n;

    double* pi_pair = new double[N];
    double* pi_impair = new double[N];

    double prem_term = alpha/N;
    double sec_term = (1.0-alpha)/N;

    #pragma omp parallel for
    for (int i = 0; i < N; i++) pi_pair[i] = 1.0/N;
    double norme = 1.0;
    double ancienne_norme = -1.0;
    int cpt = 0;

    string raison_stop = "Maximum d'itérations atteint";

    auto start = chrono::steady_clock::now();

    while (cpt <max_iter)
    {
        norme = calcul_G(P, pi_pair, pi_impair, prem_term, sec_term, alpha);
        cpt++;
        swap(pi_pair, pi_impair);

        if (norme < epsilon) {
            raison_stop = "Convergence atteinte (Epsilon)";
            break;
        }

        if (fabs(ancienne_norme - norme) < (epsilon * 0.01)) {
            raison_stop = "Stagnation détectée (Plateau)";
            break;
        }

        auto current = chrono::steady_clock::now();
        chrono::duration<double> elapsed = current - start;
        if (elapsed.count() > max_temps) {
            raison_stop = "Temps limite dépassé (Timeout)";
            break;
        }
        ancienne_norme = norme;

    }
    auto end = chrono::steady_clock::now();
    chrono::duration<double, milli> total_ms = end - start;

   
    long long id = start.time_since_epoch().count();

    extraireTop10(pi_pair, N, id);

    delete[] pi_impair;
    delete[] pi_pair;

    Result* result = new Result;
    result->id = id;
    result->iterations = cpt;
    result->time_ms = total_ms.count();
    result->raison = raison_stop;

    return result;
    
}