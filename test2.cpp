#include <iostream>
#include <cstdio>
#include <chrono>
#include <cmath>
#include <omp.h>
#include <iomanip>
#include <limits>
#include <algorithm>


using namespace std;


void multiplier_vecteurs(double* xy, int* deg, int n, double& result_x, double& result_y) {
    #pragma omp parallel for reduction(+:result_x, result_y) schedule(dynamic, 1024)
    for (int i = 0; i < n; ++i) { // On boucle directement sur le nœud
    if (deg[i] == 0) {
        result_x += xy[2*i];
        result_y += xy[2*i + 1];
    }
}    
}


// ============================================================
//  Initialisation nabla, delta et values (P = 1/deg)
//
//  G[i][j] = α·(1/deg(i)) + (α/N)·f[i] + (1-α)/N    si lien i→j
//           = (α/N)·f[i] + (1-α)/N                    si pas de lien (mais i non-dangling)
//           = 1/N                                       si i dangling
//
//  Pour le min/max sur la colonne j :
//    - Les lignes i ABSENTES de la colonne j contribuent :
//        α·0 + (α/N)·f[i] + (1-α)/N  = (α/N)·f[i] + (1-α)/N
//      Le min de ces termes est sec_term = (1-α)/N  (pour f[i]=0)
//      Le max est invN = 1/N            (pour f[i]=1, dangling)
//    - Les lignes i PRÉSENTES contribuent :
//        α/deg(i) + sec_term   (car f[i]=0 si deg(i)>0)
//
//  Donc :
//    nabla[j] = min(sec_term,  min sur présents de α/deg(i)+sec_term)
//             = sec_term   (car α/deg(i) ≥ 0)
//    delta[j] = max(invN,     max sur présents de α/deg(i)+sec_term)
//             = max sur présents  (car α/deg(i)+sec_term ≥ invN pour deg≥1 et α=0.85)
//
//  Si la colonne est vide (aucun lien entrant) : toutes les lignes sont absentes
//    → nabla[j] = sec_term, delta[j] = invN
// ============================================================
void initialiser_nabla_delta(Matrice* P, double alpha,
                              double& sum_nabla, double& sum_delta,
                              double sec_term, double a) {
    int N = P->n;
    double invN = 1.0 / (double)N;

    P->nd = new double[2*N];


    sum_nabla = 0.0;
    sum_delta = 0.0;

    #pragma omp parallel for reduction(+:sum_nabla, sum_delta) schedule(dynamic, 1024)
    for (int j = 0; j < N; j++) {
        
        double local_min = sec_term;  // nœuds sans lien vers j → sec_term
        double local_max = sec_term;  // idem

        for (long long k = P->col_ptr[j]; k < P->col_ptr[j+1]; k++) {
            int i = P->row_idx[k];
            double g_ij = alpha * P->values[k] + sec_term;
            local_min = min(local_min, g_ij);
            local_max = max(local_max, g_ij);
        }

        local_min = min(local_min, invN);
        local_max = max(local_max, invN);

        P->nd[j*2] = local_min;
        P->nd[2*j+1] = local_max;

        sum_nabla += local_min;
        sum_delta += local_max;
    }

}

// ============================================================
//  Une itération : X(k+1), Y(k+1)
//  Retourne ||X(k+1) - Y(k+1)||_1
// ============================================================
double calcul_G(Matrice* matrice, double* XY_old, double* XY_new,
                double alpha, double prem_term, double sec_term, double a) {
    int N = matrice->n;

    double xf = 0.0, yf = 0.0;

    // xf = Σ x[i]·f[i]  (contribution des dangling nodes)
    multiplier_vecteurs(XY_old, matrice->degree_Sort_Tot, N, xf, yf);

    // Partie constante de (v·G)[j] indépendante de j
    double partie_constante_x = prem_term * xf + sec_term;
    double partie_constante_y = prem_term * yf + sec_term;

    double norme_l1_x = 0.0, norme_l1_y = 0.0;

    #pragma omp parallel for reduction(+:norme_l1_x, norme_l1_y) schedule(dynamic, 1024)
    for (int j = 0; j < N; j++) {
        double somme_x = 0.0, somme_y = 0.0;
        for (long long i = matrice->col_ptr[j]; i < matrice->col_ptr[j+1]; i++) {
            int row = matrice->row_idx[i];
            somme_x += XY_old[2*row] * matrice->values[i];
            somme_y += XY_old[2*row+1] * matrice->values[i];
        }

        // (X·G)[j] = α·somme_x + partie_constante_x
        double xg_j = alpha * somme_x + partie_constante_x;
        double yg_j = alpha * somme_y + partie_constante_y;

        // Encadrement : max avec nabla, min avec delta
        int pt = 2*j;
        XY_new[pt] = max(matrice->nd[pt], xg_j);
        XY_new[pt+1] = min(matrice->nd[pt+1], yg_j);

        norme_l1_x += XY_new[pt];
        norme_l1_y += XY_new[pt+1];
    }

    // Normalisation L1
    double invNormeX = 1.0 / norme_l1_x;
    double invNormeY = 1.0 / norme_l1_y;

    double norme_totale = 0.0;
    #pragma omp parallel for reduction(+:norme_totale)
    for (int i = 0; i < N; i++) {
        int pt = 2 * i;
        XY_new[pt] *= invNormeX;
        XY_new[pt+1] *= invNormeY;
        norme_totale += fabs(XY_new[pt+1] - XY_new[pt]);
    }
    return norme_totale;
}

// ============================================================
//  Boucle principale PageRank
// ============================================================
Result* calcul_pagerank2(Matrice* P, double epsilon, double alpha,
                       int max_iter, int max_temps) {
                        
    
    int N = P->n;
    double sum_nabla = 0.0, sum_delta = 0.0;

    double prem_term = alpha / N;       // α/N
    double sec_term  = (1.0 - alpha) / N; // (1-α)/N

    initialiser_nabla_delta(P, alpha, sum_nabla, sum_delta, sec_term, 1.0);

    double* XY_old  = new double[2*N];
    double* XY_new  = new double[2*N];
    // double* normes = new double[N]();

    double invNormeNabla = 1.0 / sum_nabla;
    double invNormeDelta = 1.0 / sum_delta;


    #pragma omp parallel for
    for (int i = 0; i < 2*N; i+=2) {
        XY_old[i] = P->nd[i] * invNormeNabla;
        XY_old[i+1] = P->nd[i+1] * invNormeDelta;
    }

    string raison_stop = "Maximum d'itérations atteint";
    int cpt = 0;
    double norme = 0.0;
    double ancienne_norme = -1.0; // -1 pour éviter faux plateau dès k=0

    auto start = chrono::steady_clock::now();

    while (cpt < max_iter) {
        norme = calcul_G(P, XY_old, XY_new, alpha, prem_term, sec_term, 1.0);
        swap(XY_old, XY_new);
        cpt++;

        if (norme < epsilon) {
            raison_stop = "Convergence atteinte (Epsilon)";
            // FIX 3 : swap avant break pour que X_old/Y_old = dernière itération
            break;
        }

        if (ancienne_norme >= 0.0 && fabs(ancienne_norme - norme) < epsilon * 0.01) {
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


    // logExperiment("log.csv", matriceName, P->n, P->n_links, alpha, epsilon, cpt, total_ms.count());



    
    // double* pi = new double[N];
    // double somme_X = 0.0, somme_Y = 0.0;

    // #pragma omp parallel for reduction(+:somme_X, somme_Y)
    // for (int i = 0; i < N; i++) {
    //     somme_X += X_old[i];
    //     somme_Y += Y_old[i];
    //     pi[i] = 0.5 * (X_old[i] + Y_old[i]);
    // }

    // cout << "Somme X : " << somme_X << endl;
    // cout << "Somme Y : " << somme_Y << endl;

    // // Vérification de cohérence nabla/delta
    // int inversion_count = 0;
    // for (int i = 0; i < N; ++i) {
    //     if (P->nabla[i] > P->delta[i]) {
    //         inversion_count++;
    //         if (inversion_count <= 10)
    //             cout << "Erreur sommet #" << i+1
    //                  << " | Nabla: " << P->nabla[i]
    //                  << " | Delta: " << P->delta[i] << endl;
    //     }
    // }

    long long id = start.time_since_epoch().count();

    extraireTop10_parallel(XY_old, N, id);

    delete[] XY_old;
    delete[] XY_new;
    delete[] P->nd;

    
    



    Result* result = new Result;
    result->id = id;
    result->iterations = cpt;
    result->time_ms = total_ms.count();
    result->raison = raison_stop;

    return result;
    
}



//  main(int argc, char** argv) {
//     double epsilon  = 1e-8;   // FIX 5 : 10e-8 = 1e-7, pas 1e-8
//     double alpha    = 0.85;
//     int max_iter    = 1000;
//     int max_temps   = 300;

//     if (argc < 2) {
//         cerr << "Usage : " << argv[0] << " <fichier_matrice.mtx> [epsilon] [max_temps]" << endl;
//         cerr << "Exemple : ./pagerank web-Stanford.mtx 1e-8 300" << endl;
//         return 1;
//     }

//     matriceName = argv[1];

//     // FIX 6 : atof et non atoi pour parser epsilon
//     if (argc >= 3) epsilon   = atof(argv[2]);
//     if (argc >= 4) max_temps = atoi(argv[3]);

//     cout << "--- Initialisation du PageRank ---" << endl;
//     cout << "Fichier  : " << matriceName << endl;
//     cout << "Alpha    : " << alpha    << endl;
//     cout << "Epsilon  : " << scientific << epsilon << endl;
//     cout << "MaxIter  : " << max_iter  << endl;
//     cout << "MaxTemps : " << max_temps << " s" << endl;

//     auto start_load = chrono::steady_clock::now();
//     Matrice* myMatrix = read_mtx2(matriceName);
//     auto end_load = chrono::steady_clock::now();

//     if (!myMatrix) { cerr << "Erreur de lecture du fichier." << endl; return 1; }

//     chrono::duration<double> elapsed_load = end_load - start_load;
//     cout << "Matrice chargée avec succès." << endl;
//     cout << "Nombre de pages (N)  : " << myMatrix->n      << endl;
//     cout << "Nombre de liens (nz) : " << myMatrix->n_links << endl;
//     cout << fixed << setprecision(3);
//     cout << "==> Temps de chargement : " << elapsed_load.count() << " s" << endl;
//     cout << "-----------------------------------" << endl;

//     calcul_pagerank_G(myMatrix, epsilon, alpha, max_iter, max_temps);
//     delete_matrix(myMatrix);

//     return 0;
// }