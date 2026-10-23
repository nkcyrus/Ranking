#include <iostream>
#include <cstdio>
#include <chrono>
#include <cmath>
#include <omp.h>

#include <iomanip>

#include "ranking.cpp"
extern "C" {
    #include "mmio.c"
}

using namespace std;

struct Matrice{
    int n;
    long long n_links;
    double* values;
    long long* col_ptr;
    int* row_idx;
    int* degree_Sort_Tot;
    double* nabla;
    double* delta;
};

// Lit un entier et avance le pointeur
inline int fast_atoi(char*& p) {
    int val = 0;
    while (*p && (*p < '0' || *p > '9')) p++; // Skip non-digits
    while (*p >= '0' && *p <= '9') {
        val = val * 10 + (*p - '0');
        p++;
    }
    return val;
}

// Lit un double et avance le pointeur (si besoin)
inline double fast_atof(char*& p) {
    while (*p && *p == ' ') p++;
    char* end;
    double val = strtod(p, &end);
    p = end;
    return val;
}

Matrice* read_mtx2(const char * filename) {
    FILE* file = fopen(filename, "r");
    if (!file) return nullptr;

    char sys_buffer[1 << 20]; // 1 Mo est parfait pour ton cache L3
    setvbuf(file, sys_buffer, _IOFBF, sizeof(sys_buffer));
    
    MM_typecode matcode;
    int M, N, nz_int;
    mm_read_banner(file, &matcode);
    mm_read_mtx_crd_size(file, &M, &N, &nz_int);
    long long nz = (long long)nz_int;

    Matrice* matrice = new Matrice;
    matrice->n = M;
    matrice->n_links = nz;
    matrice->degree_Sort_Tot = new int[N]();
    int* nb_liens_par_colonne = new int[N]();

    matrice->nabla = new double[N]();
    matrice->delta = new double[N]();


    // Tableaux temporaires pour stocker les coordonnées en RAM
    int* temp_rows = new int[nz];
    int* temp_cols = new int[nz];
    double* temp_vals = nullptr;
    // double* sum_weight_out = nullptr;
    if (!mm_is_pattern(matcode)){
        temp_vals = new double[nz];
        // sum_weight_out = new double[N]();
    }
    char line[MM_MAX_LINE_LENGTH];
    
    // --- SEUL ET UNIQUE PASSAGE SUR LE DISQUE ---
    for (long long i = 0; i < nz; i++) {
        if (!fgets(line, sizeof(line), file)) break;
        char* p = line;
        int r = fast_atoi(p);
        int c = fast_atoi(p);
        // double val = 0;

        temp_rows[i] = r - 1;
        temp_cols[i] = c - 1;
        
        matrice->degree_Sort_Tot[r - 1]++;
        nb_liens_par_colonne[c - 1]++;
        // if (!mm_is_pattern(matcode)){
        //     val = fast_atof(p);
        //     sum_weight_out[r - 1] += val;
        // }
        
        
        if (temp_vals) temp_vals[i] = fast_atof(p);
    }
    fclose(file); // On libère le fichier immédiatement

    // Construction des pointeurs de colonnes (CSR)
    matrice->col_ptr = new long long[N + 1]();
    for (int j = 0; j < N; j++) {
        matrice->col_ptr[j + 1] = matrice->col_ptr[j] + nb_liens_par_colonne[j];
    }

    // Remplissage final de la structure CSR depuis la RAM
    matrice->row_idx = new int[nz];
    matrice->values = new double[nz]();

    long long* current_pos = new long long[N]();
    
    // Cette boucle est maintenant purement en RAM -> Très rapide
    for (long long i = 0; i < nz; i++) {
        int col = temp_cols[i];
        int row = temp_rows[i];
        long long pos = matrice->col_ptr[col] + current_pos[col];
        matrice->row_idx[pos] = row;
        // matrice->values[pos] = 1.0/matrice->degree_Sort_Tot[row];
        if (temp_vals) matrice->values[pos] = temp_vals[i];
        current_pos[col]++;
    }

    // Nettoyage des tableaux temporaires
    delete[] temp_rows;
    delete[] temp_cols;
    if (temp_vals) delete[] temp_vals;
    // delete[] sum_weight_out;

    delete[] current_pos;
    delete[] nb_liens_par_colonne;

    return matrice;
}

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

double calcul_G(Matrice* matrice, double* X_old, double* Y_old, double* X_new, double* Y_new, double alpha, double prem_term, double sec_term, double* normes) {
    int N = matrice->n;

    double xf = multiplier_vecteurs(X_old, matrice->degree_Sort_Tot, N);
    double yf = multiplier_vecteurs(Y_old, matrice->degree_Sort_Tot, N);

    double partie_constante_x = (prem_term*xf) + sec_term;
    double partie_constante_y = (prem_term*yf) + sec_term;


    double norme_totale = 0.0;
    double norme_l1_x = 0.0, norme_l1_y = 0.0;
    #pragma omp parallel for reduction(+:norme_totale, norme_l1_x, norme_l1_y) schedule(dynamic, 1024)
    for (int j = 0; j<N; j++){
        double somme_x =0.0, somme_y = 0.0;
        
        for(long long i = matrice->col_ptr[j]; i < matrice->col_ptr[j+1]; i++){
           int row = matrice->row_idx[i];
            
            somme_x += X_old[row] * matrice->values[i];
            somme_y += Y_old[row] * matrice->values[i];
        }
        double nouveau_score_x = (alpha*somme_x) + partie_constante_x;
        double nouveau_score_y = (alpha*somme_y) + partie_constante_y;
        
        double max_x = max(matrice->nabla[j], nouveau_score_x);
        double min_y = min(matrice->delta[j], nouveau_score_y);
        if (j < 5) {
    cout << fixed << setprecision(10);
    cout << "j=" << j << " nabla=" << matrice->nabla[j] 
         << " somme_x=" << somme_x
         << " XG=" << nouveau_score_x << endl;
}

        X_new[j] = max_x;
        Y_new[j] = min_y;
        norme_l1_x += max_x;
        norme_l1_y += min_y;
        

      

    }

    double invNormeX = 1.0 / norme_l1_x;
    double invNormeY = 1.0 / norme_l1_y;

    // #pragma omp parallel for reduction(+:norme_totale) schedule(dynamic, 1024)
    for (int i = 0; i < N; i++) {
        X_new[i] *= invNormeX; // Normalisation de X_new
        Y_new[i] *= invNormeY; // Normalisation de Y_new
        normes[i] = Y_new[i] - X_new[i];
        norme_totale += fabs(normes[i]); // Calcul de la norme L1 de la différence
    }
    return norme_totale;

}

void initialiser_nabla_delta(Matrice* P, double alpha, double& sum_nabla, double& sum_delta, double sec_term) {

    cout << "Initialisation des limites" << endl;
    int N = P->n;
    double invN = 1.0 / (double)N;

    sum_nabla = 0.0;
    sum_delta = 0.0;

    //Initialisation des limites
    // #pragma omp parallel for
    //  for (int j = 0; j < N; j++) {
    //     P->nabla[j] = numeric_limits<double>::infinity();
    //     P->delta[j] = -numeric_limits<double>::infinity();
    // }
    

    // Parcours CSC colonne par colonne
    #pragma omp parallel for reduction(+:sum_nabla, sum_delta) schedule(dynamic, 1024)
    for (int j = 0; j < N; j++) {
        double local_min = numeric_limits<double>::infinity();
        double local_max = -numeric_limits<double>::infinity();
        bool is_empty = true;
        for (long long i = P->col_ptr[j]; i < P->col_ptr[j+1]; i++) {
            is_empty = false;
            int row = P->row_idx[i];
            double deg = P->degree_Sort_Tot[row];
            P->values[i] = (deg == 0) ? 0 : (1.0 / deg);
            double g_ij = 0.0;
            if (deg == 0) g_ij = invN;
            else g_ij = alpha * P->values[i] + sec_term;
        
            local_min = min(local_min, g_ij);
            local_max = max(local_max, g_ij);
        }

        if (is_empty) {
            P->nabla[j] = invN;
            P->delta[j] = invN;
        } else {
            P->nabla[j] = local_min;
            P->delta[j] = local_max;
        }
      

        sum_nabla += P->nabla[j];
        sum_delta += P->delta[j];
    }

    cout << "Fin de l'initialisation des limites" << endl;
    
}

void print_vecteur(double* v, int n){
    cout << "pi = (";
    for(int i=0; i<n; i++){
        cout << v[i];
        if(i!=n-1) cout<< ", ";
    }
    cout << ")" << endl;
}

void calcul_pagerank_G(Matrice* P, double epsilon, double alpha, int max_iter, int max_temps){
    int N= P->n;
    double sum_nabla = 0.0, sum_delta = 0.0;

    cout << "Début de des calculs" << endl;
    
    double prem_term = alpha/N;
    double sec_term = (1.0-alpha)/N;

    initialiser_nabla_delta(P, alpha, sum_nabla, sum_delta, sec_term);

    double* X_old = new double[N];
    double* Y_old = new double[N];
    double* X_new = new double[N];
    double* Y_new = new double[N];
    double* normes = new double[N];

    double invNormeNabla = 1.0 / sum_nabla;
    double invNormeDelta = 1.0 / sum_delta;
    cout << "Somme totale de Nabla : " << sum_nabla << endl;
    #pragma omp parallel for
    for (int i = 0; i < N; i++) {
        X_old[i] = P->nabla[i] * invNormeNabla; // Initialisation de X_old à nabla
        Y_old[i] = P->delta[i] * invNormeDelta; // Initialisation de Y_old à delta
    }


    string raison_stop = "Maximum d'itérations atteint";
    int cpt = 0;
    double norme = 0.0;
    double ancienne_norme = 0.0;


    auto start = chrono::steady_clock::now();

    while (cpt <max_iter)
    {
        norme = calcul_G(P, X_old, Y_old, X_new, Y_new, alpha, prem_term, sec_term, normes);

        if (norme < epsilon) {
            raison_stop = "Convergence atteinte (Epsilon)";
            cpt++;
            break;
        }

        if (fabs(ancienne_norme - norme) < (epsilon * 0.01)) {
            raison_stop = "Stagnation détectée (Plateau)";
            cpt++;
            break;
        }

        auto current = chrono::steady_clock::now();
        chrono::duration<double> elapsed = current - start;
        if (elapsed.count() > max_temps) {
            raison_stop = "Temps limite dépassé (Timeout)";
            cpt++;
            break;
        }
        swap(X_old, X_new);
        swap(Y_old, Y_new);
        cpt++;
        ancienne_norme = norme;

    }

    cout << "Fin de des calculs" << endl;
    auto end = chrono::steady_clock::now();
    chrono::duration<double, milli> total_ms = end - start;

    cout << "Nombre d'iterations : " << cpt << endl;
    cout << "Statut : " << raison_stop << endl;
    cout << "Temps : " << total_ms.count() << " ms" << endl;
    cout << "Norme : " << scientific << setprecision(10) << norme << endl;
    // print_vecteur(pi_pair, N);
    // print_vecteur(pi_impair, N);
    // print_vecteur(f, N);

    double somme_X = 0.0;
    double somme_Y = 0.0;
    double* pi;
    pi = new double[N];
    #pragma omp parallel for reduction(+:somme_X, somme_Y)
    for(int i=0; i<N; i++) {
        somme_X += X_old[i];
        somme_Y += Y_old[i];
        pi[i] = Y_old[i] - 0.5 * normes[i]; // On peut choisir d'afficher X_old ou Y_old, ici j'ai choisi X_old
    }

    cout << "Somme X : " << somme_X << endl;
    cout << "Somme Y : " << somme_Y << endl;

    int inversion_count = 0;
    for (int i = 0; i < N; ++i) {
        if (P->nabla[i] > P->delta[i]) {
            inversion_count++;
            // Log les 10 premiers pour voir le pattern
            if (inversion_count <= 10) {
                std::cout << "Erreur sommet #" << i + 1
                        << " | Nabla: " << P->nabla[i] 
                        << " | Delta: " << P->delta[i] << std::endl;
            }
        }
    }
    std::cout << "--- BILAN INITIALISATION ---" << std::endl;
    std::cout << "Nombre total de sommets avec Nabla > Delta : " << inversion_count << std::endl;

    // Appel de ta nouvelle méthode
    extraireEtAfficherTop10(X_old, N);
    extraireEtAfficherTop10(Y_old, N);
    extraireEtAfficherTop10(pi, N);

   
   

    delete[] X_old;
    delete[] Y_old;
    delete[] X_new;
    delete[] Y_new;
    delete[] normes;
    delete[] pi;
}

// void print_vecteur(double* v, int n){
//     cout << "pi = (";
//     for(int i=0; i<n; i++){
//         cout << v[i];
//         if(i!=n-1) cout<< ", ";
//     }
//     cout << ")" << endl;
// }

void delete_matrix(Matrice* matrice) {
    if (matrice == nullptr) return;

    // 1. Libérer les tableaux de données
    if (matrice->col_ptr) {
        delete[] matrice->col_ptr;
        matrice->col_ptr = nullptr;
    }

    if (matrice->row_idx) {
        delete[] matrice->row_idx;
        matrice->row_idx = nullptr;
    }

    if (matrice->nabla) {
        delete[] matrice->nabla;
        matrice->nabla = nullptr;
    }
    if (matrice->delta) {
        delete[] matrice->delta;
        matrice->delta = nullptr;
    }

    if (matrice->values) {
        delete[] matrice->values;
        matrice->values = nullptr;
    }

    if (matrice->degree_Sort_Tot) {
        delete[] matrice->degree_Sort_Tot;
        matrice->degree_Sort_Tot = nullptr;
    }

    // 2. Libérer la structure elle-même
    delete matrice;
    
}

int main(int argc, char** argv){

    double epsilon = 10e-8;
    double alpha = 0.85;
    int max_iter = 1000, max_temps = 300;
    if (argc < 2) {
        std::cerr << "Usage : " << argv[0] << " <fichier_matrice.mtx>" << std::endl;
        std::cerr << "Exemple : ./pagerank web-Stanford.mtx" << std::endl;
        return 1; // On quitte avec un code d'erreur
    }

    const char* filename = argv[1];

    if (argc >= 3) epsilon  = atoi(argv[2]);
    // if (argc >= 4) max_temps = atoi(argv[3]);

    std::cout << "--- Initialisation du PageRank ---" << std::endl;
    std::cout << "Fichier cible : " << filename << std::endl;
    auto start_load = chrono::steady_clock::now();
    
    Matrice* myMatrix = read_mtx2(filename);
    
    auto end_load = chrono::steady_clock::now();

    if (myMatrix == nullptr) {
        return 1;
    }

    chrono::duration<double> elapsed_load = end_load - start_load;

    cout << "Matrice chargee avec succes." << endl;
    cout << "Nombre de pages (N) : " << myMatrix->n << endl;
    cout << "Nombre de liens (nz) : " << myMatrix->n_links << endl;
    cout << "1/N : " << 1.0 / (double) myMatrix->n << endl;
    cout << fixed << setprecision(2);
    cout << "==> Temps de chargement : " << elapsed_load.count() << " secondes" << endl;
    cout << "-----------------------------------" << endl;

    // if (myMatrix->n < 20) {
    //     print_matrice(myMatrix);
    // } else {
    //     std::cout << "Matrice chargee (" << myMatrix->n 
    //               << " pages). Trop grande pour un affichage complet." << std::endl;
    // }

    calcul_pagerank_G(myMatrix, epsilon, alpha, max_iter, max_temps);
    delete_matrix(myMatrix);

    return 0;
}