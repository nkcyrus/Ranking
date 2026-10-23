#include <iostream>
#include <cstdio>
#include <chrono>
#include <cmath>

#include <iomanip>



using namespace std;

// Matrice* read_mtx(const char * filename){
//     FILE* file = fopen(filename, "r");
//     if (!file) {
//         cerr << "Nous n'avons pas pu ouvrir le fichier" << filename << endl;
//         return nullptr;
//     }

//     char sys_buffer[1 << 20]; // 1 Mo
//     setvbuf(file, sys_buffer, _IOFBF, sizeof(sys_buffer));
    
//     MM_typecode matcode;
//     int M, N, nz_int;
    

//     mm_read_banner(file, &matcode);
//     mm_read_mtx_crd_size(file, &M, &N, &nz_int);
//     long long nz= (long long) nz_int;

//     if (M != N) {
//         cerr << "La matrice n'est pas carrée" << endl;
//         fclose(file);
//         return nullptr;
//     }


//     Matrice* matrice = new Matrice;
//     matrice->n = M;
//     matrice->n_links = nz;

//     int* nb_liens_par_colonne = new int[N]();
//     matrice->col_ptr = new long long[N + 1]();
//     matrice->row_idx = new int[nz];
//     matrice->degree_Sort_Tot = new int[N]();

//     char line[MM_MAX_LINE_LENGTH];
    

    
//     for (long long i = 0; i < nz; i++) {

//         fgets(line, sizeof(line), file);

//         char* p = line;
//         int row = fast_atoi(p);
//         int col = fast_atoi(p);

//         matrice->degree_Sort_Tot[row - 1]++;
//         nb_liens_par_colonne[col-1]++;
//     }

//     matrice->col_ptr[0] = 0;

//     for (int j = 0; j < N; j++) {
//         matrice->col_ptr[j + 1] = matrice->col_ptr[j] + nb_liens_par_colonne[j];
//     }

//     long long* current_position = new long long[N];
//     #pragma omp parallel for
//     for (int i = 0; i < N; i++) current_position[i] = 0;

//     rewind(file);
//     mm_read_banner(file, &matcode);
//     mm_read_mtx_crd_size(file, &M, &N, &nz_int);

//     if (!mm_is_pattern(matcode)) {
//         matrice->values = new double[nz];
//     } else {
//         matrice->values = nullptr;
//     }

//     for (long long i = 0; i < nz; i++) {

//         fgets(line, sizeof(line), file);
//         char* p = line;
//         int row = fast_atoi(p);
//         int col = fast_atoi(p);


//         int dest = col - 1;
//         long long pos = matrice->col_ptr[dest] + current_position[dest];
//         matrice->row_idx[pos] = row - 1;
//         if (!mm_is_pattern(matcode)) {
//             matrice->values[pos] = fast_atof(p);
//         }
        
//         current_position[dest]++;
//     }


//     delete[] current_position;
//     delete[] nb_liens_par_colonne;
//     fclose(file);
//     return matrice;
// }

// Matrice* read_mtx(const char * filename) {
//     FILE* file = fopen(filename, "r");
//     if (!file) return nullptr;

//     char sys_buffer[1 << 20]; // 1 Mo est parfait pour ton cache L3
//     setvbuf(file, sys_buffer, _IOFBF, sizeof(sys_buffer));
    
//     MM_typecode matcode;
//     int M, N, nz_int;
//     mm_read_banner(file, &matcode);
//     mm_read_mtx_crd_size(file, &M, &N, &nz_int);
//     long long nz = (long long)nz_int;

//     Matrice* matrice = new Matrice;
//     matrice->n = M;
//     matrice->n_links = nz;
//     matrice->degree_Sort_Tot = new int[N]();
//     int* nb_liens_par_colonne = new int[N]();

//     // Tableaux temporaires pour stocker les coordonnées en RAM
//     int* temp_rows = new int[nz];
//     int* temp_cols = new int[nz];
//     double* temp_vals = new double[nz];

//     char line[MM_MAX_LINE_LENGTH];
    
//     // --- SEUL ET UNIQUE PASSAGE SUR LE DISQUE ---
//     for (long long i = 0; i < nz; i++) {
//         if (!fgets(line, sizeof(line), file)) break;
//         char* p = line;
//         int r = fast_atoi(p);
//         int c = fast_atoi(p);

//         temp_rows[i] = r - 1;
//         temp_cols[i] = c - 1;
        
//         matrice->degree_Sort_Tot[r - 1]++;
//         nb_liens_par_colonne[c - 1]++;
        
//         if (temp_vals) temp_vals[i] = fast_atof(p);
//         else temp_vals[i] = 1;
//     }
//     fclose(file); // On libère le fichier immédiatement

//     // Construction des pointeurs de colonnes (CSR)
//     matrice->col_ptr = new long long[N + 1]();
//     for (int j = 0; j < N; j++) {
//         matrice->col_ptr[j + 1] = matrice->col_ptr[j] + nb_liens_par_colonne[j];
//     }

//     // Remplissage final de la structure CSR depuis la RAM
//     matrice->row_idx = new int[nz];
//     if (temp_vals) matrice->values = new double[nz];
//     else matrice->values = nullptr;

//     long long* current_pos = new long long[N]();
    
//     // Cette boucle est maintenant purement en RAM -> Très rapide
//     for (long long i = 0; i < nz; i++) {
//         int col = temp_cols[i];
//         long long pos = matrice->col_ptr[col] + current_pos[col];
//         matrice->row_idx[pos] = temp_rows[i];
//         if (temp_vals) matrice->values[pos] = temp_vals[i];
//         current_pos[col]++;
//     }

//     // Nettoyage des tableaux temporaires
//     delete[] temp_rows;
//     delete[] temp_cols;
//     if (temp_vals) delete[] temp_vals;
//     delete[] current_pos;
//     delete[] nb_liens_par_colonne;

//     return matrice;
// }


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




// int main(int argc, char** argv){

//     double epsilon = 1e-8;
//     double alpha = 0.85;
//     int max_iter = 200, max_temps = 30;
//     if (argc < 2) {
//         std::cerr << "Usage : " << argv[0] << " <fichier_matrice.mtx>" << std::endl;
//         std::cerr << "Exemple : ./pagerank web-Stanford.mtx" << std::endl;
//         return 1; // On quitte avec un code d'erreur
//     }

//     const char* filename = argv[1];

//     if (argc >= 3) epsilon = atof(argv[2]);
//     // if (argc >= 3) max_iter  = atoi(argv[2]);
//     // if (argc >= 4) max_temps = atoi(argv[3]);

//     std::cout << "--- Initialisation du PageRank ---" << std::endl;
//     std::cout << "Fichier cible : " << filename << std::endl;
//     auto start_load = chrono::steady_clock::now();
    
//     Matrice* myMatrix = read_mtx2(filename);
    
//     auto end_load = chrono::steady_clock::now();

//     if (myMatrix == nullptr) {
//         return 1;
//     }

//     chrono::duration<double> elapsed_load = end_load - start_load;

//     cout << "Matrice chargee avec succes." << endl;
//     cout << "Nombre de pages (N) : " << myMatrix->n << endl;
//     cout << "Nombre de liens (nz) : " << myMatrix->n_links << endl;
//     cout << fixed << setprecision(2);
//     cout << "==> Temps de chargement : " << elapsed_load.count() << " secondes" << endl;
//     cout << "-----------------------------------" << endl;

//     // if (myMatrix->n < 20) {
//     //     print_matrice(myMatrix);
//     // } else {
//     //     std::cout << "Matrice chargee (" << myMatrix->n 
//     //               << " pages). Trop grande pour un affichage complet." << std::endl;
//     // }

//     // calcul_pagerank_classique(myMatrix, epsilon, alpha);
//     calcul_pagerank_G(myMatrix, epsilon, alpha, max_iter, max_temps);
//     delete_matrix(myMatrix);

//     return 0;
// }