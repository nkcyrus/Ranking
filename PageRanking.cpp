#include <iostream>
#include <cstdio>
#include <chrono>

#include <iomanip>


extern "C" {
    #include "mmio.h"
    #include "mmio.c"
}

using namespace std;

struct Element {
    int row_index;
    double value;
    struct Element* next;
};

struct Matrice {
    int n;
    long long n_links;
    Element** columns;
};

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

Matrice* read_mtx(const char * filename, int*& degree_Sort_Tot){
    FILE* file = fopen(filename, "r");
    if (!file) {
        cerr <<  "Une erreur s'est produite lors de la lecture du fichier" << endl;
        return nullptr;
    }

    char sys_buffer[1 << 20];
    setvbuf(file, sys_buffer, _IOFBF, sizeof(sys_buffer));
    
    MM_typecode matcode;
    int M, N, nz_int;

    mm_read_banner(file, &matcode);
    mm_read_mtx_crd_size(file, &M, &N, &nz_int);
    long long nz = (long long)nz_int;

     if (M!=N) {
        cerr << "Matrice pas carré" << endl;
        return nullptr;
    }
    
    Matrice* matrice = new Matrice;
    matrice->n = M;
    matrice->n_links = nz;
    matrice->columns = new Element*[M];

    degree_Sort_Tot = new int[N];
 
    #pragma omp parallel for
    for(int i=0; i<M; i++){
        matrice->columns[i] = nullptr;
        degree_Sort_Tot[i] = 0 ;
    }
    
    char line[MM_MAX_LINE_LENGTH];
    for (long long cpt = 0; cpt < nz; cpt++){
        // int i, j;
        double val = 1.0;

        fgets(line, sizeof(line), file);

        char* p = line;
        int i = fast_atoi(p); // Utilise ton pointeur p qui avance
        int j = fast_atoi(p);
        if (!mm_is_pattern(matcode)) {
            val = fast_atof(p);
        }

        // if (mm_is_pattern(matcode)) {
        //     // Fichier Pattern : on ne lit que deux entiers
        //     fscanf(file, "%d %d\n", &i, &j);
        // } else {
        //     // Fichier Real/Integer : on lit les trois colonnes
        //     fscanf(file, "%d %d %lg\n", &i, &j, &val);
        // }
    
        int src = i - 1;
        int dest = j - 1;

        Element* elem = new Element;
        elem->row_index = src;
        elem->value = val;
        elem->next = matrice->columns[dest];
        matrice->columns[dest] = elem;

        degree_Sort_Tot[src] += val;
    }
    #pragma omp parallel for
    for (int j=0; j<matrice->n; j++) {
        Element* curr = matrice->columns[j];
        while (curr!= nullptr)
        {
            int i = curr->row_index;
            curr->value = curr->value/degree_Sort_Tot[i];
            curr = curr->next;
        }
    }

    // delete[] degree_Sort_Tot;
    fclose(file);
    return matrice;

}

void print_matrice(Matrice* matrice){
    for (int j = 0; j < 20; j++) {
        Element* curr = matrice->columns[j];
        if (!curr) std::cout << "vide";
        
        while (curr != nullptr) {
            std::cout << curr->row_index + 1 
                      << " -> " << j + 1 << " : " << std::fixed << setprecision(2) << curr->value << endl;
            curr = curr->next;
        }
    }

}

void print_vecteur(double* v, int n){
    cout << "pi = (";
    for(int i=0; i<n; i++){
        cout << v[i];
        if(i!=n-1) cout<< ", ";
    }
    cout << ")" << endl;
}

void print_vecteur(int* v, int n){
    cout << "f = (";
    for(int i=0; i<n; i++){
        cout << v[i];
        if(i!=n-1) cout<< ", ";
    }
    cout << ")" << endl;
}

void multiplier_vecteur_matrice(Matrice* matrice, double* vecteur, double* result, double alpha){
    int N = matrice->n;
    for (int j = 0; j<N; j++){
        double somme =0.0;
        Element* curr = matrice->columns[j];
        while (curr!=nullptr)
        {
           int i = curr->row_index;
           somme += vecteur[i]*curr->value;
           curr = curr->next;
        }
        result[j] =  alpha*somme;
    }
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

double calcul_G(Matrice* matrice, double* vecteur, double* result, int* f, double prem_term, double sec_term, double alpha){
    int N = matrice->n;

    double xf = multiplier_vecteurs(vecteur, f, N);
    double partie_constante = (prem_term*xf) + sec_term;

    double norme_totale = 0.0;

    #pragma omp parallel for reduction(+:norme_totale) schedule(dynamic, 512)
    for (int j = 0; j<N; j++){
        double somme =0.0;
        Element* curr = matrice->columns[j];
        while (curr!=nullptr)
        {
           somme += vecteur[curr->row_index]*curr->value;
           curr = curr->next;
        }
        double nouveau_score = (alpha*somme) + partie_constante;

        double diff = nouveau_score - vecteur[j];
        norme_totale += (diff < 0) ? -diff : diff;

        result[j] = nouveau_score;
    }

    return norme_totale;
}

double calcul_norme(double* v1, double* v2, int n){
    double somme = 0.0;

    #pragma omp parallel for reduction(+:somme)
    for (int i = 0; i < n; i++) {
        double ecart = v1[i] - v2[i];
        if (ecart < 0) ecart = -ecart;
        somme += ecart;
    }
    return somme;
}

// double calcul_norme2(double* v1, double* v2, int n) {
//     double somme_carres = 0.0;
//     for (int i = 0; i < n; i++) {
//         double diff = v1[i] - v2[i];
//         somme_carres += diff * diff;
//     }
//     return sqrt(somme_carres);
// }

void calcul_pagerank_classique(Matrice* P, double epsilon, double alpha){
    int N= P->n;

    double* pi_pair = new double[N];
    double* pi_impair = new double[N];
    
    for (int i = 0; i < N; i++) pi_pair[i] = 1.0/N;
    double norme = 1.0;
    int cpt = 0;

    while (norme > epsilon)
    {
        if(cpt % 2 == 0){
            multiplier_vecteur_matrice(P, pi_pair, pi_impair, alpha);
        }
        else{
            multiplier_vecteur_matrice(P, pi_impair, pi_pair, alpha);
        }

        // if(cpt% 7 == 0) 
        norme = calcul_norme(pi_impair, pi_pair, N);
        cpt++;

    }

    cout << "Nombre d'iterations : " << cpt << " | Norme : " << norme << endl;
    // print_vecteur(pi_impair, N);

    delete[] pi_impair;
    delete[] pi_pair;
    
}

void calcul_pagerank_G(Matrice* P, int* f, double epsilon, double alpha, int max_iter, int max_temps){
    int N= P->n;

    double* pi_pair = new double[N];
    double* pi_impair = new double[N];

    double prem_term = alpha/N;
    double sec_term = (1.0-alpha)/N;

    #pragma omp parallel for
    for (int i = 0; i < N; i++) pi_pair[i] = 1.0/N;
    double norme = 1.0;
    double ancienne_norme = 1.0;
    int cpt = 0;

    string raison_stop = "Maximum d'itérations atteint";

    auto start = chrono::steady_clock::now();

    while (cpt <max_iter)
    {
        if(cpt % 2 == 0){
            norme = calcul_G(P, pi_pair, pi_impair, f, prem_term, sec_term, alpha);
        }
        else{
            norme = calcul_G(P, pi_impair, pi_pair, f, prem_term, sec_term, alpha);
        }

        // norme = calcul_norme(pi_impair, pi_pair, N);
        if (norme < epsilon) {
            raison_stop = "Convergence atteinte (Epsilon)";
            cpt++;
            break;
        }

        if (abs(ancienne_norme - norme) < (epsilon * 0.01)) {
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

        cpt++;

    }
    auto end = chrono::steady_clock::now();
    chrono::duration<double, milli> total_ms = end - start;

    cout << "Nombre d'iterations : " << cpt << endl;
    cout << "Statut : " << raison_stop << endl;
    cout << "Temps : " << total_ms.count() << " ms" << endl;
    cout << "Norme : " << scientific << setprecision(10) << norme << endl;
    // print_vecteur(pi_pair, N);
    // print_vecteur(pi_impair, N);
    // print_vecteur(f, N);

    delete[] pi_impair;
    delete[] pi_pair;
    
}

void delete_matrix(Matrice* matrix) {
    if (!matrix) return;

    for (int i = 0; i < matrix->n; i++) {
        Element* curr = matrix->columns[i];
        while (curr != nullptr) {
            Element* toDelete = curr;
            curr = curr->next;
            delete toDelete; // On supprime chaque maillon un par un
        }
    }

    delete[] matrix->columns; // On supprime le tableau de pointeurs
    // delete[] matrix->out_weights; // Si tu as gardé le tableau de poids
    delete matrix; // On supprime la structure elle-même
    
}

int main(int argc, char** argv){

    double epsilon = 10e-8;
    double alpha = 0.85;
    int max_iter = 200, max_temps = 30;
    if (argc < 2) {
        std::cerr << "Usage : " << argv[0] << " <fichier_matrice.mtx>" << std::endl;
        std::cerr << "Exemple : ./pagerank web-Stanford.mtx" << std::endl;
        return 1; // On quitte avec un code d'erreur
    }

    const char* filename = argv[1];

    if (argc >= 3) max_iter  = atoi(argv[2]);
    if (argc >= 4) max_temps = atoi(argv[3]);

    std::cout << "--- Initialisationeee du PageRank ---" << std::endl;
    std::cout << "Fichier cible : " << filename << std::endl;
    int* degree_Sort_Tot;
    auto start_load = chrono::steady_clock::now();
    
    Matrice* myMatrix = read_mtx(filename, degree_Sort_Tot);
    
    auto end_load = chrono::steady_clock::now();

    if (myMatrix == nullptr) {
        return 1;
    }

    chrono::duration<double> elapsed_load = end_load - start_load;

    cout << "Matrice chargee avec succes." << endl;
    cout << "Nombre de pages (N) : " << myMatrix->n << endl;
    cout << "Nombre de liens (nz) : " << myMatrix->n_links << endl;
    cout << fixed << setprecision(2);
    cout << "==> Temps de chargement : " << elapsed_load.count() << " secondes" << endl;
    cout << "-----------------------------------" << endl;

    // if (myMatrix->n < 20) {
    //     print_matrice(myMatrix);
    // } else {
    //     std::cout << "Matrice chargee (" << myMatrix->n 
    //               << " pages). Trop grande pour un affichage complet." << std::endl;
    // }

    // calcul_pagerank_classique(myMatrix, epsilon, alpha);
    calcul_pagerank_G(myMatrix, degree_Sort_Tot, epsilon, alpha, max_iter, max_temps);
    delete_matrix(myMatrix);

    return 0;
}