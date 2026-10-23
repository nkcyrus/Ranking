#include "structs.cpp"
extern "C" {
    #include "mmio.c"
}

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

Matrice* read_mtx(char* filename){
    FILE* file = fopen(filename, "r");
    if (!file) return nullptr;
 
    char sys_buffer[1 << 20];
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
 
    int* temp_rows = new int[nz];
    int* temp_cols = new int[nz];
    bool vals = !mm_is_pattern(matcode);
    double* temp_vals = nullptr;
    if (vals) temp_vals = new double[nz];
 
    char line[MM_MAX_LINE_LENGTH];
    for (long long i = 0; i < nz; i++) {
        if (!fgets(line, sizeof(line), file)) break;
        char* p = line;
        int r = fast_atoi(p);
        int c = fast_atoi(p);
        temp_rows[i] = r - 1;
        temp_cols[i] = c - 1;
        nb_liens_par_colonne[c - 1]++;

        if (vals) {
            temp_vals[i] = fast_atof(p);
            matrice->degree_Sort_Tot[r - 1] += temp_vals[i];
        }
        else matrice->degree_Sort_Tot[r - 1]++;
    }
    fclose(file);
 
    // Construction des pointeurs CSC
    matrice->col_ptr = new long long[N + 1]();
    for (int j = 0; j < N; j++)
        matrice->col_ptr[j + 1] = matrice->col_ptr[j] + nb_liens_par_colonne[j];
 
    matrice->row_idx = new int[nz];
    matrice->values = new double[nz]();
 
    long long* current_pos = new long long[N]();
    for (long long i = 0; i < nz; i++) {
        int col = temp_cols[i];
        int row = temp_rows[i];
        long long pos = matrice->col_ptr[col] + current_pos[col];
        matrice->row_idx[pos] = row;
        int deg = matrice->degree_Sort_Tot[row];

        if (vals){
            if (deg !=0 ) matrice->values[pos] = temp_vals[i] / (double) deg;
        }
        else matrice->values[pos] = 1.0 / (double) deg;
        current_pos[col]++;
    }
 
    delete[] temp_rows;
    delete[] temp_cols;
    if (vals) delete[] temp_vals;
    delete[] current_pos;
    delete[] nb_liens_par_colonne;
 
    return matrice;
}

void delete_matrix(Matrice* matrice) {
    if (matrice == nullptr) return;

    delete[] matrice->col_ptr;
    delete[] matrice->row_idx;
    delete[] matrice->values;
    delete[] matrice->degree_Sort_Tot;

    delete matrice;
    
}