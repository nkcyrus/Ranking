#include <string>
struct Matrice{
    int n;
    long long n_links;
    double* values;
    long long* col_ptr;
    int* row_idx;
    int* degree_Sort_Tot;
    double* nabla;
    double* delta;
    double* nd;
};

struct Result {
    long long id;
    int iterations;
    double time_ms;
    std::string raison;

};

struct Page {
    int id;
    double score;
};