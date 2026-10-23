#include "ranking.hpp"
#include <iostream>
#include <algorithm>
#include <iomanip>


#include "writeTop.cpp"

const int top = 20;

bool comparerPages(const Page& a, const Page& b) {
    return a.score > b.score;
}

void extraireTop10(double* pi, int N, long long id) {
    std::vector<Page> topPages(N);

    // Initialisation parallèle du vecteur de tri
    #pragma omp parallel for
    for (int i = 0; i < N; ++i) {
        topPages[i].id = i;
        topPages[i].score = pi[i];
    }

    std::nth_element(topPages.begin(), topPages.begin() + top, topPages.end(), comparerPages);

    std::sort(topPages.begin(), topPages.begin() + top, comparerPages);

    writeTop10("classements.csv", id, topPages, top);

}

void extraireTop10(double* x, double* y, int N, long long id) {
    std::vector<Page> topPages(N);

    // Initialisation parallèle du vecteur de tri
    #pragma omp parallel for
    for (int i = 0; i < N; ++i) {
        topPages[i].id = i;
        topPages[i].score = (x[i] + y[i]) * 0.5;
    }

    std::nth_element(topPages.begin(), topPages.begin() + top, topPages.end(), comparerPages);

    std::sort(topPages.begin(), topPages.begin() + top, comparerPages);

    writeTop10("classements.csv", id, topPages, top);

}

void extraireTop10_parallel(double* pi, int N, long long id) {
    std::vector<Page> topPages(N);

    // Initialisation parallèle du vecteur de tri
    int i;
    #pragma omp parallel for
    for (i = 0; i < N; ++i) {
        int j = i*2;
        topPages[i].id = i;
        topPages[i].score = (pi[j] + pi[j+1]) * 0.5;
    }

    std::nth_element(topPages.begin(), topPages.begin() + top, topPages.end(), comparerPages);

    std::sort(topPages.begin(), topPages.begin() + top, comparerPages);

    writeTop10("classements.csv", id, topPages, top);
    
}