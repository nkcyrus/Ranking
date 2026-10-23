#include <iostream>
#include <vector>
#include <fstream>
#include <algorithm>
#include <cmath>
#include <iomanip>
#include <sstream>

struct Matrice {
    int n;
    std::vector<long long> col_ptr; // CSR : indices de début de colonne
    std::vector<int> row_idx;       // CSR : indices de ligne (sources)
    std::vector<int> out_degree;    // Nombre de liens sortants par nœud
};

// Chargeur de matrice au format Matrix Market
Matrice charger_mtx(const std::string& filename) {
    std::ifstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Erreur : " << filename << " introuvable." << std::endl;
        exit(1);
    }
    std::string line;
    while (std::getline(file, line) && line[0] == '%');
    int r, c; long long nz;
    std::stringstream ss(line); ss >> r >> c >> nz;
    
    Matrice m; m.n = c;
    m.out_degree.assign(m.n, 0);
    m.col_ptr.assign(m.n + 1, 0);
    
    struct Edge { int u, v; };
    std::vector<Edge> edges;
    for(long long i = 0; i < nz; ++i) {
        int u, v; file >> u >> v; u--; v--;
        if (u >= 0 && u < m.n && v >= 0 && v < m.n) {
            edges.push_back({u, v});
            m.out_degree[u]++;
            m.col_ptr[v + 1]++;
        }
    }
    // Construction de la structure de stockage par colonnes (pour multiplier par G)
    for (int i = 0; i < m.n; i++) m.col_ptr[i+1] += m.col_ptr[i];
    m.row_idx.resize(edges.size());
    std::vector<long long> current_pos = m.col_ptr;
    for (auto& e : edges) m.row_idx[current_pos[e.v]++] = e.u;
    return m;
}

int main() {
    Matrice m = charger_mtx("harvard500.mtx");
    int N = m.n;
    double alpha = 0.85;
    double invN = 1.0 / N;
    double epsilon = 1e-10;

    // 1. X(0) = nabla
    std::vector<double> nabla(N, (1.0 - alpha) * invN);
    // 2. Y(0) = delta (initialisé avec le max possible par ligne de G)
    std::vector<double> delta(N, (1.0 - alpha) * invN);
    for (int j = 0; j < N; j++) {
        double max_val = 0;
        for (long long i = m.col_ptr[j]; i < m.col_ptr[j+1]; i++) {
            max_val = std::max(max_val, alpha / (double)m.out_degree[m.row_idx[i]]);
        }
        delta[j] += max_val;
    }

    std::vector<double> X = nabla, Y = delta;
    std::vector<double> Xnew(N), Ynew(N);

    std::cout << "Demarrage Squeeze Strict (ton algorithme)..." << std::endl;

    for (int k = 0; k < 500; k++) {
        double normX = 0, normY = 0;
        double danglingX = 0, danglingY = 0;

        for (int i = 0; i < N; i++) {
            normX += X[i];
            normY += Y[i];
            if (m.out_degree[i] == 0) {
                danglingX += X[i];
                danglingY += Y[i];
            }
        }

        double diffL1 = 0;
        for (int j = 0; j < N; j++) {
            // Calcul du produit vecteur-matrice : X(k) * G
            // G_ij = alpha * P_ij + (1-alpha)/N
            // Pour les dangling nodes : P_ij = 1/N
            double XG_j = 0, YG_j = 0;
            
            // Partie alpha * P (liens explicites)
            for (long long i = m.col_ptr[j]; i < m.col_ptr[j+1]; i++) {
                int u = m.row_idx[i];
                XG_j += X[u] * (alpha / m.out_degree[u]);
                YG_j += Y[u] * (alpha / m.out_degree[u]);
            }
            
            // Partie Dangling + Saut Aleatoire (le reste de G)
            double commonX = (alpha * danglingX + (1.0 - alpha)) * invN;
            double commonY = (alpha * danglingY + (1.0 - alpha)) * invN;
            XG_j += commonX;
            YG_j += commonY;

            // (a) X(k+1) = max(X(k), X(k)G + nabla(1 - ||X(k)||1))
            Xnew[j] = std::max(X[j], XG_j + nabla[j] * (1.0 - normX));

            // (b) Y(k+1) = min(Y(k), Y(k)G + delta(1 - ||Y(k)||1))
            Ynew[j] = std::min(Y[j], YG_j + delta[j] * (1.0 - normY));

            diffL1 += std::abs(Ynew[j] - Xnew[j]);
        }

        X = Xnew; 
        Y = Ynew;

        if (k % 10 == 0 || diffL1 < epsilon) {
            std::cout << "Iter " << std::setw(3) << k 
                      << " | L1: " << std::scientific << std::setprecision(4) << diffL1 
                      << " | NormX: " << std::fixed << std::setprecision(4) << normX << std::endl;
        }

        if (diffL1 < epsilon) break;
    }

    return 0;
}