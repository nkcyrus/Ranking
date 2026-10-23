/*
 * Distribution stationnaire d'une chaîne de Markov — PageRank
 *
 * Charge une matrice d'adjacence creuse au format Matrix Market (.mtx),
 * construit la matrice de Google :
 *
 *   G = α·P + α·(1/N)·f·eᵀ + (1−α)·(1/N)·e·eᵀ
 *
 * où :
 *   P   = matrice stochastique par ligne (normalisée depuis l'adjacence brute)
 *   f   = vecteur dangling nodes : f[i] = 1 si la ligne i de P est nulle, 0 sinon
 *   e   = vecteur colonne de 1
 *   α   = facteur d'amortissement (0.85 par défaut)
 *   N   = nombre de nœuds
 *
 * Puis applique l'algorithme d'encadrement X(k)/Y(k) → distribution stationnaire π.
 *
 * Usage :
 *   ./markov_stationary <fichier.mtx> [alpha=0.85] [eps=1e-10] [maxIter=10000]
 *
 * Compilation :
 *   g++ -O2 -std=c++17 -o markov_stationary markov_stationary.cpp
 */

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <cmath>
#include <limits>
#include <numeric>
#include <stdexcept>
#include <iomanip>
#include <algorithm>

// ============================================================
//  Types
// ============================================================
using Matrix = std::vector<std::vector<double>>;
using Vector = std::vector<double>;

// ============================================================
//  Utilitaires vectoriels
// ============================================================
int nRows(const Matrix& M) { return (int)M.size(); }
int nCols(const Matrix& M) { return M.empty() ? 0 : (int)M[0].size(); }

double norm1(const Vector& v) {
    double s = 0.0;
    for (double x : v) s += std::abs(x);
    return s;
}

double norm1Diff(const Vector& a, const Vector& b) {
    double s = 0.0;
    for (int i = 0; i < (int)a.size(); i++) s += std::abs(a[i] - b[i]);
    return s;
}

Vector normalize(const Vector& v) {
    double n = norm1(v);
    if (n == 0.0) throw std::runtime_error("Normalisation par zero.");
    Vector r(v.size());
    for (int i = 0; i < (int)v.size(); i++) r[i] = v[i] / n;
    return r;
}

// Produit vecteur-ligne x matrice dense : (v*M)[j] = sum_i v[i]*M[i][j]
Vector vecMatMul(const Vector& v, const Matrix& M) {
    int rows = nRows(M), cols = nCols(M);
    if ((int)v.size() != rows)
        throw std::runtime_error("Dimensions incompatibles pour le produit v*M.");
    Vector r(cols, 0.0);
    for (int i = 0; i < rows; i++)
        if (v[i] != 0.0)
            for (int j = 0; j < cols; j++)
                r[j] += v[i] * M[i][j];
    return r;
}

Vector vecMax(const Vector& a, const Vector& b) {
    Vector r(a.size());
    for (int i = 0; i < (int)a.size(); i++) r[i] = std::max(a[i], b[i]);
    return r;
}

Vector vecMin(const Vector& a, const Vector& b) {
    Vector r(a.size());
    for (int i = 0; i < (int)a.size(); i++) r[i] = std::min(a[i], b[i]);
    return r;
}

void printVector(const std::string& label, const Vector& v, int maxShow = 10) {
    int n = (int)v.size();
    std::cout << label << " (dim=" << n << ") = [";
    int show = std::min(n, maxShow);
    for (int i = 0; i < show; i++) {
        std::cout << std::fixed << std::setprecision(8) << v[i];
        if (i + 1 < show) std::cout << ", ";
    }
    if (n > maxShow) std::cout << ", ...";
    std::cout << "]\n";
}

// ============================================================
//  Chargement Matrix Market creux (.mtx)
//  Format : lignes de commentaire (%), puis "rows cols nnz",
//           puis triplets "i j [val]" (1-indexes).
//  On construit une matrice dense N×N.
// ============================================================
Matrix loadMTX(const std::string& filename) {
    std::ifstream f(filename);
    if (!f.is_open())
        throw std::runtime_error("Impossible d'ouvrir : " + filename);

    std::string line;
    bool hasHeader  = false;
    bool isPattern  = false;   // format "pattern" -> valeurs toutes 1
    bool isSymmetric = false;  // format "symmetric" -> symetrie implicite

    // Lire la banniere %%MatrixMarket et les commentaires
    while (std::getline(f, line)) {
        if (line.empty()) continue;
        if (line[0] == '%') {
            std::string low = line;
            std::transform(low.begin(), low.end(), low.begin(), ::tolower);
            if (low.find("pattern")   != std::string::npos) isPattern   = true;
            if (low.find("symmetric") != std::string::npos) isSymmetric = true;
            continue;
        }
        // Premiere ligne non-commentaire : dimensions
        hasHeader = true;
        break;
    }
    if (!hasHeader)
        throw std::runtime_error("Fichier MTX invalide : en-tete manquant.");

    int rows, cols, nnz;
    {
        std::istringstream ss(line);
        ss >> rows >> cols >> nnz;
    }
    if (rows != cols)
        throw std::runtime_error("La matrice doit etre carree (N x N).");

    int N = rows;
    std::cout << "MTX charge : N=" << N << "  nnz=" << nnz;
    if (isPattern)   std::cout << "  [pattern]";
    if (isSymmetric) std::cout << "  [symmetric]";
    std::cout << "\n";

    Matrix A(N, Vector(N, 0.0));

    // Auto-detection du format pattern sur la premiere ligne de donnees :
    // si le parseur ne trouve pas de 3e token (valeur), on bascule en mode pattern.
    bool patternAutoDetected = false;
    int r, c; double val;
    while (std::getline(f, line)) {
        if (line.empty() || line[0] == '%') continue;
        std::istringstream ss(line);
        if (!(ss >> r >> c)) continue;   // ligne invalide, on saute
        if (isPattern) {
            val = 1.0;
        } else {
            if (!(ss >> val)) {
                // Pas de 3e colonne : format "i j" pur -> mode pattern
                val = 1.0;
                if (!patternAutoDetected) {
                    isPattern = patternAutoDetected = true;
                    std::cout << "  (format 'i j' detecte automatiquement : valeurs = 1)\n";
                }
            }
        }
        // MTX est 1-indexe
        A[r-1][c-1] = val;
        if (isSymmetric && r != c)
            A[c-1][r-1] = val;
    }

    return A;
}

// ============================================================
//  Construction de la matrice de Google
//
//  G = alpha*P + alpha*(1/N)*f*e^T + (1-alpha)*(1/N)*e*e^T
//
//  P[i][j] = A[i][j] / rowSum[i]  (si rowSum[i] > 0, sinon 0)
//  f[i]    = 1 si rowSum[i] == 0 (dangling node), 0 sinon
//
//  Developpe terme a terme :
//  G[i][j] = alpha*P[i][j] + alpha*f[i]/N + (1-alpha)/N
// ============================================================
Matrix buildGoogleMatrix(const Matrix& A, double alpha) {
    int N = nRows(A);
    double invN = 1.0 / N;

    // Sommes de lignes et identification des dangling nodes
    Vector rowSum(N, 0.0);
    for (int i = 0; i < N; i++)
        for (int j = 0; j < N; j++)
            rowSum[i] += A[i][j];

    Vector f(N, 0.0);
    int nDangling = 0;
    for (int i = 0; i < N; i++)
        if (rowSum[i] == 0.0) { f[i] = 1.0; nDangling++; }

    std::cout << "Dangling nodes : " << nDangling << " / " << N << "\n";

    // G[i][j] = alpha*P[i][j] + alpha*f[i]/N + (1-alpha)/N
    Matrix G(N, Vector(N));
    for (int i = 0; i < N; i++) {
        // Terme commun a toute la ligne i (independant de j)
        double base = alpha * f[i] * invN + (1.0 - alpha) * invN;
        for (int j = 0; j < N; j++) {
            double Pij = (rowSum[i] > 0.0) ? A[i][j] / rowSum[i] : 0.0;
            G[i][j] = alpha * Pij + base;
        }
    }

    return G;
}

// ============================================================
//  Algorithme d'encadrement -> distribution stationnaire
//
//  nabla[j] = min_i G[i][j]
//  delta[j] = max_i G[i][j]
//
//  X(0) = nabla / ||nabla||_1
//  Y(0) = delta / ||delta||_1
//
//  Iterations :
//    X <- max(X*G, nabla) / ||...||_1
//    Y <- min(Y*G, delta) / ||...||_1
//  jusqu'a ||X - Y||_1 < eps
// ============================================================
Vector stationaryDistribution(const Matrix& G, double eps, int maxIter) {
    int N = nCols(G);

    Vector nabla(N,  std::numeric_limits<double>::infinity());
    Vector delta(N, -std::numeric_limits<double>::infinity());
    for (int j = 0; j < N; j++)
        for (int i = 0; i < N; i++) {
            nabla[j] = std::min(nabla[j], G[i][j]);
            delta[j] = std::max(delta[j], G[i][j]);
        }

    Vector X = normalize(nabla);
    Vector Y = normalize(delta);

    std::cout << "\n=== Algorithme de distribution stationnaire ===\n";
    printVector("nabla (min col)", nabla);
    printVector("delta (max col)", delta);
    printVector("X(0)", X);
    printVector("Y(0)", Y);
    std::cout << "\n";

    for (int k = 0; k < maxIter; k++) {
        // (a) X <- max(X*G, nabla),  puis normaliser
        Vector Xnew = normalize(vecMax(vecMatMul(X, G), nabla));
        // (b) Y <- min(Y*G, delta),  puis normaliser
        Vector Ynew = normalize(vecMin(vecMatMul(Y, G), delta));

        X = Xnew;
        Y = Ynew;

        double diff = norm1Diff(X, Y);

        if (k < 10 || k % 100 == 0)
            std::cout << "k=" << std::setw(6) << k+1
                      << "  ||X-Y||_1 = " << std::scientific << std::setprecision(4) << diff << "\n";

        if (diff < eps) {
            std::cout << "\nConvergence a k=" << k+1
                      << "  (||X-Y||_1=" << diff << " < eps=" << eps << ")\n\n";
            break;
        }
        if (k == maxIter - 1)
            std::cout << "\nAttention : maxIter=" << maxIter << " atteint sans convergence.\n\n";
    }

    printVector("X final (borne inf)", X);
    printVector("Y final (borne sup)", Y);

    // Distribution stationnaire : milieu de l'encadrement
    Vector pi(N);
    for (int j = 0; j < N; j++) pi[j] = 0.5 * (X[j] + Y[j]);
    printVector("pi ~ (X+Y)/2      ", pi);

    return pi;
}

// ============================================================
//  Main
// ============================================================
int main(int argc, char* argv[]) {
    if (argc < 2) {
        std::cerr << "Usage : " << argv[0]
                  << " <fichier.mtx> [alpha=0.85] [eps=1e-10] [maxIter=10000]\n";
        return 1;
    }

    std::string mtxFile = argv[1];
    double alpha   = (argc >= 3) ? std::stod(argv[2]) : 0.85;
    double eps     = (argc >= 4) ? std::stod(argv[3]) : 1e-10;
    int    maxIter = (argc >= 5) ? std::stoi(argv[4]) : 10000;

    if (alpha <= 0.0 || alpha >= 1.0)
        throw std::runtime_error("alpha doit etre dans ]0, 1[.");

    std::cout << "=== PageRank / Distribution stationnaire ===\n"
              << "Fichier : " << mtxFile << "\n"
              << "alpha   = " << alpha   << "\n"
              << "eps     = " << eps     << "\n"
              << "maxIter = " << maxIter << "\n\n";

    // 1. Charger la matrice d'adjacence brute
    Matrix A = loadMTX(mtxFile);
    int N = nRows(A);

    // 2. Construire la matrice de Google G
    std::cout << "\nConstruction de G (N=" << N << ", alpha=" << alpha << ")...\n";
    Matrix G = buildGoogleMatrix(A, alpha);

    // 3. Calculer la distribution stationnaire
    Vector pi = stationaryDistribution(G, eps, maxIter);

    // 4. Verification : pi*G ~= pi
    Vector piG = vecMatMul(pi, G);
    std::cout << "=== Verification : pi*G ~= pi ===\n"
              << "||pi*G - pi||_1 = " << std::scientific << norm1Diff(pi, piG) << "\n";

    // 5. Top-10 noeuds par score PageRank
    int top = std::min(N, 10);
    std::cout << "\n=== Top-" << top << " noeuds (PageRank) ===\n";
    std::vector<int> idx(N);
    std::iota(idx.begin(), idx.end(), 0);
    std::partial_sort(idx.begin(), idx.begin() + top, idx.end(),
                      [&pi](int a, int b){ return pi[a] > pi[b]; });
    for (int k = 0; k < top; k++)
        std::cout << "  #" << std::setw(2) << k+1
                  << "  noeud " << std::setw(6) << idx[k]+1   // affiche en 1-indexe
                  << "  pi=" << std::fixed << std::setprecision(10) << pi[idx[k]] << "\n";

    return 0;
}
