#include <cmath>
#include <omp.h>
#include <iomanip>
#include <limits>
#include <algorithm>
#include <vector>
#include <fstream>
#include <string>
#include <iostream>



#include "matrice.cpp"
#include "ranking.cpp"
#include "test0.cpp"
#include "test.cpp"
#include "test2.cpp"
#include "writeResults.cpp"


using namespace std;

char* matriceName;

void printProgressBar(int current, int total, int width = 50) {
    float progress = (float)current / total;
    int pos = width * progress;

    std::cout << "[";
    for (int i = 0; i < width; ++i) {
        if (i < pos) std::cout << "=";
        else if (i == pos) std::cout << ">";
        else std::cout << " ";
    }
    
    // On utilise \r pour revenir au début de la ligne
    std::cout << "] " << int(progress * 100.0) << " %\r";
    
    // Important : il faut vider le tampon pour que l'affichage soit immédiat
    std::cout.flush();
}

int main(int argc, char** argv) {

    vector<double> epsilons = {1e-6, 1e-7, 1e-8, 1e-9, 1e-10, 1e-11, 1e-12};
    vector<double> alphas = {0.9, 0.85, 0.8, 0.75, 0.7, 0.65, 0.6, 0.55, 0.5, 0.45, 0.4, 0.35, 0.3, 0.25, 0.2, 0.15, 0.1};
    vector<int> max_temps_values = {100, 50, 25, 10, 5, 2, 1};

    // int max_temps = 300;
    int max_iter = 1000;

    // double alpha = 0.85;
    // double epsilon = 1e-6;



     if (argc < 2) {
        cerr << "Usage : " << argv[0] << " <fichier_matrice.mtx>" << endl;
        cerr << "Exemple : ./pagerank web-Stanford.mtx" << endl;
        return 1;
    }

    matriceName = argv[1];
    Matrice* myMatrice = read_mtx(matriceName);

    // if (argc >= 3) max_temps = atoi(argv[2]);

    int step = 0;
    int total_steps = epsilons.size() * alphas.size() * max_temps_values.size();

    printProgressBar(step, total_steps);

    for (int max_temps : max_temps_values) {
        

        for (const auto& alpha : alphas) {
            for (const auto& epsilon : epsilons) {
                Result* result1 = calcul_pagerank(myMatrice, epsilon, alpha, max_iter, max_temps);
                Result* result2 = calcul_pagerank1(myMatrice, epsilon, alpha, max_iter, max_temps);
                Result* result3 = calcul_pagerank2(myMatrice, epsilon, alpha, max_iter, max_temps);
                logExperiment(
                    "results.csv",
                    matriceName,
                    max_temps,
                    alpha,
                    epsilon,
                    result1,
                    result2,
                    result3
                );
                delete result1;
                delete result2;
                delete result3;
                step++;
            }
            printProgressBar(step, total_steps);
        }

        
    }

    // Result* result1 = calcul_pagerank(myMatrice, epsilon, alpha, max_iter, max_temps);
    // Result* result2 = calcul_pagerank1(myMatrice, epsilon, alpha, max_iter, max_temps);
    // Result* result3 = calcul_pagerank2(myMatrice, epsilon, alpha, max_iter, max_temps);

    // cout << result1->iterations << " - " << result1->raison << " - " << result1->time_ms << endl;
    // cout << result2->iterations << " - " << result2->raison << " - " << result2->time_ms << endl;
    // cout << result3->iterations << " - " << result3->raison << " - " << result3->time_ms << endl;
    

    delete_matrix(myMatrice);
    return 0;

}