void logExperiment(const std::string& filename,
                   const std::string& matrixName,
                   int max_temps,
                   double alpha, double epsilon,
                   Result* result1, Result* result2, Result* result3) {

    // Ouverture en mode 'append' (std::ios::app)
    std::ofstream file(filename, std::ios::app);

    if (!file.is_open()) {
        std::cerr << "Erreur : Impossible d'ouvrir le log " << filename << std::endl;
        return;
    }

    // Si le fichier est vide, on écrit les en-têtes (facultatif mais propre)
    // On peut vérifier la taille du fichier pour savoir si on écrit l'entête.
    file.seekp(0, std::ios::end);
    if (file.tellp() == 0) {
        file << "Matrice;Max_temps;Alpha;Epsilon;Id1;Iterations_1;Temps_1;Raison_1;Id2;Iterations_2;Temps_2;Raison_2;Id3;Iterations_3;Temps_3;Raison_3\n";
    }

    // Écriture des données sur une seule ligne
    file << matrixName << ";"
         << max_temps << ";"
         << alpha << ";"
         << epsilon << ";"
         << result1->id << ";"
         << result1->iterations << ";"
         << result1->time_ms << ";"
         << result1->raison << ";"
         << result2->id << ";"
         << result2->iterations << ";"
         << result2->time_ms << ";"
         << result2->raison << ";"
         << result3->id << ";"
         << result3->iterations << ";"
         << result3->time_ms << ";"
         << result3->raison << "\n";


    file.close();

}