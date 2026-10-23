void writeTop10(const std::string& filename, long long id, const std::vector<Page>& topPages, int top) {
    std::ofstream file(filename, std::ios::app);

    if (!file.is_open()) {
        std::cerr << "Erreur : Impossible d'ouvrir le log " << filename << std::endl;
        return;
    }

    file.seekp(0, std::ios::end);
    if (file.tellp() == 0) {
        file << "ID;Position;Page;Score\n";

    }

    for (int i = 0; i < top; ++i) {
            file << id << ";"
                 << i+1 << ";"
                 << topPages[i].id + 1 << ";"
                 << topPages[i].score << "\n";
        }

        


    file.close();
}