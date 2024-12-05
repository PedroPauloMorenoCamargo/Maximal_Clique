#include <iostream>
#include <vector>
#include <algorithm>
#include <fstream>
#include <chrono>

using namespace std;
using namespace std::chrono;

// Function prototypes
vector<vector<int>> LerGrafo(const string& nomeArquivo, int& numVertices);

int main(int argc, char* argv[]) {
    if (argc < 2) {
        cerr << "Uso: " << argv[0] << " <nome_do_arquivo>" << endl;
        return 1;
    }

    string nomeArquivo = argv[1]; // Recebe o nome do arquivo da linha de comando
    int numVertices;

    // Lê o grafo a partir do arquivo fornecido
    vector<vector<int>> graph = LerGrafo(nomeArquivo, numVertices);

    // Calcula o grau de cada vértice
    vector<pair<int, int>> degreeVertexPairs; // Pair of (degree, vertex)
    for (int i = 0; i < numVertices; ++i) {
        int degree = 0;
        for (int j = 0; j < numVertices; ++j) {
            degree += graph[i][j];
        }
        degreeVertexPairs.push_back({degree, i});
    }

    // Ordena os vértices por grau decrescente
    sort(degreeVertexPairs.rbegin(), degreeVertexPairs.rend());

    // Prepara a lista de candidatos
    vector<int> candidates;
    for (const auto &pair : degreeVertexPairs) {
        candidates.push_back(pair.second);
    }

    // Inicia a contagem de tempo
    auto start_time = high_resolution_clock::now();

    // Implementa a heurística gulosa para encontrar uma clique
    vector<int> clique;
    for (int v : candidates) {
        bool canAdd = true;
        for (int u : clique) {
            if (graph[u][v] == 0) {
                canAdd = false;
                break;
            }
        }
        if (canAdd) {
            clique.push_back(v);
        }
    }

    // Finaliza a contagem de tempo
    auto end_time = high_resolution_clock::now();
    auto duration = duration_cast<std::chrono::duration<double>>(end_time - start_time);

    // Exibe o tamanho da clique encontrada e os vértices
    cout << "Tamanho da Clique encontrada pela heurística gulosa: " << clique.size() << endl;
    cout << "Vértices na Clique: ";
    for (int vertex : clique) {
        cout << vertex + 1 << " "; // Ajusta o índice para ser 1-based
    }
    cout << endl;

    // Exibe o tempo gasto
    cout << "Tempo gasto para encontrar a clique: " << duration.count() << " segundos" << endl;

    return 0;
}

vector<vector<int>> LerGrafo(const string& nomeArquivo, int& numVertices) {
    ifstream arquivo(nomeArquivo);
    if (!arquivo.is_open()) {
        throw runtime_error("Erro ao abrir o arquivo: " + nomeArquivo);
    }

    int numArestas;
    arquivo >> numVertices >> numArestas;

    vector<vector<int>> grafo(numVertices, vector<int>(numVertices, 0));

    for (int i = 0; i < numArestas; ++i) {
        int u, v;
        arquivo >> u >> v;
        grafo[u - 1][v - 1] = 1;
        grafo[v - 1][u - 1] = 1; 
    }

    arquivo.close();

    return grafo;
}
