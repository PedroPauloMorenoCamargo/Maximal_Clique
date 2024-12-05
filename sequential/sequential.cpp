#include <iostream>
#include <vector>
#include <fstream>
#include <string>
#include <algorithm>
#include <chrono>

using namespace std;

// Função para ler o grafo a partir do arquivo de entrada
vector<vector<int>> LerGrafo(const string& nomeArquivo, int& numVertices) {
    ifstream arquivo(nomeArquivo);
    int numArestas;
    arquivo >> numVertices >> numArestas;

    vector<vector<int>> grafo(numVertices, vector<int>(numVertices, 0));

    for (int i = 0; i < numArestas; ++i) {
        int u, v;
        arquivo >> u >> v;
        grafo[u - 1][v - 1] = 1;
        grafo[v - 1][u - 1] = 1;  // O grafo é não direcionado
    }

    arquivo.close();

    return grafo;
}

// Função para verificar se um conjunto de vértices forma uma clique
bool isClique(const vector<int>& vertices, const vector<vector<int>>& grafo) {
    for (size_t i = 0; i < vertices.size(); ++i) {
        for (size_t j = i + 1; j < vertices.size(); ++j) {
            if (grafo[vertices[i]][vertices[j]] == 0) {
                return false;
            }
        }
    }
    return true;
}

// Função para encontrar a maior clique no grafo
vector<int> encontrarCliqueMaxima(const vector<vector<int>>& grafo, int numVertices) {
    int maxSize = 0;
    vector<int> cliqueMaxima;

    for (int subset = 1; subset < (1 << numVertices); ++subset) {
        vector<int> clique;
        for (int i = 0; i < numVertices; ++i) {
            if (subset & (1 << i)) {
                clique.push_back(i);
            }
        }

        if (isClique(clique, grafo) && clique.size() > maxSize) {
            maxSize = clique.size();
            cliqueMaxima = clique;
        }
    }

    return cliqueMaxima;
}

int main(int argc, char* argv[]) {
    // Recebe arquivo via terminal
    if (argc < 2) {
        cerr << "Uso: " << argv[0] << " <nome_do_arquivo>" << endl;
        return 1;
    }

    // Nome do Arquivo
    string nomeArquivo = argv[1];
    // Número de Vértices do grafo
    int numVertices;

    try {
        vector<vector<int>> grafo = LerGrafo(nomeArquivo, numVertices);

        // Inicia o cronômetro
        auto inicio = chrono::high_resolution_clock::now();

        // Encontra a clique máxima no grafo
        vector<int> cliqueMaxima = encontrarCliqueMaxima(grafo, numVertices);

        // Para o cronômetro
        auto fim = chrono::high_resolution_clock::now();
        chrono::duration<double> duracao = fim - inicio;

        // Exibe o resultado no formato solicitado
        cout << "Clique máxima encontrada: [";
        for (size_t i = 0; i < cliqueMaxima.size(); ++i) {
            cout << "'" << (cliqueMaxima[i] + 1) << "'";
            if (i != cliqueMaxima.size() - 1) {
                cout << ", ";
            }
        }
        cout << "]" << endl;

        // Exibe o tempo de execução
        cout << "Tempo de execução: " << duracao.count() << " segundos" << endl;

    } catch (const exception& e) {
        cerr << "Erro: " << e.what() << endl;
        return 1;
    }

    return 0;
}


