#include <iostream>
#include <vector>
#include <algorithm>
#include <fstream>
#include <chrono>
#include <random>
#include <omp.h>

using namespace std;
using namespace std::chrono;

vector<vector<int>> LerGrafo(const string& nomeArquivo, int& numVertices);
vector<int> ConstrucaoClique(const vector<vector<int>>& graph, const vector<int>& candidates);

int main(int argc, char* argv[]) {
    if (argc < 3) {
        cerr << "Uso: " << argv[0] << " <nome_do_arquivo> <numero_de_iteracoes>" << endl;
        return 1;
    }

    // Recebe o nome do arquivo da linha de comando
    string nomeArquivo = argv[1]; 
    // Número de iterações para a busca local
    int numIteracoes = stoi(argv[2]); 
    
    int numVertices;
    // Lê o grafo a partir do arquivo fornecido
    vector<vector<int>> graph = LerGrafo(nomeArquivo, numVertices);

    // Calcula o grau de cada vértice
    vector<pair<int, int>> degreeVertexPairs;
    for (int i = 0; i < numVertices; ++i) {
        int degree = 0;
        for (int j = 0; j < numVertices; ++j) {
            degree += graph[i][j];
        }
        degreeVertexPairs.push_back({degree, i});
    }

    // Ordena os vértices por grau decrescente
    sort(degreeVertexPairs.rbegin(), degreeVertexPairs.rend());

    // Prepara a lista de candidatos inicial
    vector<int> candidates;
    for (const auto &pair : degreeVertexPairs) {
        candidates.push_back(pair.second);
    }

    // Inicia a contagem de tempo
    auto start_time = high_resolution_clock::now();

    // Constrói a clique usando a ordem inicial (primeira iteração)
    vector<int> bestClique = ConstrucaoClique(graph, candidates);

    // Preparação para o paralelismo
    // Vamos criar um gerador de números aleatórios por thread
    int nThreads;
    #pragma omp parallel
    {
        #pragma omp single
        nThreads = omp_get_num_threads();
    }

    vector<mt19937> gens(nThreads);
    {
        // Seed global para cada thread
        random_device rd;
        for (int i = 0; i < nThreads; i++) {
            gens[i].seed(rd());
        }
    }

    // Execução das iterações de hill climbing em paralelo
    #pragma omp parallel
    {
        int tid = omp_get_thread_num();
        mt19937 &gen = gens[tid];
        uniform_int_distribution<> dist(0, (int)candidates.size()-1);

        #pragma omp for
        for (int it = 0; it < numIteracoes; ++it) {
            // Cria uma cópia do vetor de candidatos para perturbar
            vector<int> perturbedCandidates = candidates;
            
            // Perturbação simples: troca a posição de dois vértices aleatórios
            int idx1 = dist(gen);
            int idx2 = dist(gen);
            while (idx2 == idx1) {
                idx2 = dist(gen);
            }
            swap(perturbedCandidates[idx1], perturbedCandidates[idx2]);

            // Constrói a clique com a nova ordenação
            vector<int> newClique = ConstrucaoClique(graph, perturbedCandidates);

            // Se melhorou, atualiza a melhor solução (região crítica)
            #pragma omp critical
            {
                if ((int)newClique.size() > (int)bestClique.size()) {
                    bestClique = newClique;
                    candidates = perturbedCandidates;
                }
            }
        }
    }

    // Finaliza a contagem de tempo
    auto end_time = high_resolution_clock::now();
    auto duration = duration_cast<std::chrono::duration<double>>(end_time - start_time);

    // Exibe o tamanho da melhor clique encontrada e os vértices
    cout << "Tamanho da melhor clique encontrada: " << bestClique.size() << endl;
    cout << "Vértices na Clique: ";
    for (int vertex : bestClique) {
        cout << vertex + 1 << " ";
    }
    cout << endl;


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

vector<int> ConstrucaoClique(const vector<vector<int>>& graph, const vector<int>& candidates) {
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
    return clique;
}
