#include <iostream>
#include <vector>
#include <algorithm>
#include <fstream>
#include <chrono>
#include <stdexcept>
#include <bitset>

using namespace std;
using namespace std::chrono;

/*Ajuste o valor de MAX_VERTICES conforme necessário para o tamanho máximo do grafo sendo que MAX_VERTICES deve ser múltiplo de 64.
Uma vez que iremos trabalhar com  bitsets. Além do mais, o valor de MAX_VERTICES deve ser maior ou igual ao número de vértices do grafo.*/
const int MAX_VERTICES = 2048; 
// Checa se MAX_VERTICES é múltiplo de 64
static_assert(MAX_VERTICES % 64 == 0, "MAX_VERTICES should be multiple of 64");

// O grafo é representado por um vetor de bitsets, onde cada bitset representa a adjacência de um vértice.
static bitset<MAX_VERTICES> adjacency[MAX_VERTICES];

// Variável global para armazenar o tamanho da maior clique encontrada
int bestCliqueSizeGlobal = 0; 

void findMaximumClique(const vector<int> &candidates, vector<int> &currentClique, vector<int> &bestClique);

vector<vector<int>> LerGrafo(const string& nomeArquivo, int &numVertices);

int main(int argc, char* argv[]) {
    // Verifica se o nome do arquivo foi fornecido
    if (argc < 2) {
        cerr << "Uso: " << argv[0] << " <nome_do_arquivo>" << endl;
        return 1;
    }
    // Recebe o nome do arquivo da linha de comando
    string nomeArquivo = argv[1];
    // Número de vértices no grafo
    int numVertices;

    // Lê o grafo
    vector<vector<int>> graph = LerGrafo(nomeArquivo, numVertices);
    // Verifica se o número de vértices não excede MAX_VERTICES
    if (numVertices > MAX_VERTICES) {
        cerr << "Number of vertices exceeds MAX_VERTICES. Increase MAX_VERTICES." << endl;
        return 1;
    }

    // Converte a matriz de adjacência para bitsets
    for (int i = 0; i < numVertices; i++) {
        // Zera o bitset
        adjacency[i].reset();
        for (int j = 0; j < numVertices; j++) {
            // Se houver uma aresta entre os vértices i e j, seta o bit correspondente
            if (graph[i][j] == 1) {
                adjacency[i].set(j);
            }
        }
    }

    // Vetor de pares (grau, vértice)
    vector<pair<int, int>> degreeVertexPairs;
    // Reserva espaço para o vetor
    degreeVertexPairs.reserve(numVertices);

    // Calcula o grau e ordena os vértices por grau decrescente
    for (int i = 0; i < numVertices; ++i) {
        // Conta o número de bits setados no bitset
        int degree = (int)adjacency[i].count();
        // Adiciona o par (grau, vértice) ao vetor
        degreeVertexPairs.push_back({degree, i});
    }

    // Ordena o vetor por grau decrescente
    sort(degreeVertexPairs.rbegin(), degreeVertexPairs.rend());

    // Vetor de candidatos
    vector<int> candidates;
    // Reserva espaço para o vetor
    candidates.reserve(numVertices);
    // Adiciona os vértices ordenados por grau decrescente
    for (auto &p : degreeVertexPairs) {
        candidates.push_back(p.second);
    }
    // Vetor para armazenar a clique atual
    vector<int> currentClique;
    // Reserva espaço para o vetor
    currentClique.reserve(numVertices);
    // Vetor para armazenar a maior clique encontrada
    vector<int> bestClique;
    // Reserva espaço para o vetor
    bestClique.reserve(numVertices);

    // Contagem de tempo
    auto start_time = high_resolution_clock::now();
    // Encontra a maior clique
    findMaximumClique(candidates, currentClique, bestClique);
    auto end_time = high_resolution_clock::now();

    auto duration = duration_cast<std::chrono::duration<double>>(end_time - start_time);

    cout << "Maximum Clique Size: " << bestClique.size() << endl;
    cout << "Vertices in the Maximum Clique: ";
    for (int vertex : bestClique) {
        cout << vertex + 1 << " ";
    }
    cout << endl;
    cout << "Time taken: " << duration.count() << " seconds" << endl;

    return 0;
}

void findMaximumClique(const vector<int> &candidates, vector<int> &currentClique, vector<int> &bestClique) {
    // Se não houver mais candidatos, a clique atual é a maior encontrada
    if (candidates.empty()) {
        // Atualiza a maior clique encontrada caso a clique atual seja maior
        if ((int)currentClique.size() > bestCliqueSizeGlobal) {
            bestCliqueSizeGlobal = (int)currentClique.size();
            bestClique = currentClique;
        }
        return;
    }

    // Podagem: Se o tamanho da clique atual mais o número de candidatos restantes não puder exceder o tamanho da maior clique global, retorna
    if ((int)currentClique.size() + (int)candidates.size() <= bestCliqueSizeGlobal) {
        return;
    }

    // Itera sobre os candidatos em ordem reversa
    for (int i = (int)candidates.size() - 1; i >= 0; --i) {
        // Vértice atual
        int v = candidates[i];

        // Segunda podagem: Se o tamanho da clique atual mais o número de candidatos restantes não puder exceder o tamanho da maior clique global, retorna
        if ((int)currentClique.size() + i + 1 <= bestCliqueSizeGlobal) {
            return;
        }

        // Checa se o vértice v pode ser adicionado à clique atual (é adjacente a todos os vértices da clique atual)
        for (int u : currentClique) {
            if (!adjacency[u].test(v)) {
                goto skip_vertex;
            }
        }

        // Adiciona o vértice v à clique atual
        currentClique.push_back(v);

        // Constrói novos candidatos por interseção dos candidatos anteriores com a adjacência de v e chama recursivamente a função
        {
            // Cria uma máscara de bits para a relação de adjacência de v com todos os vértices da clique atual
            bitset<MAX_VERTICES> mask;
            // Começa com a adjacência de v
            mask = adjacency[v];
            // Intersecta com a adjacência de todos os vértices da clique atual
            for (int u : currentClique) {
                // Ignora o caso em que u é o próprio v
                if (u == v) continue;
                // Interseção bit a bit (AND)
                mask &= adjacency[u];
            }

            // Cria novos candidatos
            vector<int> newCandidates;
            // Reserva espaço para o vetor
            newCandidates.reserve(i); 
            // Considera apenas os vértices que aparecem na máscara e estão em 'candidates[0..i-1]'
            for (int j = 0; j < i; j++) {
                int c = candidates[j];
                // Se o vértice c está na máscara (é 1), adiciona-o aos novos candidatos
                if (mask.test(c)) {
                    newCandidates.push_back(c);
                }
            }
            // Chama recursivamente a função
            findMaximumClique(newCandidates, currentClique, bestClique);
        }

        // Backtrack da função recursiva
        currentClique.pop_back();
        // Pula para o próximo vértice
        skip_vertex: ;
    }
}

vector<vector<int>> LerGrafo(const string& nomeArquivo, int &numVertices) {
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
