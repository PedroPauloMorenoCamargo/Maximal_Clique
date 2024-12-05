#include <iostream>
#include <vector>
#include <algorithm>
#include <fstream>
#include <chrono>
#include <stdexcept>
#include <bitset>

using namespace std;
using namespace std::chrono;

// Adjust this if you know the maximum number of vertices
// For example if you know maximum vertices <= 2000, use 2048 or the nearest multiple of 64.
// For demonstration, let's pick a large fixed size:
const int MAX_VERTICES = 2048; 
static_assert(MAX_VERTICES % 64 == 0, "MAX_VERTICES should be multiple of 64");

// We'll represent the graph as a vector of bitsets (or an array of bitsets).
// Each vertex has a bitset representing edges to other vertices.
static bitset<MAX_VERTICES> adjacency[MAX_VERTICES];

int bestCliqueSizeGlobal = 0; // Keep track of best clique size globally for pruning

// Function prototypes
void findMaximumClique(const vector<int> &candidates,
                       vector<int> &currentClique,
                       vector<int> &bestClique);

vector<vector<int>> LerGrafo(const string& nomeArquivo, int &numVertices);

int main(int argc, char* argv[]) {
    if (argc < 2) {
        cerr << "Uso: " << argv[0] << " <nome_do_arquivo>" << endl;
        return 1;
    }

    string nomeArquivo = argv[1];
    int numVertices;

    // Lê o grafo
    vector<vector<int>> graph = LerGrafo(nomeArquivo, numVertices);
    if (numVertices > MAX_VERTICES) {
        cerr << "Number of vertices exceeds MAX_VERTICES. Increase MAX_VERTICES." << endl;
        return 1;
    }

    // Build bitset graph representation
    for (int i = 0; i < numVertices; i++) {
        adjacency[i].reset();
        for (int j = 0; j < numVertices; j++) {
            if (graph[i][j] == 1) {
                adjacency[i].set(j);
            }
        }
    }

    // Calcula o grau e ordena os vértices por grau decrescente
    vector<pair<int, int>> degreeVertexPairs;
    degreeVertexPairs.reserve(numVertices);
    for (int i = 0; i < numVertices; ++i) {
        int degree = (int)adjacency[i].count();
        degreeVertexPairs.push_back({degree, i});
    }

    sort(degreeVertexPairs.rbegin(), degreeVertexPairs.rend());

    // Prepara a lista de candidatos
    vector<int> candidates;
    candidates.reserve(numVertices);
    for (auto &p : degreeVertexPairs) {
        candidates.push_back(p.second);
    }

    vector<int> currentClique;
    currentClique.reserve(numVertices);
    vector<int> bestClique;
    bestClique.reserve(numVertices);

    // Tempo
    auto start_time = high_resolution_clock::now();
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

void findMaximumClique(const vector<int> &candidates,
                       vector<int> &currentClique,
                       vector<int> &bestClique) {
    // If candidates is empty, check if current clique is better
    if (candidates.empty()) {
        if ((int)currentClique.size() > bestCliqueSizeGlobal) {
            bestCliqueSizeGlobal = (int)currentClique.size();
            bestClique = currentClique;
        }
        return;
    }

    // Pruning: If currentClique.size() + candidates.size() <= bestCliqueSizeGlobal, no need to continue
    if ((int)currentClique.size() + (int)candidates.size() <= bestCliqueSizeGlobal) {
        return;
    }

    // We'll pick candidates in reverse order (like original code)
    for (int i = (int)candidates.size() - 1; i >= 0; --i) {
        int v = candidates[i];

        // More pruning: if even including v and all candidates before i cannot exceed bestCliqueSizeGlobal, stop
        if ((int)currentClique.size() + i + 1 <= bestCliqueSizeGlobal) {
            return;
        }

        // Check adjacency with current clique:
        // Since we only pick candidates that are already filtered, normally we would do:
        // But we can confirm adjacency quickly:
        // If needed, but since we always filter candidates at each step, all should be adjacent.
        // For safety, we could check:
        for (int u : currentClique) {
            if (!adjacency[u].test(v)) {
                goto skip_vertex;
            }
        }

        // Add v to currentClique
        currentClique.push_back(v);

        // Now filter new candidates:
        // We intersect adjacency masks of the currentClique vertices.
        // Instead of intersecting all, we maintain an intersection bitset:
        {
            bitset<MAX_VERTICES> mask;
            // Start with adjacency of v
            mask = adjacency[v];
            // Intersect with previous clique members
            for (int u : currentClique) {
                if (u == v) continue;
                mask &= adjacency[u];
            }

            // Build newCandidates from mask and from previously available candidates
            vector<int> newCandidates;
            newCandidates.reserve(i); 
            // Only consider vertices that appear in the mask and are in 'candidates[0..i-1]'
            for (int j = 0; j < i; j++) {
                int c = candidates[j];
                if (mask.test(c)) {
                    newCandidates.push_back(c);
                }
            }

            findMaximumClique(newCandidates, currentClique, bestClique);
        }

        // Backtrack
        currentClique.pop_back();

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
