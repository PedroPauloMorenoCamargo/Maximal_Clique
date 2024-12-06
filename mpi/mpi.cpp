#include <iostream>
#include <vector>
#include <algorithm>
#include <fstream>
#include <chrono>
#include <stdexcept>
#include <mpi.h>
#include <cstdint>

using namespace std;
using namespace std::chrono;

/*  
A matriz de adjacência será representada usando um bitgraph. 
Cada linha do grafo é um vector<uint64_t>, onde cada bit corresponde a uma aresta entre dois vértices.
O número de blocos de 64 bits necessários é dado por (numVertices + 63) / 64.
*/



void findMaximumClique(const vector<vector<uint64_t>> &bitGraph, vector<int> &currentClique, vector<int> &maxClique, const vector<int> &candidates,
                       vector<uint64_t> &intersectionMask,
                       int &globalBest) {
    // Se não houver mais candidatos, a clique atual é a maior encontrada
    if (candidates.empty()) {
        int size = (int)currentClique.size();
        if (size > globalBest) {
            globalBest = size;
            maxClique = currentClique;
        }
        return;
    }

    // Podagem: Se o tamanho da clique atual mais o número de candidatos restantes não puder exceder o tamanho da maior clique global, retorna
    if ((int)currentClique.size() + (int)candidates.size() <= globalBest) {
        return;
    }

    // Itera sobre os candidatos em ordem reversa
    for (int i = (int)candidates.size() - 1; i >= 0; i--) {
        // Vértice atual
        int v = candidates[i];
        // Segunda podagem: Se o tamanho da clique atual mais o número de candidatos restantes não puder exceder o tamanho da maior clique global, retorna
        if ((int)currentClique.size() + i + 1 <= globalBest) {
            return;
        }

        currentClique.push_back(v);

        // Constrói nova máscara de interseção
        vector<uint64_t> newMask = intersectionMask;
        { 
            // Intersecção da máscara com a linha do bitgraph correspondente ao vértice v
            const vector<uint64_t> &vRow = bitGraph[v];
            for (size_t b = 0; b < newMask.size(); b++) {
                newMask[b] &= vRow[b];
            }
        }

        // Filtra os candidatos
        vector<int> newCandidates;
        newCandidates.reserve(candidates.size());
        for (int j = i - 1; j >= 0; j--) {
            int u = candidates[j];
            // Verifica se u está na nova máscara
            if ((newMask[u >> 6] & (1ULL << (u & 63))) != 0ULL) {
                newCandidates.push_back(u);
            }
        }

        // Chama recursivamente a função
        findMaximumClique(bitGraph, currentClique, maxClique, newCandidates, newMask, globalBest);

        // Backtrack
        currentClique.pop_back();
    }
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

// Constrói o bitgraph a partir da matriz de adjacência
static vector<vector<uint64_t>> buildBitGraph(const vector<vector<int>> &graph) {
    // Número de vértices
    int n = (int)graph.size();
    // Número de blocos de 64 bits necessários
    int numBlocks = (n + 63) / 64;
    // Inicializa o bitgraph
    vector<vector<uint64_t>> bitGraph(n, vector<uint64_t>(numBlocks, 0ULL));

    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            // Se houver uma aresta entre os vértices i e j, seta o bit correspondente
            if (graph[i][j]) {
                // Atribui 1 ao bit (j % 64) do bloco (j / 64) da linha i
                bitGraph[i][j >> 6] |= (1ULL << (j & 63));
            }
        }
    }

    return bitGraph;
}

int main(int argc, char* argv[]) {
    MPI_Init(&argc, &argv);

    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    int numVertices = 0;
    vector<vector<int>> graph;

    // Se o rank for 0, lê o grafo do arquivo
    if (rank == 0) {
        if (argc < 2) {
            cerr << "Usage: " << argv[0] << " <input_file>" << endl;
            MPI_Abort(MPI_COMM_WORLD, 1);
        }

        string nomeArquivo = argv[1];
        graph = LerGrafo(nomeArquivo, numVertices);
    }

    // Broadcast do número de vértices
    MPI_Bcast(&numVertices, 1, MPI_INT, 0, MPI_COMM_WORLD);
    // Se o número de vértices for 0, encerra o programa
    if (numVertices == 0) {
        MPI_Finalize();
        return 0;
    }

    // Transfere a matriz de adjacência para um vetor para broadcast
    vector<int> flatGraph;
    if (rank == 0) {
        flatGraph.reserve(numVertices * numVertices);
        for (int i = 0; i < numVertices; ++i) {
            for (int j = 0; j < numVertices; ++j) {
                flatGraph.push_back(graph[i][j]);
            }
        }
    } else {
        flatGraph.resize(numVertices * numVertices);
    }
    // Broadcast da matriz de adjacência
    MPI_Bcast(flatGraph.data(), numVertices * numVertices, MPI_INT, 0, MPI_COMM_WORLD);
    // Reconstrói a matriz de adjacência para os processos diferentes de 0
    if (rank != 0) {
        graph.assign(numVertices, vector<int>(numVertices, 0));
        for (int i = 0; i < numVertices; i++) {
            for (int j = 0; j < numVertices; j++) {
                graph[i][j] = flatGraph[i * numVertices + j];
            }
        }
    }

    // Constrói o bitgraph a partir da matriz de adjacência
    vector<vector<uint64_t>> bitGraph = buildBitGraph(graph);
    int numBlocks = (numVertices + 63) / 64;

    // O mestre (rank 0) ordena os vértices por grau decrescente
    vector<int> candidates;
    if (rank == 0) {
        candidates.reserve(numVertices);
        vector<pair<int, int>> degreeVertexPairs;
        degreeVertexPairs.reserve(numVertices);
        for (int i = 0; i < numVertices; ++i) {
            int degree = 0;
            for (int j = 0; j < numVertices; ++j) {
                degree += graph[i][j];
            }
            degreeVertexPairs.push_back({degree, i});
        }
        sort(degreeVertexPairs.rbegin(), degreeVertexPairs.rend());

        for (auto &p : degreeVertexPairs) {
            candidates.push_back(p.second);
        }
    }

    // Distribui os candidatos entre os processos
    int totalCandidates = 0;
    if (rank == 0) {
        totalCandidates = (int)candidates.size();
    }
    // Broadcast do total de candidatos
    MPI_Bcast(&totalCandidates, 1, MPI_INT, 0, MPI_COMM_WORLD);

    // Divide os candidatos entre os processos
    int baseCount = totalCandidates / size;
    // Resto da divisão
    int remainder = totalCandidates % size;
    // Número de candidatos do processo atual
    int myCount = baseCount + (rank < remainder ? 1 : 0);
    // Deslocamento dos candidatos do processo atual  
    vector<int> sendCounts, displs;
    sendCounts.resize(size, 0);
    displs.resize(size, 0);
    // Calcula os deslocamentos e o número de candidatos de cada processo
    if (rank == 0) {
        for (int r = 0; r < size; r++) {
            sendCounts[r] = baseCount + (r < remainder ? 1 : 0);
        }
        for (int r = 1; r < size; r++) {
            displs[r] = displs[r - 1] + sendCounts[r - 1];
        }
    }
    // Distribui os candidatos
    vector<int> myCandidates(myCount, -1);
    // Envia os candidatos para os processos
    MPI_Scatterv((rank == 0 ? candidates.data() : nullptr),
                 sendCounts.data(),
                 displs.data(),
                 MPI_INT,
                 myCandidates.data(),
                 myCount,
                 MPI_INT,
                 0,
                 MPI_COMM_WORLD);

    MPI_Barrier(MPI_COMM_WORLD);
    auto start_time = high_resolution_clock::now();

    // Cada processo executa a busca localmente
    vector<int> currentClique;
    currentClique.reserve(numVertices);
    vector<int> maxCliqueLocal;
    vector<int> maxCliqueGlobal;

    // Inicializa a maior clique global
    int globalBestLocal = 0;
    // Máscara completa com todos os bits inicialmente habilitados (todos os vértices disponíveis).
    vector<uint64_t> fullMask(numBlocks, ~0ULL);
    // Bits extras no último bloco
    int extraBits = (numBlocks * 64) - numVertices;
    // Se houver bits extras, ajusta a máscara
    if (extraBits > 0) {
        // Máscara para os bits extras
        uint64_t mask = (1ULL << (64 - extraBits)) - 1ULL;
        // Zera os bits extras
        fullMask[numBlocks - 1] &= mask;
    }
    // Para cada candidato do processo
    for (int c : myCandidates) {
        currentClique.clear();
        currentClique.push_back(c);

        // Cria uma máscara de bits para a relação de adjacência de c com todos os vértices da clique atual
        vector<uint64_t> cMask = bitGraph[c];
        // Filtra os candidatos que são adjacentes a c
        vector<int> filteredCandidates;
        filteredCandidates.reserve(totalCandidates);
        for (int u : candidates) {
            if (u != c && (cMask[u >> 6] & (1ULL << (u & 63))) != 0ULL) {
                filteredCandidates.push_back(u);
            }
        }

        // Chama a função recursiva
        findMaximumClique(bitGraph, currentClique, maxCliqueLocal, filteredCandidates, cMask, globalBestLocal);
    }

    int localMaxSize = (int)maxCliqueLocal.size();

    vector<int> allSizes(size, 0);
    // Coleta o tamanho da maior clique de cada processo
    MPI_Gather(&localMaxSize, 1, MPI_INT, allSizes.data(), 1, MPI_INT, 0, MPI_COMM_WORLD);

    int maxSizeGlobal = 0;
    int bestRank = 0;
    if (rank == 0) {
        for (int r = 0; r < size; r++) {
            // Atualiza o tamanho da maior clique global
            if (allSizes[r] > maxSizeGlobal) {
                maxSizeGlobal = allSizes[r];
                bestRank = r;
            }
        }
    }
    // Broadcast do rank do processo com a maior clique
    MPI_Bcast(&bestRank, 1, MPI_INT, 0, MPI_COMM_WORLD);
    // Broadcast do tamanho da maior clique global
    MPI_Bcast(&maxSizeGlobal, 1, MPI_INT, 0, MPI_COMM_WORLD);
    // Se o processo atual tem a maior clique, envia a clique para o processo 0
    if (rank == bestRank) {
        MPI_Send(maxCliqueLocal.data(), maxSizeGlobal, MPI_INT, 0, 0, MPI_COMM_WORLD);
    }
    // O processo 0 recebe a maior clique global
    if (rank == 0) {
        maxCliqueGlobal.resize(maxSizeGlobal);
        MPI_Recv(maxCliqueGlobal.data(), maxSizeGlobal, MPI_INT, bestRank, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }

    MPI_Barrier(MPI_COMM_WORLD);
    auto end_time = high_resolution_clock::now();
    // O processo 0 exibe o resultado
    if (rank == 0) {
        auto duration = duration_cast<std::chrono::duration<double>>(end_time - start_time);
        cout << "Maximum Clique Size: " << maxCliqueGlobal.size() << endl;
        cout << "Vertices in the Maximum Clique: ";
        for (int vertex : maxCliqueGlobal) {
            cout << vertex + 1 << " "; // Convert to 1-based indexing
        }
        cout << endl;
        cout << "Time taken: " << duration.count() << " seconds" << endl;
    }

    MPI_Finalize();
    return 0;
}
