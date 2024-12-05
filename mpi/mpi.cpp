#include <iostream>
#include <vector>
#include <algorithm>
#include <fstream>
#include <chrono>
#include <stdexcept>
#include <mpi.h>
#include <cstdint>

// Uncomment to enable debugging messages
//#define DEBUG_OUTPUT

using namespace std;
using namespace std::chrono;

// We'll represent adjacency using bitsets. 
// Each row of the graph is a vector<uint64_t> where each bit corresponds to an edge.
// For numVertices, we need (numVertices+63)/64 64-bit blocks.

static inline bool isAdjacent(const vector<uint64_t> &row, int v) {
    return (row[v >> 6] & (1ULL << (v & 63))) != 0ULL;
}

// Intersect candidate sets by checking adjacency to a particular vertex v
// Using a bitset intersection can be done if we keep track of candidate sets as bitsets.
static inline void filterCandidates(const vector<vector<uint64_t>> &bitGraph,
                                    const vector<uint64_t> &currentMask,
                                    const vector<int> &candidates,
                                    vector<int> &newCandidates,
                                    int numBlocks) {
    newCandidates.clear();
    newCandidates.reserve(candidates.size());
    
    // For each candidate, check if it's adjacent to all in currentClique.
    // Instead of checking individually, we know currentMask is a mask of allowed vertices.
    // If we keep a mask of allowed vertices, we can just check if candidate is in that mask.
    // However, since we prune candidates directly, we can just check adjacency to all 
    // in the current clique. But to optimize further, we store a mask representing
    // intersection of adjacency of all vertices in currentClique.
    // For simplicity, we'll do the old-fashioned check here, but using bitsets.
    
    // currentMask is a combined mask of valid vertices. candidates are chosen from those vertices.
    for (int u : candidates) {
        // Check if u is in currentMask
        if ((currentMask[u >> 6] & (1ULL << (u & 63))) != 0ULL) {
            newCandidates.push_back(u);
        }
    }
}

// Builds a mask representing the intersection of adjacency among all vertices in currentClique.
// This mask will be used to quickly filter candidates.
static inline void buildIntersectionMask(const vector<vector<uint64_t>> &bitGraph,
                                         const vector<int> &currentClique,
                                         vector<uint64_t> &mask) {
    if (currentClique.empty()) {
        return;
    }
    // Start with the adjacency of the first vertex
    mask = bitGraph[currentClique[0]];
    // Intersect with others
    for (size_t i = 1; i < currentClique.size(); i++) {
        const vector<uint64_t> &row = bitGraph[currentClique[i]];
        for (size_t b = 0; b < mask.size(); b++) {
            mask[b] &= row[b];
        }
    }
}

// Recursive maximum clique search
// We'll pass a global reference to the best clique size to allow more aggressive pruning.
void findMaximumClique(const vector<vector<uint64_t>> &bitGraph,
                       vector<int> &currentClique,
                       vector<int> &maxClique,
                       const vector<int> &candidates,
                       vector<uint64_t> &intersectionMask,
                       int &globalBest) {
    if (candidates.empty()) {
        int size = (int)currentClique.size();
        if (size > globalBest) {
            globalBest = size;
            maxClique = currentClique;
        }
        return;
    }

    // Pruning: If even taking all candidates can't beat globalBest, stop.
    if ((int)currentClique.size() + (int)candidates.size() <= globalBest) {
        return;
    }

    // We'll pick candidates in reverse order to mimic original behavior
    for (int i = (int)candidates.size() - 1; i >= 0; i--) {
        int v = candidates[i];
        // Prune check again (as we remove candidates)
        if ((int)currentClique.size() + i + 1 <= globalBest) {
            return;
        }

        // Check adjacency to all in currentClique is trivial now because we used intersectionMask
        // We know candidates are already filtered. No need to check adjacency again.
        currentClique.push_back(v);

        // Build intersection mask for new clique
        vector<uint64_t> newMask = intersectionMask;
        { 
            // Intersect with adjacency of v
            const vector<uint64_t> &vRow = bitGraph[v];
            for (size_t b = 0; b < newMask.size(); b++) {
                newMask[b] &= vRow[b];
            }
        }

        // Filter candidates with newMask
        vector<int> newCandidates;
        newCandidates.reserve(candidates.size());
        for (int j = i - 1; j >= 0; j--) {
            int u = candidates[j];
            if ((newMask[u >> 6] & (1ULL << (u & 63))) != 0ULL) {
                newCandidates.push_back(u);
            }
        }

        // Recursive call
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

// Convert adjacency matrix to a bitset representation
// Each row: vector<uint64_t>, where the j-th bit indicates if there's an edge.
static vector<vector<uint64_t>> buildBitGraph(const vector<vector<int>> &graph) {
    int n = (int)graph.size();
    int numBlocks = (n + 63) / 64;

    vector<vector<uint64_t>> bitGraph(n, vector<uint64_t>(numBlocks, 0ULL));

    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            if (graph[i][j]) {
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

    // Only the master reads the input file
    if (rank == 0) {
        if (argc < 2) {
            cerr << "Usage: " << argv[0] << " <input_file>" << endl;
            MPI_Abort(MPI_COMM_WORLD, 1);
        }

        string nomeArquivo = argv[1];
        graph = LerGrafo(nomeArquivo, numVertices);
    }

    // Broadcast numVertices to all processes
    MPI_Bcast(&numVertices, 1, MPI_INT, 0, MPI_COMM_WORLD);

    if (numVertices == 0) {
        MPI_Finalize();
        return 0;
    }

    // Flatten the graph and broadcast
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

    MPI_Bcast(flatGraph.data(), numVertices * numVertices, MPI_INT, 0, MPI_COMM_WORLD);

    if (rank != 0) {
        graph.assign(numVertices, vector<int>(numVertices, 0));
        for (int i = 0; i < numVertices; i++) {
            for (int j = 0; j < numVertices; j++) {
                graph[i][j] = flatGraph[i * numVertices + j];
            }
        }
    }

    // Build bitGraph
    vector<vector<uint64_t>> bitGraph = buildBitGraph(graph);
    int numBlocks = (numVertices + 63) / 64;

    // The master computes the degrees and prepares candidates
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

    // Distribute candidates among processes
    int totalCandidates = 0;
    if (rank == 0) {
        totalCandidates = (int)candidates.size();
    }
    MPI_Bcast(&totalCandidates, 1, MPI_INT, 0, MPI_COMM_WORLD);

    int baseCount = totalCandidates / size;
    int remainder = totalCandidates % size;
    int myCount = baseCount + (rank < remainder ? 1 : 0);

    vector<int> sendCounts, displs;
    sendCounts.resize(size, 0);
    displs.resize(size, 0);
    if (rank == 0) {
        for (int r = 0; r < size; r++) {
            sendCounts[r] = baseCount + (r < remainder ? 1 : 0);
        }
        for (int r = 1; r < size; r++) {
            displs[r] = displs[r - 1] + sendCounts[r - 1];
        }
    }

    vector<int> myCandidates(myCount, -1);
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

    // Each process searches for max clique starting from assigned candidates
    vector<int> currentClique;
    currentClique.reserve(numVertices);
    vector<int> maxCliqueLocal;
    vector<int> maxCliqueGlobal;

    // Global best clique size for pruning
    int globalBestLocal = 0;

    vector<uint64_t> fullMask(numBlocks, ~0ULL);
    // Adjust fullMask because we have exactly numVertices bits
    // Clear unused bits in the last block if necessary
    int extraBits = (numBlocks * 64) - numVertices;
    if (extraBits > 0) {
        uint64_t mask = (1ULL << (64 - extraBits)) - 1ULL;
        fullMask[numBlocks - 1] &= mask;
    }

    for (int c : myCandidates) {
        currentClique.clear();
        currentClique.push_back(c);

        // Build mask for candidates adjacent to c
        vector<uint64_t> cMask = bitGraph[c];
        // Filter the global candidates based on cMask
        vector<int> filteredCandidates;
        filteredCandidates.reserve(totalCandidates);
        for (int u : candidates) {
            if (u != c && (cMask[u >> 6] & (1ULL << (u & 63))) != 0ULL) {
                filteredCandidates.push_back(u);
            }
        }

        // Run recursive search
        findMaximumClique(bitGraph, currentClique, maxCliqueLocal, filteredCandidates, cMask, globalBestLocal);
    }

    int localMaxSize = (int)maxCliqueLocal.size();

    vector<int> allSizes(size, 0);
    MPI_Gather(&localMaxSize, 1, MPI_INT, allSizes.data(), 1, MPI_INT, 0, MPI_COMM_WORLD);

    int maxSizeGlobal = 0;
    int bestRank = 0;
    if (rank == 0) {
        for (int r = 0; r < size; r++) {
            if (allSizes[r] > maxSizeGlobal) {
                maxSizeGlobal = allSizes[r];
                bestRank = r;
            }
        }
    }

    MPI_Bcast(&bestRank, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&maxSizeGlobal, 1, MPI_INT, 0, MPI_COMM_WORLD);

    if (rank == bestRank) {
        MPI_Send(maxCliqueLocal.data(), maxSizeGlobal, MPI_INT, 0, 0, MPI_COMM_WORLD);
    }

    if (rank == 0) {
        maxCliqueGlobal.resize(maxSizeGlobal);
        MPI_Recv(maxCliqueGlobal.data(), maxSizeGlobal, MPI_INT, bestRank, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }

    MPI_Barrier(MPI_COMM_WORLD);
    auto end_time = high_resolution_clock::now();

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
