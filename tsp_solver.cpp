#include <iostream>
#include <vector>
#include <fstream>
#include <sstream>
#include <cmath>
#include <algorithm>
#include <climits>
#include <chrono>
#include <queue>
#include <unordered_map>
#include <random>
#include <limits>

using namespace std;
using namespace std::chrono;

struct Point {
    double x, y;
    Point(double x = 0, double y = 0) : x(x), y(y) {}
};

class TSPSolver {
private:
    vector<Point> cities;
    vector<vector<double>> dist;
    int n;

public:
    TSPSolver() : n(0) {}
    

    bool loadTSP(const string& filename) {
        ifstream file(filename);
        if (!file.is_open()) {
            cout << "Error: Cannot open file " << filename << endl;
            return false;
        }
        
        string line;
        bool readingCoords = false;
        cities.clear();
        
        while (getline(file, line)) {
            if (line.find("NODE_COORD_SECTION") != string::npos) {
                readingCoords = true;
                continue;
            }
            if (line.find("EOF") != string::npos) break;
            
            if (readingCoords && !line.empty()) {
                istringstream iss(line);
                int id;
                double x, y;
                if (iss >> id >> x >> y) {
                    cities.push_back(Point(x, y));
                }
            }
        }
        
        n = cities.size();
        cout << "Loaded " << n << " cities from " << filename << endl;
        

        buildDistanceMatrix();
        return true;
    }
    
    void buildDistanceMatrix() {
        dist.assign(n, vector<double>(n, 0));
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                if (i != j) {
                    double dx = cities[i].x - cities[j].x;
                    double dy = cities[i].y - cities[j].y;
                    dist[i][j] = sqrt(dx * dx + dy * dy);
                }
            }
        }
    }
    
    double calculateTourCost(const vector<int>& tour) {
        double cost = 0;
        for (int i = 0; i < tour.size(); i++) {
            int next = (i + 1) % tour.size();
            cost += dist[tour[i]][tour[next]];
        }
        return cost;
    }
    

    vector<int> mstApproximation() {
        cout << "Running MST-based 2-approximation..." << endl;
        

        vector<bool> inMST(n, false);
        vector<vector<int>> mst(n);
        priority_queue<pair<double, pair<int, int>>, 
                      vector<pair<double, pair<int, int>>>, 
                      greater<pair<double, pair<int, int>>>> pq;
        
        inMST[0] = true;
        for (int i = 1; i < n; i++) {
            pq.push({dist[0][i], {0, i}});
        }
        
        while (!pq.empty()) {
            auto edge = pq.top();
            pq.pop();
            
            double weight = edge.first;
            int u = edge.second.first;
            int v = edge.second.second;
            
            if (inMST[v]) continue;
            
            inMST[v] = true;
            mst[u].push_back(v);
            mst[v].push_back(u);
            
            for (int i = 0; i < n; i++) {
                if (!inMST[i]) {
                    pq.push({dist[v][i], {v, i}});
                }
            }
        }
        

        vector<int> tour;
        vector<bool> visited(n, false);
        
        function<void(int)> dfs = [&](int u) {
            visited[u] = true;
            tour.push_back(u);
            for (int v : mst[u]) {
                if (!visited[v]) {
                    dfs(v);
                }
            }
        };
        
        dfs(0);
        return tour;
    }
    
    pair<vector<int>, double> heldKarp() {
        cout << "Running Held-Karp algorithm..." << endl;
        
        if (n > 20) {
            cout << "Warning: Held-Karp may be too slow for " << n << " cities" << endl;
            return {{}, -1};
        }
        
        
        int maxMask = 1 << (n - 1); 
        vector<vector<double>> dp(maxMask, vector<double>(n, 1e9));
        vector<vector<int>> parent(maxMask, vector<int>(n, -1));
        
    
        for (int i = 1; i < n; i++) {
            int mask = 1 << (i - 1); 
            dp[mask][i] = dist[0][i];
            parent[mask][i] = 0;
        }
        
        // Fill DP table
        for (int mask = 1; mask < maxMask; mask++) {
            int bitCount = __builtin_popcount(mask);
            if (bitCount < 2) continue;
            
            for (int u = 1; u < n; u++) {
                int uBit = 1 << (u - 1);
                if (!(mask & uBit)) continue;
                
                int prevMask = mask ^ uBit; 
                
                for (int v = 1; v < n; v++) {
                    int vBit = 1 << (v - 1);
                    if (!(prevMask & vBit)) continue; 
                    if (dp[prevMask][v] >= 1e9) continue; 
                    
                    double cost = dp[prevMask][v] + dist[v][u];
                    if (cost < dp[mask][u]) {
                        dp[mask][u] = cost;
                        parent[mask][u] = v;
                    }
                }
            }
        }
        
    
        int finalMask = maxMask - 1;
        double minTourCost = 1e9;
        int lastCity = -1;
        
        for (int i = 1; i < n; i++) {
            if (dp[finalMask][i] < 1e9) {
                double cost = dp[finalMask][i] + dist[i][0];
                if (cost < minTourCost) {
                    minTourCost = cost;
                    lastCity = i;
                }
            }
        }
        
        if (lastCity == -1 || minTourCost >= 1e9) {
            cout << "Error: No valid tour found in Held-Karp" << endl;
            cout << "Debug: finalMask = " << finalMask << ", checking costs:" << endl;
            for (int i = 1; i < n; i++) {
                cout << "  dp[" << finalMask << "][" << i << "] = " << dp[finalMask][i] << endl;
            }
            return {{}, -1};
        }
        

        vector<int> tour;
        int mask = finalMask;
        int curr = lastCity;
        
        while (curr != 0 && mask != 0) {
            tour.push_back(curr);
            int prev = parent[mask][curr];
            if (prev == -1) break;
            
            if (prev != 0) {
                mask ^= (1 << (curr - 1)); 
            } else {
                break; 
            }
            curr = prev;
        }
        
        tour.push_back(0); 
        reverse(tour.begin(), tour.end());
        
    
        if (tour.size() != n) {
            cout << "Warning: Reconstructed tour has " << tour.size() << " cities instead of " << n << endl;
            cout << "Tour: ";
            for (int city : tour) cout << city << " ";
            cout << endl;
        }
        
        return {tour, minTourCost};
    }


    vector<int> sqrtLookaheadAlgo() {
        cout << "Running Sqrt-N Lookahead algorithm..." << endl;
        
        int lookaheadSize = max(1, (int)sqrt(n));
        cout << "  Lookahead size: " << lookaheadSize << " cities" << endl;
        
        vector<int> tour;
        vector<bool> visited(n, false);
        

        int currentCity = 0;
        tour.push_back(currentCity);
        visited[currentCity] = true;
        

        for (int step = 1; step < n; step++) {
            vector<pair<double, int>> candidates;
            
     
            for (int i = 0; i < n; i++) {
                if (!visited[i]) {
                    candidates.push_back({dist[currentCity][i], i});
                }
            }
            

            sort(candidates.begin(), candidates.end());
            int actualLookahead = min(lookaheadSize, (int)candidates.size());
            
           
            double bestScore = numeric_limits<double>::max();
            int bestCity = -1;
            
            for (int i = 0; i < actualLookahead; i++) {
                int candidateCity = candidates[i].second;
                double immediateCost = candidates[i].first;
                
        
                double futureCost = calculateFutureCostHeuristic(candidateCity, visited);
                
      
                double totalScore = immediateCost + 0.3 * futureCost; 
                
                if (totalScore < bestScore) {
                    bestScore = totalScore;
                    bestCity = candidateCity;
                }
            }
            

            if (bestCity != -1) {
                tour.push_back(bestCity);
                visited[bestCity] = true;
                currentCity = bestCity;
            }
        }
        
        return tour;
    }
    

    double calculateFutureCostHeuristic(int fromCity, const vector<bool>& visited) {

        double totalDist = 0;
        int unvisitedCount = 0;
        
        for (int i = 0; i < n; i++) {
            if (!visited[i] && i != fromCity) {
                totalDist += dist[fromCity][i];
                unvisitedCount++;
            }
        }
        
        if (unvisitedCount == 0) {

            return dist[fromCity][0];
        }
        

        return totalDist / unvisitedCount;
    }

    void runExperiments() {
        cout << "\n=== TSP Experiments ===" << endl;
        
        auto start = high_resolution_clock::now();
        vector<int> mstTour = mstApproximation();
        auto end = high_resolution_clock::now();
        double mstTime = duration_cast<microseconds>(end - start).count() / 1000.0;
        double mstCost = calculateTourCost(mstTour);
        
        cout << "\nMST Approximation:" << endl;
        cout << "  Cost: " << mstCost << endl;
        cout << "  Time: " << mstTime << " ms" << endl;
        
        double hkCost = numeric_limits<double>::max();
        double hkTime = 0;
        vector<int> hkTour;
        
        if (n <= 15) {
            start = high_resolution_clock::now();
            auto hkResult = heldKarp();
            end = high_resolution_clock::now();
            hkTime = duration_cast<microseconds>(end - start).count() / 1000.0;
            hkTour = hkResult.first;
            hkCost = hkResult.second;
            
            cout << "\nHeld-Karp:" << endl;
            if (hkCost != -1) {
                cout << "  Cost: " << hkCost << endl;
                cout << "  Time: " << hkTime << " ms" << endl;

                double verifiedCost = calculateTourCost(hkTour);
                cout << "  Verified cost: " << verifiedCost << endl;
                hkCost = verifiedCost;
            } else {
                cout << "  Failed to find solution" << endl;
                hkCost = numeric_limits<double>::max();
            }
        } else {
            cout << "\nHeld-Karp: Skipped (too large for n=" << n << ")" << endl;
        }
        
        start = high_resolution_clock::now();
        vector<int> lookaheadTour = sqrtLookaheadAlgo();
        end = high_resolution_clock::now();
        double lookaheadTime = duration_cast<microseconds>(end - start).count() / 1000.0;
        double lookaheadCost = calculateTourCost(lookaheadTour);

        cout << "\nSqrt-N Lookahead Algorithm:" << endl;
        cout << "  Cost: " << lookaheadCost << endl;
        cout << "  Time: " << lookaheadTime << " ms" << endl;
        

        cout << "\n=== Summary ===" << endl;
        cout << "Dataset size: " << n << " cities" << endl;
        
        double bestCost;
        if (hkCost != numeric_limits<double>::max()) {
            bestCost = min({mstCost, lookaheadCost, hkCost});
        } else {
            bestCost = min(mstCost, lookaheadCost);
        }
        
        cout << "MST Approximation ratio: " << (mstCost / bestCost) << endl;
        cout << "Sqrt-N Lookahead ratio: " << (lookaheadCost / bestCost) << endl;
        
        if (hkCost != numeric_limits<double>::max()) {
            cout << "Held-Karp ratio: " << (hkCost / bestCost) << " (optimal)" << endl;
        }
        
        cout << "Best solution cost: " << bestCost << endl;
        

        cout << "\n=== Performance Analysis ===" << endl;
        cout << "Fastest algorithm: ";
        if (hkTime > 0 && lookaheadTime <= mstTime && lookaheadTime <= hkTime) {
            cout << "Sqrt-N Lookahead (" << lookaheadTime << " ms)" << endl;
        } else if (hkTime > 0 && mstTime <= lookaheadTime && mstTime <= hkTime) {
            cout << "MST (" << mstTime << " ms)" << endl;
        } else if (hkTime > 0) {
            cout << "Held-Karp (" << hkTime << " ms)" << endl;
        } else if (lookaheadTime <= mstTime) {
            cout << "Sqrt-N Lookahead (" << lookaheadTime << " ms)" << endl;
        } else {
            cout << "MST (" << mstTime << " ms)" << endl;
        }
        
        cout << "Best quality: ";
        if (hkCost != numeric_limits<double>::max() && hkCost <= lookaheadCost && hkCost <= mstCost) {
            cout << "Held-Karp (optimal)" << endl;
        } else if (lookaheadCost <= mstCost) {
            cout << "Sqrt-N Lookahead Algorithm" << endl;
        } else {
            cout << "MST Approximation" << endl;
        }
    }
};

int main(int argc, char *argv[]) {
    TSPSolver solver;
    vector<string> testFiles;

    if (argc <= 1) {
        // Default test files if no arguments provided
        testFiles = {
	    "datasets/test.tsp",
	    "datasets/xql662.tsp",
	    "datasets/a280.tsp",
	    "datasets/kz9976.tsp",
	    "datasets/mona-lisa100K.tsp"
	};
    } else {
        // Convert command line arguments to strings
        testFiles.assign(argv + 1, argv + argc);
    }

    for (const string& filename : testFiles) {
        cout << "\n" << string(50, '=') << endl;
        cout << "Testing with: " << filename << endl;
        cout << string(50, '=') << endl;

        if (solver.loadTSP(filename)) {
            solver.runExperiments();
        } else {
            cout << "Skipping " << filename << " (file not found)" << endl;
        }
    }

    return 0;
}

