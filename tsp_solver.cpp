#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <queue>
#include <unordered_map>
#include <vector>
#include <cmath>
#include <functional>
#include <limits>

using namespace std;
using namespace std::chrono;

struct Point {
	double x, y;
	explicit Point(double x = 0, double y = 0) : x(x), y(y) {}
};

class TSPSolver {
private:
	int num_cities;
	vector<Point> cities;
	vector<vector<double>> dist;

public:
	TSPSolver() : num_cities(0) {}
	
	auto loadTSP(const string &filename) -> bool {
		ifstream file(filename);
		if (!file.is_open()) {
			cout << "Error: Cannot open file " << filename << endl;
			return false;
		}
		
		cities.clear();
		string line;
		
		while (std::getline(file, line) && !line.contains("NODE_COORD_SECTION"));
		while (std::getline(file, line) && !line.contains("EOF")) {
			istringstream iss(line);
			int _, x, y;
			iss >> _ >> x >> y;
			cities.emplace_back(x, y);
		}
		
		num_cities = static_cast<int>(cities.size());
		cout << "Loaded " << num_cities << " cities from " << filename << endl;
		
		buildDistanceMatrix();
		return true;
	}
	
	void buildDistanceMatrix() {
		dist.assign(num_cities, vector<double>(num_cities, 0.f));
		for (int i = 0; i < num_cities; i++) {
			for (int j = 0; j < num_cities; j++) {
				if (i != j) {
					double dx = cities[i].x - cities[j].x;
					double dy = cities[i].y - cities[j].y;
					dist[i][j] = sqrt(dx * dx + dy * dy);
				}
			}
		}
	}
	
	auto calculateTourCost(const vector<int> &tour) -> double {
		double cost = 0;
		for (int i = 0; i < tour.size(); i++) {
			int next = (i + 1) % tour.size();
			cost += dist[tour[i]][tour[next]];
		}
		return cost;
	}
	
	auto mstApproximation() -> vector<int> {
		cout << "Running MST-based 2-approximation..." << endl;
		
		vector<bool> inMST(num_cities, false);
		vector<vector<int>> mst(num_cities);
		priority_queue<pair<double, pair<int, int>>, vector<pair<double, pair<int, int>>>, greater<>> pq;
		
		inMST[0] = true;
		for (int i = 1; i < num_cities; i++) {
			pq.push({dist[0][i], {0, i}});
		}
		
		while (!pq.empty()) {
			auto edge = pq.top();
			pq.pop();
			
			int u = edge.second.first;
			int v = edge.second.second;
			
			if (inMST[v]) continue;
			
			inMST[v] = true;
			mst[u].push_back(v);
			mst[v].push_back(u);
			
			for (int i = 0; i < num_cities; i++) {
				if (!inMST[i]) pq.push({dist[v][i], {v, i}});
			}
		}
		
		vector<int> tour;
		vector<bool> visited(num_cities, false);
		
		function<void(int)> dfs = [&](int u) {
			visited[u] = true;
			tour.emplace_back(u);
			for (int v: mst[u]) {
				if (!visited[v]) dfs(v);
			}
		};
		
		dfs(0);
		return tour;
	}
	
	auto heldKarp() -> pair<vector<int>, double> {
		cout << "Running Held-Karp algorithm..." << endl;
		
		if (num_cities > 20) {
			cout << "Warning: Held-Karp may be too slow for " << num_cities << " cities" << endl;
			return {{}, -1};
		}
		
		int maxMask = 1 << (num_cities - 1);
		vector<vector<double>> dp(maxMask, vector<double>(num_cities, 1e9));
		vector<vector<int>> parent(maxMask, vector<int>(num_cities, -1));
		
		for (int i = 1; i < num_cities; i++) {
			int mask = 1 << (i - 1);
			dp[mask][i] = dist[0][i];
			parent[mask][i] = 0;
		}
		
		// Fill DP table
		for (int mask = 1; mask < maxMask; mask++) {
			int bitCount = __builtin_popcount(mask);
			if (bitCount < 2) continue;
			
			for (int u = 1; u < num_cities; u++) {
				int uBit = 1 << (u - 1);
				if (!(mask & uBit)) continue;
				
				int prevMask = mask ^ uBit;
				
				for (int v = 1; v < num_cities; v++) {
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
		
		for (int i = 1; i < num_cities; i++) {
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
			for (int i = 1; i < num_cities; i++) {
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
		std::ranges::reverse(tour);
		
		if (tour.size() != num_cities) {
			cout << "Warning: Reconstructed tour has " << tour.size() << " cities instead of " << num_cities << endl;
			cout << "Tour: ";
			for (int city: tour) cout << city << " ";
			cout << endl;
		}
		
		return {tour, minTourCost};
	}
	
	auto sqrtLookaheadAlgo() -> vector<int> {
		cout << "Running Sqrt-N Lookahead algorithm..." << endl;
		
		int lookaheadSize = max(1, (int) sqrt(num_cities));
		cout << "  Lookahead size: " << lookaheadSize << " cities" << endl;
		
		vector<int> tour;
		vector<bool> visited(num_cities, false);
		
		int currentCity = 0;
		tour.push_back(currentCity);
		visited[currentCity] = true;
		
		for (int step = 1; step < num_cities; step++) {
			vector<pair<double, int>> candidates;
			
			for (int i = 0; i < num_cities; i++) {
				if (!visited[i]) {
					candidates.emplace_back(dist[currentCity][i], i);
				}
			}
			
			sort(candidates.begin(), candidates.end());
			int actualLookahead = min(lookaheadSize, (int) candidates.size());
			
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
	
	auto calculateFutureCostHeuristic(int fromCity, const vector<bool> &visited) -> double {
		
		double totalDist = 0;
		int unvisitedCount = 0;
		
		for (int i = 0; i < num_cities; i++) {
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
		
		if (num_cities <= 15) {
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
			cout << "\nHeld-Karp: Skipped (too large for n=" << num_cities << ")" << endl;
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
		cout << "Dataset size: " << num_cities << " cities" << endl;
		
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
		double minTime = std::numeric_limits<double>::max();
		std::string fastest;
		
		if (hkTime > 0 && hkTime < minTime) minTime = hkTime, fastest = "Held-Karp";
		if (lookaheadTime < minTime) minTime = lookaheadTime, fastest = "Sqrt-N Lookahead";
		if (mstTime < minTime) minTime = mstTime, fastest = "MST";
		std::cout << "Fastest algorithm: " << fastest << " (" << minTime << " ms)" << std::endl;
		
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

auto main(int argc, char *argv[]) -> int {
	TSPSolver solver;
	vector<string> testFiles;
	
	if (argc <= 1) {
		testFiles = {
				"datasets/test.tsp",
				"datasets/xql662.tsp",
				"datasets/a280.tsp",
				"datasets/kz9976.tsp",
				"datasets/mona-lisa100K.tsp"
		};
	} else {
		testFiles.assign(argv + 1, argv + argc);
	}
	
	for (const string &filename: testFiles) {
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

