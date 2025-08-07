#include <iostream>
#include <vector>
#include <iomanip>
#include <cmath>
#include <chrono>
#include <algorithm>  // Added for min/max

using namespace std;
using Clock = chrono::high_resolution_clock;

// Parameters (unchanged)
const int max_days = 36500;
const int max_damages = 11;
const double C = 5.0;
const double R = 100.0;
const double gamma = 0.99;
const double tolerance = 1e-6;
const int max_iterations = 100000;

double get_p(int day) {
    double p = 0.01 + (0.19 * static_cast<double>(day)) / max_days;
    return min(max(p, 0.0), 1.0);
}

inline int state_index(int d, int dmg) {
    return d * (max_damages + 1) + dmg;
}

int main() {
    const int n_states = (max_days + 1) * (max_damages + 1);
    vector<double> V(n_states, 0.0);

    auto start_time = Clock::now();

    for (int iter = 0; iter < max_iterations; ++iter) {
        double max_diff = 0.0;

        // Process states in REVERSE order of days
        for (int d = max_days; d >= 0; --d) {
            for (int dmg = 0; dmg <= max_damages; ++dmg) {
                int idx = state_index(d, dmg);
                double old_value = V[idx];
                double new_value;

                if (dmg == max_damages) {
                    new_value = R + gamma * V[state_index(0, 0)];
                }
                else {
                    // Maintenance cost (depends on reset state d=0)
                    double cost_maintenance = C + gamma * V[state_index(0, dmg)];

                    // No maintenance cost (depends on future states d+1)
                    double p_dmg = get_p(d);
                    double cost_no_maintenance;

                    if (dmg + 1 < max_damages) {
                        int next1 = state_index(min(d + 1, max_days), dmg);
                        int next2 = state_index(min(d + 1, max_days), dmg + 1);
                        cost_no_maintenance = gamma * ((1.0 - p_dmg) * V[next1] + p_dmg * V[next2]);
                    }
                    else {
                        int next1 = state_index(min(d + 1, max_days), dmg);
                        cost_no_maintenance = gamma * ((1.0 - p_dmg) * V[next1] + p_dmg * (R + gamma * V[state_index(0, 0)]));
                    }

                    new_value = min(cost_maintenance, cost_no_maintenance);
                }

                double diff = abs(new_value - old_value);
                max_diff = max(max_diff, diff);
                V[idx] = new_value;  // In-place update
            }
        }

        if (max_diff < tolerance) {
            cout << "Converged in " << iter << " iterations.\n";
            break;
        }
    }

    auto end_time = Clock::now();
    chrono::duration<double> elapsed = end_time - start_time;
    cout << fixed << setprecision(4);
    cout << "\nTime elapsed during optimized Gauss-Seidel: " << elapsed.count() << " seconds\n";

    return 0;
}