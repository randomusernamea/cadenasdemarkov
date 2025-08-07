#include <iostream>
#include <vector>
#include <iomanip>
#include <cmath>
#include <chrono>

using namespace std;
using Clock = chrono::high_resolution_clock;

// Parameters
const int max_days = 36500;
const int max_damages = 23;
const double C = 5.0;            // Maintenance cost
const double R = 100.0;          // Repair cost
const double gamma = 0.99;       // Discount factor
const double tolerance = 1e-6;
const int max_iterations = 100000;

// Get damage probability for a given day since maintenance
double get_p(int day) {
    double p = 0.01 + (0.19 * day) / max_days; // Linear from 0.01 to 0.2
    return min(max(p, 0.0), 1.0);
}

// Convert (days, dmg) into a 1D index
inline int state_index(int d, int dmg) {
    return d * (max_damages + 1) + dmg;
}

int main() {
    const int n_states = (max_days + 1) * (max_damages + 1);
    vector<double> V(n_states, 0.0);
    vector<double> V_new(n_states, 0.0);

    auto start_time = Clock::now();

    for (int iter = 0; iter < max_iterations; ++iter) {
        double max_diff = 0.0;

        for (int d = 0; d <= max_days; ++d) {
            for (int dmg = 0; dmg <= max_damages; ++dmg) {
                int idx = state_index(d, dmg);

                if (dmg == max_damages) {
                    // Forced repair
                    V_new[idx] = R + gamma * V[state_index(0, 0)];
                    continue;
                }

                // Action 1: Do maintenance
                double cost_maintenance = C + gamma * V[state_index(0, dmg)];

                // Action 2: Do nothing
                double p_dmg = get_p(d);
                double cost_no_maintenance;

                if (dmg + 1 < max_damages) {
                    int next1 = state_index(min(d + 1, max_days), dmg);
                    int next2 = state_index(min(d + 1, max_days), dmg + 1);
                    double expected = (1.0 - p_dmg) * V[next1] + p_dmg * V[next2];
                    cost_no_maintenance = gamma * expected;
                }
                else {
                    int next1 = state_index(min(d + 1, max_days), dmg);
                    double expected = (1.0 - p_dmg) * V[next1] + p_dmg * (R + gamma * V[state_index(0, 0)]);
                    cost_no_maintenance = gamma * expected;
                }

                V_new[idx] = min(cost_maintenance, cost_no_maintenance);
                max_diff = max(max_diff, abs(V_new[idx] - V[idx]));
            }
        }

        if (max_diff < tolerance) {
            cout << "Converged in " << iter << " iterations.\n";
            break;
        }

        V = V_new;
    }

    auto end_time = Clock::now();
    chrono::duration<double> elapsed = end_time - start_time;
    cout << fixed << setprecision(4);
    cout << "\nTime elapsed during Jacobi method: " << elapsed.count() << " seconds\n";



    return 0;
}
