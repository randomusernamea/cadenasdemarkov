#include <iostream>
#include <vector>
#include <iomanip>
#include <cmath>
#include <chrono>
#include <algorithm>

using namespace std;
using Clock = chrono::high_resolution_clock;

// Parameters
const int max_days = 365;
const int max_damages = 10;  // Damage levels: 0 to 10; 11 = failure
const double C = 5.0;
const double R = 100.0;
const double gamma = 0.8;
const double tolerance = 1e-6;
const int max_iterations = 10000;
const double omega = 1;  // SOR parameter

double get_p(int day) {
    double p = 0.01 + (0.19 * static_cast<double>(day)) / max_days;
    return min(max(p, 0.0), 1.0);
}

int main() {
    // Create a 2D value function matrix: V[day][damage]
    vector<vector<double>> V(max_days + 1, vector<double>(max_damages + 1, 0.0));

    auto start_time = Clock::now();

    for (int iter = 0; iter < max_iterations; ++iter) {
        cout << "doing iteration " << iter << "\n";
        double max_diff = 0.0;

        for (int d = max_days; d >= 0; --d) {
            for (int dmg = 0; dmg <= max_damages; ++dmg) {
                double new_value;

                if (dmg == max_damages) {
                    // Failure state: reset to (0, 0)
                    new_value = R + gamma * V[0][0];
                }
                else {
                    // Maintenance: reset day, same damage
                    double cost_maintenance = C + gamma * V[0][dmg];

                    // No maintenance
                    double p_dmg = get_p(d);
                    double cost_no_maintenance;

                    if (dmg + 1 < max_damages) {
                        double val_same = V[min(d + 1, max_days)][dmg];
                        double val_more = V[min(d + 1, max_days)][dmg + 1];
                        cost_no_maintenance = gamma * ((1.0 - p_dmg) * val_same + p_dmg * val_more);
                    }
                    else {
                        double val_same = V[min(d + 1, max_days)][dmg];
                        double fail_value = R + gamma * V[0][0];
                        cost_no_maintenance = gamma * ((1.0 - p_dmg) * val_same + p_dmg * fail_value);
                    }

                    new_value = min(cost_maintenance, cost_no_maintenance);
                }

                double diff = abs(new_value - V[d][dmg]);
                V[d][dmg] += omega * (new_value - V[d][dmg]);  // SOR update
                max_diff = max(max_diff, diff);

                if (!isfinite(V[d][dmg])) {
                    cerr << "Diverged at iteration " << iter << " (d=" << d << ", dmg=" << dmg << ")\n";
                    return 1;
                }
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
    cout << "\nTime elapsed during SOR Value Iteration: " << elapsed.count() << " seconds\n";

    // Output part of the value function
    cout << "\nValue Function Table (rows = days, columns = damage count):\n";
    cout << "Days";
    for (int dmg = 0; dmg < max_damages; ++dmg)
        cout << setw(12) << dmg << "dmg";
    cout << "\n" << string(12 + 12 * max_damages, '-') << "\n";

    for (int d = 0; d <= 360; d+=10) {
        cout << setw(5) << d;
        for (int dmg = 0; dmg < max_damages; ++dmg) {
            cout << setw(12) << V[d][dmg];
        }
        cout << "\n";
    }

    return 0;
}