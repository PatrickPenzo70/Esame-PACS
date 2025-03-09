#include <iostream>
#include <fstream>
#include <sstream>
#include <unordered_map>
#include <cmath>
#include <vector>
#include <random>
#include <algorithm>
#include <iomanip>

// Function to read simulation parameters from CSV file
bool readParameters(const std::string& filename, std::unordered_map<std::string, double>& params) {
    std::ifstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Error: Unable to open file " << filename << "\n";
        return false;
    }

    std::string line, key;
    double value;
    
    std::getline(file, line); // Skip header
    while (std::getline(file, line)) {
        std::stringstream ss(line);
        if (std::getline(ss, key, ',') && ss >> value) {
            params[key] = value;
        }
    }
    file.close();
    return true;
}

// Function to calculate fluid velocity
void fluid_velocity(double x, double y, double v_avg, double H, double& vx, double& vy) {
    vx = 1.5 * v_avg * (1 - (y * y) / (H * H));
    vy = 0;
}

// Function to calculate particle force
void particle_force(double x, double y, double g, double m, double& fx, double& fy) {
    fx = 0;
    fy = -m * g;
}

int main() {
    std::unordered_map<std::string, double> params;
    if (!readParameters("Data.csv", params)) {
        return 1;
    }

    // Assign values from the map
    double v_avg = params["v_avg"];
    int Np = static_cast<int>(params["Np"]);
    double dt = params["dt"];
    double T = params["T"];
    double g = params["g"];
    double m = params["m"];
    double H = params["H"];
    double KT = params["KT"];
    double mob = 1.0;
    double D = KT * mob;
    int Nt = std::ceil(T / dt);

    std::vector<double> x(Np), y(Np), vx(Np), vy(Np), fx(Np), fy(Np), dxb(Np), dyb(Np);

    std::random_device rd;
    std::mt19937 gen(rd());
    std::normal_distribution<> normal_x(0.0, 1.0);
    std::normal_distribution<> normal_y(0.0, H / 3.0);
    std::normal_distribution<> disbrownian(0.0, std::sqrt(2 * D * dt));

    // Initialize particle positions
    for (int n = 0; n < Np; ++n) {
        x[n] = normal_x(gen);
        y[n] = normal_y(gen);
    }

    double t = 0;
    int iff = 0;

    // Main simulation loop
    for (int it = 0; it < Nt; ++it) {
        std::ostringstream oss;
        oss << "particle_position_" << std::setw(5) << std::setfill('0') << iff++ << ".csv";
        std::ofstream outFile(oss.str());
        
        if (!outFile.is_open()) {
            std::cerr << "Error opening output file." << std::endl;
            return 1;
        }
        
        outFile << "Time,Ball_X,Ball_Y\n";

        for (int n = 0; n < Np; ++n) {
            fluid_velocity(x[n], y[n], v_avg, H, vx[n], vy[n]);
            particle_force(x[n], y[n], g, m, fx[n], fy[n]);

            dxb[n] = disbrownian(gen);
            dyb[n] = disbrownian(gen);

            x[n] += (vx[n] + mob * fx[n]) * dt + dxb[n];
            y[n] += (vy[n] + mob * fy[n]) * dt + dyb[n];

            // Reflecting boundary condition
            if (y[n] < -H) y[n] = -H - (y[n] + H);
            if (y[n] > H) y[n] = H - (y[n] - H);
            
            outFile << t << "," << x[n] << "," << y[n] << "\n";
        }
        
        t += dt;
        outFile.close();
    }

    std::cout << "Simulation complete. Data saved to individual CSV files.\n";
    return 0;
}

