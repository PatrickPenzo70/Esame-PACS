// MPM-style static melting demo: one Cr2O3 particle in hot Ar+H2 gas
// Zero velocities (gas and particle). Only heat diffusion + latent heat + simple phase kinetics.
// C++17 single file, finite-difference diffusion on grid, MPM used to represent the particle geometry
// and to exchange fields between particles and the Eulerian grid.
// -----------------------------------------------------------------------------
// DISCLAIMER: Educational prototype. Calibrate material properties and kinetics
// before using for quantitative predictions.

#include <iostream>
#include <vector>
#include <algorithm>
#include <fstream> 
#include <cmath>
#include <random>

using namespace std;

struct Vec2 {
    double x = 0, y = 0;

    // Constructors
    Vec2() {}
    Vec2(double X, double Y) : x(X), y(Y) {}

    // Operator overloads
    Vec2 operator+(const Vec2& o) const { return {x + o.x, y + o.y}; }
    Vec2 operator-(const Vec2& o) const { return {x - o.x, y - o.y}; }
    Vec2 operator*(double s) const { return {x * s, y * s}; }

    // Optional extras
    double dot(const Vec2& o) const { return x * o.x + y * o.y; }
    double norm() const { return std::sqrt(x * x + y * y); }
    Vec2 normalized() const {
        double n = norm();
        return n > 0 ? Vec2(x / n, y / n) : Vec2(0, 0);
    }
};

// Clamp an integer 'a' between [lo, hi]
static inline int clampi(int a, int lo, int hi) {
    return std::max(lo, std::min(a, hi));
}

// Convert 2D indices (i, j) to a 1D array index
static inline int I(int i, int j, int NX) {
    return j * NX + i;
}

inline Vec2 operator*(double s, const Vec2& v) { return {v.x * s, v.y * s}; }

struct Params{
    // Grid
    int NX=160, NY=160; // cells
    double dx=1e-5; // 10 microns per cell
    double dt=2e-8; // seconds
    int steps=10000; // time steps
    int output_every=100; // CSV cadence

    // Gas (Ar+H2 75/25) – placeholders; tune as needed
    double rho_g=0.2;       // kg/m^3
    double cp_g=1500.0;     // J/kg-K
    double k_g=10.0;        // W/m-K
    double T_inf=15000.0;   // K (ambient plasma temperature)

    // Chromium oxide (Cr2O3)
    double rho_cr=5200.0;   // kg/m^3
    double cp_s=800.0;      // J/kg-K
    double cp_l=1000.0;     // J/kg-K
    double k_s=5.0;         // W/m-K
    double k_l=3.0;         // W/m-K
    double Tm=2435.0;       // K
    double Lm=4e5;          // J/kg

    // Phase-field kinetics (no gradient term; local kinetics for simplicity)
    double Kphi=1e-3;       // 1/s/K   (phi rate ~ Kphi * max(T-Tm,0))

    // Particle
    Vec2 center{0.0008,0.0008}; // domain ~ 1.6 mm square
    double R = 100e-6;   // 100 micron radius
    int nParticles = 10000;  // material points sampling the disc

    // Numerics
    bool dirichlet_gas=true; // enforce T=T_inf in gas each step (gas at rest, infinite reservoir)
};

struct Cell{
    double T=300.0;     // temperature [K]
    double rhoCp=0.0;   // volumetric heat capacity [J/m^3-K]
    double k=0.0;       // thermal conductivity [W/m-K]
    double phi=0.0;     // melt fraction in this cell (only where particle exists)
    double mask=0.0;    // 0 gas, 1 particle (Cr2O3) fraction [0,1]
};

struct Particle{
    Vec2 x; double m; double vol; double T; double phi; // 0..1
};

struct Sim{
    Params P;
    vector<Cell> grid; 
    vector<Particle> pts;

    Sim(const Params& p):P(p){ grid.resize(P.NX*P.NY); }

    void seed_particle(){
        // Poisson disk-ish random fill inside circle
        std::mt19937 rng(42);
        std::uniform_real_distribution<double> U(0,1);
        int tries = 0; pts.clear();
        double area = M_PI*P.R*P.R;
        double V = area*P.dx; // pseudo-thickness = dx (2D slab)
        double mpt = P.rho_cr * V / (double)P.nParticles; // average mass per point
        while((int)pts.size()<P.nParticles && tries < P.nParticles*50){
            tries++;
            double r = P.R*sqrt(U(rng));
            double th = 2*M_PI*U(rng);
            Vec2 pos{ P.center.x + r*cos(th), P.center.y + r*sin(th) };
            Particle p; p.x=pos; p.vol=V/(double)P.nParticles; p.m=mpt; p.T=300.0; p.phi=0.0; pts.push_back(p);
        }
    }

    void reset_material_fields(){
        // Zero fields
        for(auto &c: grid){ 
        
                c.rhoCp=0; 
                c.k=0; 
                c.mask=0; 
                c.phi=0; 
                
        }
        // Gas defaults
        for(auto &c: grid){ 
        
                c.T = P.T_inf; 
                
        }
        // Accumulate particle to grid (nearest cell)
        for(const auto &p: pts){
            int i = clampi((int)floor(p.x.x/P.dx),0,P.NX-1);
            int j = clampi((int)floor(p.x.y/P.dx),0,P.NY-1);
            Cell &C = grid[I(i,j,P.NX)];
            C.mask += 1.0; // count points
            C.T = min(C.T, p.T); // colder particle imposes initial lower T locally
            C.phi += p.phi;
        }
        // Convert counts to fractions and set properties
        for(auto &c: grid){
            if(c.mask>0){
                double f = min(1.0, c.mask/10.0); // crude fractional fill (points per cell / 10)
                c.mask = f;
                c.phi /= max(1.0, c.mask*10.0);
                // mixture inside cell (simple blend):
                double rhoCp_cr = P.rho_cr*((1.0-c.phi)*P.cp_s + c.phi*P.cp_l);
                double k_cr = (1.0-c.phi)*P.k_s + c.phi*P.k_l;
                // volumetric averaging with gas by mask fraction
                c.rhoCp = f * rhoCp_cr + (1.0-f) * P.rho_g*P.cp_g;
                c.k     = f * k_cr     + (1.0-f) * P.k_g;
                // start from colder temperature if particle present
                c.T = min(c.T, 300.0);
            } else {
                c.rhoCp = P.rho_g*P.cp_g;
                c.k     = P.k_g;
                c.T     = P.T_inf; // gas is hot reservoir
            }
        }
    }

    void diffuse_and_melt_step(){
        // Compute dphi for particle cells (local kinetics), and latent source S [W/m^3]
        vector<double> dphi(grid.size(),0.0); 
        vector<double> S(grid.size(),0.0);
        for(int j=0;j<P.NY;++j){
            for(int i=0;i<P.NX;++i){
                Cell &c = grid[I(i,j,P.NX)];
                if(c.mask>0.0){
                    if(c.T > P.Tm && c.phi < 1.0){
                        double rate = P.Kphi * (c.T - P.Tm) * (1.0 - c.phi); // 0..∞
                        // clamp so phi in [0,1]
                        rate = min(rate, (1.0 - c.phi)/P.dt);
                        dphi[I(i,j,P.NX)] = rate;
                        // latent heat sink (melting consumes energy): S = -rho*Lm*dphi
                        double rho_cr_vol = P.rho_cr * c.mask; // crude fraction to volumetric
                        S[I(i,j,P.NX)] = - rho_cr_vol * P.Lm * rate; 
                    }
                }
            }
        }
        // Diffusion: ∂(ρCp T)/∂t = ∇·(k ∇T) + S  => explicit Euler on T
        vector<double> Tnew(grid.size());
        auto k_face = [&](double kL,double kR){ // harmonic average
            return 2.0*kL*kR / max(1e-12, (kL+kR));
        };
        for(int j=0;j<P.NY;++j){
            for(int i=0;i<P.NX;++i){
                int id = I(i,j,P.NX);
                Cell &c = grid[id];
                double flux = 0.0;
                // 4-neighbor variable-k Laplacian
                if(i>0){
                    int il = I(i-1,j,P.NX); double kf = k_face(c.k, grid[il].k);
                    flux += kf * (grid[il].T - c.T) / (P.dx*P.dx);
                }
                if(i<P.NX-1){
                    int ir = I(i+1,j,P.NX); double kf = k_face(c.k, grid[ir].k);
                    flux += kf * (grid[ir].T - c.T) / (P.dx*P.dx);
                }
                if(j>0){
                    int jd = I(i,j-1,P.NX); double kf = k_face(c.k, grid[jd].k);
                    flux += kf * (grid[jd].T - c.T) / (P.dx*P.dx);
                }
                if(j<P.NY-1){
                    int ju = I(i,j+1,P.NX); double kf = k_face(c.k, grid[ju].k);
                    flux += kf * (grid[ju].T - c.T) / (P.dx*P.dx);
                }
                double rhs = flux + S[id];
                double Tn = c.T + P.dt * rhs / max(1e-12, c.rhoCp);
                Tnew[id] = Tn;
            }
        }
        // Update T, phi; enforce gas boundary if needed
        for(int j=0;j<P.NY;++j){
            for(int i=0;i<P.NX;++i){
                int id = I(i,j,P.NX);
                Cell &c = grid[id];
                c.T = Tnew[id];
                if(c.mask>0.0){
                    c.phi = min(1.0, c.phi + P.dt * dphi[id]);
                } else if(P.dirichlet_gas){
                    c.T = P.T_inf; // gas held at hot temperature
                }
            }
        }
        // Map grid back to particles (sample nearest cell)
        for(auto &p: pts){
            int i = clampi((int)floor(p.x.x/P.dx),0,P.NX-1);
            int j = clampi((int)floor(p.x.y/P.dx),0,P.NY-1);
            const Cell &C = grid[I(i,j,P.NX)];
            p.T = C.T;
            p.phi = C.phi;
        }
    }

    void write_csv(int s){
        string base = string("static_out_") + to_string(s);
        ofstream fg(base+"_grid.csv");
        fg << "i,j,x,y,T,phi,mask\n";
        for(int j=0;j<P.NY;++j){
            for(int i=0;i<P.NX;++i){
                double x=(i+0.5)*P.dx, y=(j+0.5)*P.dx;
                const Cell &c = grid[I(i,j,P.NX)];
                fg<<i<<","<<j<<","<<x<<","<<y<<","<<c.T<<","<<c.phi<<","<<c.mask<<"\n";
            }
        }
        ofstream fp(base+"_particles.csv");
        fp << "id,x,y,T,phi\n";
        for(size_t k=0;k<pts.size();++k){
            fp<<k<<","<<pts[k].x.x<<","<<pts[k].x.y<<","<<pts[k].T<<","<<pts[k].phi<<"\n";
        }
        cerr<<"Wrote "<<base<<"_grid.csv and _particles.csv\n";
    }

    void run(){
        seed_particle();
        reset_material_fields();
        for(int s=0;s<P.steps;++s){
            diffuse_and_melt_step();
            if(s%P.output_every==0){
                double phim=0,Tavg=0; for(auto &p: pts){ phim+=p.phi; Tavg+=p.T; }
                phim/=max(1.0,(double)pts.size()); Tavg/=max(1.0,(double)pts.size());
                cout<<"step "<<s<<" phi_avg="<<phim<<" T_avg="<<Tavg<<"\n";
                write_csv(s);
            }
        }
    }
};

int main(){
    Params P; 
    Sim sim(P);
    sim.run();
    return 0;
}

