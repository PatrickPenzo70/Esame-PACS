#include <mpi.h>
#include <cmath>
#include <iostream>
#include <map>
#include <random>
#include <string>
#include <vector>
#include <algorithm>
#include <fstream>
#include <limits>

#include "quadgrid_cpp.h"
#include "particles.h"

using dist_vec    = std::vector<double>;
using grid_t      = quadgrid_t<dist_vec>;
using idx_t       = grid_t::idx_t;
using particles_t = ::particles_t;

/* ---------------- Costanti materiali / modello ---------------- */

// Ambiente / plasma
constexpr double T_AMB    = 300.0;     // K
constexpr double T_PLASMA = 15000.0;   // K

// Ti-6Al-4V
constexpr double T_M = 1928.0;         // K (fusione)

// Ti solid (usato anche per liquido in questo modello semplificato)
constexpr double RHO_TI_SOLID = 4420.0;
constexpr double CP_TI_SOLID  = 560.0;
constexpr double K_TI_SOLID   = 7.2;

// Argon (gas) -> NON usato (psi1=0 sempre)
constexpr double RHO_AR = 0.35;
constexpr double CP_AR  = 520.0;
constexpr double K_AR   = 0.18;

// Phase-field
constexpr double DW_SCALE     = 0.05;
constexpr double LATENT_SCALE = 0.002;
constexpr double LATENT_HEAT  = 3.64e5;   // J/kg
constexpr double G_PHI        = 1.0e4;
constexpr double M_PHI        = 1.0e-7;
constexpr double XI_PHI2      = 1.0e-12;

// CFL
constexpr double CFL_SAFE = 0.20;

// Scambio termico
constexpr double H_PLASMA = 1.0e6;     // W/m^2/K
constexpr double H_AIR    = 20.0;      // W/m^2/K

// Radiazione
constexpr double EMISS    = 0.7;
constexpr double SIGMA_SB = 5.670374419e-8;

// Clamp temperatura
constexpr double T_MAX_CLAMP = 5.0 * T_PLASMA;

// 2D estruso: spessore fisico
constexpr double THICKNESS = 1.0e-7;   // m

// Soglie numeriche
constexpr double EPS_W   = 1e-16;
constexpr double EPS_MCP = 1e-18;

// Maschera diffusione
constexpr double FILL_MIN = 0.10;  // diffondi solo se fill>10%

// G2P additivo: soglia robusta sul peso
constexpr double WP_MIN = 0.20;    // se wp<0.2 => fallback su T_old

/* ---------------- Funzioni di potenziale PF ---------------- */

inline double dF_phi(double phi) { return phi*phi*phi - phi; }

inline double dP_phi(double phi) {
    double phi2 = phi * phi;
    double phi4 = phi2 * phi2;
    return -(1.0 / 16.0) * (15.0 * phi4 - 30.0 * phi2 + 15.0);
}

/* ---------------- Proiezione entalpia solido/liquido (NO gas) ----------------
   Convenzione:
   - psi1 = 0 sempre
   - psi2 = frazione solida
   - psi3 = frazione liquida
   - mcpT su particella = ENERGIA totale E [J] (sensibile + latente)
*/
inline void project_enthalpy_solid_liquid(particles_t &pts, idx_t p)
{
    const double mp = pts.dp("mp", p);
    double &T       = pts.dp("T", p);
    double &E       = pts.dp("mcpT", p); // energia totale [J]

    const double mcp = mp * CP_TI_SOLID;
    pts.dp("mcp", p) = mcp;

    // NO gas sempre
    pts.dp("psi1", p) = 0.0;

    // limiti energia
    E = std::max(E, 0.0);

    // soglie energia per il plateau a T_M
    const double E_solid  = mcp * T_M;              // T = T_M, fL=0
    const double E_liquid = E_solid + mp*LATENT_HEAT; // T = T_M, fL=1

    double fL = 0.0;

    if (E <= E_solid) {
        fL = 0.0;
        T  = E / std::max(mcp, 1e-30);
    } else if (E < E_liquid) {
        T  = T_M;
        fL = (E - E_solid) / std::max(mp*LATENT_HEAT, 1e-30);
    } else {
        fL = 1.0;
        T  = T_M + (E - E_liquid) / std::max(mcp, 1e-30);
    }

    fL = std::clamp(fL, 0.0, 1.0);
    T  = std::clamp(T, T_AMB, T_MAX_CLAMP);

    pts.dp("psi3", p) = fL;        // liquido
    pts.dp("psi2", p) = 1.0 - fL;  // solido
    pts.dp("rho",  p) = RHO_TI_SOLID;

    // rendi E coerente col clamp
    if (T < T_M) {
        E = mcp * T;
    } else if (fL < 1.0) {
        E = E_solid + fL * mp * LATENT_HEAT;
    } else {
        E = E_liquid + mcp * (T - T_M);
    }
}

inline void enthalpy_to_T_fL_node(
    double E, double mcp,
    double &T, double &fL
){
    E   = std::max(E, 0.0);
    mcp = std::max(mcp, 1e-30);

    const double mp = mcp / CP_TI_SOLID;

    const double E_solid  = mcp * T_M;
    const double E_liquid = E_solid + mp * LATENT_HEAT;

    if (E <= E_solid) {
        fL = 0.0;
        T  = E / mcp;
    } else if (E < E_liquid) {
        T  = T_M;
        fL = (E - E_solid) / std::max(mp * LATENT_HEAT, 1e-30);
    } else {
        fL = 1.0;
        T  = T_M + (E - E_liquid) / mcp;
    }

    fL = std::clamp(fL, 0.0, 1.0);
    T  = std::clamp(T,  T_AMB, T_MAX_CLAMP);
}


/* ---------------- Laplaciano su QuadGrid ---------------- */

double laplacian(const dist_vec &f, const grid_t &grid, idx_t idx) {
    idx_t r  = grid.gind2row(idx);
    idx_t c  = grid.gind2col(idx);
    idx_t nr = grid.num_rows();
    idx_t nc = grid.num_cols();
    if (r == 0 || c == 0 || r == nr - 1 || c == nc - 1) return 0.0;

    auto s = [&](idx_t rr, idx_t cc) { return grid.sub2gind(rr, cc); };
    double hx = grid.hx(), hy = grid.hy();
    double c0 = f[idx];
    double L = f[s(r, c - 1)];
    double R = f[s(r, c + 1)];
    double D = f[s(r - 1, c)];
    double U = f[s(r + 1, c)];
    return (L - 2*c0 + R) / (hx*hx) + (D - 2*c0 + U) / (hy*hy);
}

/* ---------------- Init: disco 2D diametro 50 µm ---------------- */

void init_particle_disk_50um(particles_t &pts, double cx, double cy, double R) {
    idx_t Np = pts.num_particles;

    // Volume totale (2D area * thickness) e volume per particella
    double area_disk = M_PI * R * R;
    double Vtot      = area_disk * THICKNESS;
    double Vp0       = Vtot / double(Np);

    std::mt19937_64 rng(42);
    std::uniform_real_distribution<double> U(0.0, 1.0);

    for (idx_t p = 0; p < Np; ++p) {
        double theta = 2.0 * M_PI * U(rng);
        double r     = R * std::sqrt(U(rng));  // uniforme in area
        pts.x[p] = cx + r * std::cos(theta);
        pts.y[p] = cy + r * std::sin(theta);

        pts.dp("phi", p) = 1.0;

        // stato iniziale: solido a T_AMB
        pts.dp("T",    p) = T_AMB;
        pts.dp("psi1", p) = 0.0;   // gas OFF sempre
        pts.dp("psi3", p) = 0.0;   // liquido
        pts.dp("psi2", p) = 1.0;   // solido

        pts.dp("rho",  p) = RHO_TI_SOLID;
        pts.dp("w",    p) = 1.0;

        pts.dp("one",  p) = 1.0;
        pts.dp("wpE",  p) = 0.0;

        // volume e massa
        pts.dp("Vp", p) = Vp0;
        double mp = RHO_TI_SOLID * Vp0;   // massa costante per particella
        pts.dp("mp", p) = mp;

        // energia totale iniziale E = mcp*T
        double mcp = mp * CP_TI_SOLID;
        pts.dp("mcp",  p) = mcp;
        pts.dp("mcpT", p) = mcp * T_AMB;  // E [J]

        pts.dp("in_jet", p) = 0.0;
        
        // marca particelle che stanno sulla corona esterna e lato sinistro (jet)
        double dx = pts.x[p] - cx;
        double dy = pts.y[p] - cy;
        double rr = std::sqrt(dx*dx + dy*dy);

        double dr = 0.5e-6; // usa hx/2 (qui hx=1e-6), metto fisso per init
        double theta_j = std::atan2(dy, dx); // [-pi, pi]

        // lato sinistro: cos(theta) < 0  (dx<0) e anche dentro semicerchio sinistro
        bool left = (dx < 0.0);

        // superficie: corona esterna (tolleranza più larga)
        bool surf = (rr > R - 2.0*dr && rr < R + 2.0*dr);

        pts.dp("in_jet", p) = (left && surf) ? 1.0 : 0.0;
   
    }

    pts.init_particle_mesh();
    pts.build_mass();
}

/* ---------------- Writer VTK griglia (.vts) ---------------- */

void write_vtk_vts(
    const grid_t &grid,
    const std::map<std::string, dist_vec> &gvars,
    int step
) {
    std::string fname = "grid_" + std::to_string(step) + ".vts";
    std::ofstream f(fname);
    if (!f) return;

    idx_t nx = grid.num_cols();
    idx_t ny = grid.num_rows();
    double hx = grid.hx(), hy = grid.hy();

    f << "<?xml version=\"1.0\"?>\n";
    f << "<VTKFile type=\"StructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n";
    f << "  <StructuredGrid WholeExtent=\"0 " << nx-1 << " 0 " << ny-1 << " 0 0\">\n";
    f << "    <Piece Extent=\"0 " << nx-1 << " 0 " << ny-1 << " 0 0\">\n";

    f << "      <Points>\n";
    f << "        <DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">\n";
    for (idx_t r = 0; r < ny; ++r)
        for (idx_t c = 0; c < nx; ++c)
            f << c*hx << " " << r*hy << " 0 ";
    f << "\n        </DataArray>\n";
    f << "      </Points>\n";

    f << "      <PointData Scalars=\"scalars\">\n";
    auto dump = [&](const std::string &name){
        auto it = gvars.find(name);
        if (it == gvars.end()) return;
        f << "        <DataArray type=\"Float64\" Name=\"" << name
          << "\" NumberOfComponents=\"1\" format=\"ascii\">\n";
        for (double v : it->second) f << v << " ";
        f << "\n        </DataArray>\n";
    };

    dump("phi");
    dump("T");
    dump("mcp");
    dump("mcpT");
    dump("w");

    f << "      </PointData>\n";
    f << "    </Piece>\n";
    f << "  </StructuredGrid>\n";
    f << "</VTKFile>\n";

    std::cout << "  🔷 VTK grid exported: " << fname << "\n";
}

/* ---------------- Writer VTK particelle (.vtp) ---------------- */

void write_particles_vtp(const particles_t &pts, int step) {
    std::string fname = "particles_" + std::to_string(step) + ".vtp";
    std::ofstream f(fname);
    if (!f) return;

    const idx_t N = pts.num_particles;

    f << "<?xml version=\"1.0\"?>\n";
    f << "<VTKFile type=\"PolyData\" version=\"0.1\" byte_order=\"LittleEndian\">\n";
    f << "  <PolyData>\n";
    f << "    <Piece NumberOfPoints=\"" << N << "\" NumberOfVerts=\"" << N << "\">\n";

    f << "      <Points>\n";
    f << "        <DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">\n";
    for (idx_t p = 0; p < N; ++p)
        f << pts.x[p] << " " << pts.y[p] << " 0 ";
    f << "\n        </DataArray>\n";
    f << "      </Points>\n";

    f << "      <Verts>\n";
    f << "        <DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">\n";
    for (idx_t p = 0; p < N; ++p) f << p << " ";
    f << "\n        </DataArray>\n";
    f << "        <DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">\n";
    for (idx_t p = 0; p < N; ++p) f << (p+1) << " ";
    f << "\n        </DataArray>\n";
    f << "      </Verts>\n";

    f << "      <PointData Scalars=\"scalars\">\n";

    auto dump_prop = [&](const std::string &name){
        f << "        <DataArray type=\"Float64\" Name=\"" << name << "\" format=\"ascii\">\n";
        for (idx_t p = 0; p < N; ++p) f << pts.dp(name, p) << " ";
        f << "\n        </DataArray>\n";
    };

    dump_prop("T");
    dump_prop("phi");
    dump_prop("psi2");
    dump_prop("psi3");
    dump_prop("mcp");
    dump_prop("mcpT");

    // phase_id: 0=solido, 1=liquido (psi3>0.5)
    f << "        <DataArray type=\"Int32\" Name=\"phase_id\" format=\"ascii\">\n";
    for (idx_t p = 0; p < N; ++p) {
        int pid = (pts.dp("psi3", p) > 0.5) ? 1 : 0;
        f << pid << " ";
    }
    f << "\n        </DataArray>\n";

    f << "      </PointData>\n";
    f << "    </Piece>\n";
    f << "  </PolyData>\n";
    f << "</VTKFile>\n";

    std::cout << "  🟠 VTK particles exported: " << fname << "\n";
}

/* ---------------- BC convettiva + radiativa su "superficie" ----------------
   Area totale (2πR*THICKNESS) distribuita sulle particelle di superficie nel jet.
   Aggiorna SOLO energia (mcpT=E), poi proietta entalpia.
*/

void apply_plasma_bc_to_surface_points(
    particles_t &pts,
    double dt,
    double cx, double cy,
    double R,
    double hx,
    double time,
    double plasma_on_time
){
    const bool plasma_on = (time <= plasma_on_time);
    const double Tinf = plasma_on ? T_PLASMA : T_AMB;
    const double h    = plasma_on ? H_PLASMA : H_AIR;

    const double dr = 0.5 * hx;

    constexpr double JET_HALF_ANGLE_DEG = 90.0;
    const double cos0 = std::cos(JET_HALF_ANGLE_DEG * M_PI / 180.0);

    auto in_jet = [&](double dx, double r)->bool{
        double nx = dx / (r + 1e-30);   // cos(theta) rispetto a +x
        return (nx <= -cos0);           // lato sinistro (~-x)
    };

    idx_t Nsurf = 0;
    for (idx_t p = 0; p < pts.num_particles; ++p) {
        double dxp = pts.x[p] - cx;
        double dyp = pts.y[p] - cy;
        double rr  = std::sqrt(dxp*dxp + dyp*dyp);

        if (rr < R - dr || rr > R + 1e-12) continue;
        if (!in_jet(dxp, rr)) continue;
        Nsurf++;
    }
    if (Nsurf == 0) return;

    static int dbg_counter = 0;
    if ((dbg_counter++ % 200) == 0) {
        std::cout << "  [BC DBG] time=" << time
                  << " plasma_on=" << plasma_on
                  << " Nsurf=" << Nsurf << "\n";
    }

    const double theta0  = JET_HALF_ANGLE_DEG * M_PI / 180.0;
    const double frac    = theta0 / M_PI;
    const double A_total = frac * (2.0 * M_PI * R) * THICKNESS;
    const double A_p     = A_total / double(Nsurf);

    for (idx_t p = 0; p < pts.num_particles; ++p) {
        double dxp = pts.x[p] - cx;
        double dyp = pts.y[p] - cy;
        double rr  = std::sqrt(dxp*dxp + dyp*dyp);

        if (rr < R - dr || rr > R + 1e-12) continue;
        if (!in_jet(dxp, rr)) continue;

        const double mp  = pts.dp("mp", p);
        const double mcp = mp * CP_TI_SOLID;
        pts.dp("mcp", p) = mcp;

        double &T = pts.dp("T", p);
        double &E = pts.dp("mcpT", p);

        // (consigliato) riallinea T con E prima di usare T nei flussi
        project_enthalpy_solid_liquid(pts, p);

        // flussi [W/m^2]
        const double qconv = h * (Tinf - T);
        const double qrad  = EMISS * SIGMA_SB * (std::pow(Tinf,4) - std::pow(T,4));
        double qin = qconv + qrad;

        constexpr double QIN_MAX = 1.0e8;
        qin = std::clamp(qin, -QIN_MAX, QIN_MAX);

        const double dE = qin * A_p * dt;
        E += dE;

        static int dbg2 = 0;
        if ((dbg2++ % 200) == 0) {
            std::cout << "  [BC DBG2] t=" << time
                      << " plasma_on=" << plasma_on
                      << " p=" << p
                      << " T_before=" << T
                      << " qin=" << qin
                      << " A_p=" << A_p
                      << " dE=" << dE
                      << " E_after=" << E
                      << "\n";
        }

        project_enthalpy_solid_liquid(pts, p);
    }
}


/* ---------------- 1 passo: P2G -> PF/termica -> G2P -> BC ---------------- */

void do_one_mpm_timestep(
    const grid_t &grid,
    particles_t &pts,
    std::map<std::string, dist_vec> &gvars,
    double dt,
    double time,
    double cx, double cy, double R,
    double plasma_on_time
){
    const idx_t nn = grid.num_global_nodes();
    const idx_t ny = grid.num_rows();
    const idx_t nx = grid.num_cols();

    const double Vcv = grid.hx() * grid.hy() * THICKNESS;

    // 0) particelle: NO gas. Aggiorna SOLO mcp e clamp T.
    //    NON ricostruire mcpT qui: è energia conservata (include BC)!
    for (idx_t p = 0; p < pts.num_particles; ++p) {
      pts.dp("psi1", p) = 0.0;

      const double mp  = pts.dp("mp", p);
      const double mcp = mp * CP_TI_SOLID;
      pts.dp("mcp", p) = mcp;
    
      // clamp solo diagnostico (ma NON toccare mcpT)
      pts.dp("T", p) = std::clamp(pts.dp("T", p), T_AMB, T_MAX_CLAMP);

      pts.dp("rho", p) = RHO_TI_SOLID;
      pts.dp("one", p) = 1.0;
    }
    
    // 0b) plasma BC SUBITO, così entra nel P2G di questo timestep
    apply_plasma_bc_to_surface_points(
         pts, dt, cx, cy, R, grid.hx(), time, plasma_on_time
    );
    
    // 1) azzera griglia
    for (auto &kv : gvars) kv.second.assign(nn, 0.0);

    // 2) P2G: quantità conservative (energia via mcp/mcpT)
    pts.p2g(
        gvars,
        {"phi","rho","w","mcp","mcpT"},
        {"phi","rho","w","mcp","mcpT"}
    );

    dist_vec &phi  = gvars["phi"];
    dist_vec &rho  = gvars["rho"];
    dist_vec &w    = gvars["w"];
    dist_vec &mcp  = gvars["mcp"];
    dist_vec &mcpT = gvars["mcpT"];
    dist_vec &T    = gvars["T"];

    // 3) normalizzazione nodale + ricostruzione T = mcpT/mcp
    for (idx_t i = 0; i < nn; ++i) {
        if (w[i] > EPS_W) {
            phi[i] /= w[i];
            rho[i] /= w[i];
        }

        if (mcp[i] > EPS_MCP) {
    double fL_i = 0.0;
    enthalpy_to_T_fL_node(mcpT[i], mcp[i], T[i], fL_i);
    } else {
          T[i] = T_AMB;
    }
          T[i] = std::clamp(T[i], T_AMB, T_MAX_CLAMP);
    }

    dist_vec phi_new(nn);
    dist_vec mcpT_new = mcpT;   // energia nuova

    // 4) update PF + termica (termica: diffusione su T con proprietà costanti Ti)
    for (idx_t r = 0; r < ny; ++r) {
        for (idx_t c = 0; c < nx; ++c) {
            idx_t i = grid.sub2gind(r, c);

            const bool has_geom = (w[i] > EPS_W);

            // PF
            if (has_geom) {
                const double lap_phi = laplacian(phi, grid, i);
                const double chemical = DW_SCALE * G_PHI * dF_phi(phi[i]);
                const double latent   = LATENT_SCALE * LATENT_HEAT * (T[i] - T_M) / T_M * dP_phi(phi[i]);
                const double rhs_phi  = -M_PHI * (chemical + latent - XI_PHI2 * lap_phi);
                phi_new[i] = std::clamp(phi[i] + dt * rhs_phi, -1.0, 1.0);
            } else {
                phi_new[i] = phi[i];
            }

            // TERMICA
            const double k_eff   = K_TI_SOLID;
            const double rho_eff = RHO_TI_SOLID;
            const double cp_eff  = CP_TI_SOLID;

            const double mcp_full = std::max(rho_eff * cp_eff * Vcv, 1e-30);
            const double fill     = (mcp_full > 0.0) ? (mcp[i] / mcp_full) : 0.0;

            if (fill > FILL_MIN) {
              const double lap_T = laplacian(T, grid, i);

            // dE/dt = k * V * fill * lap(T)
              const double dE_dt = (k_eff * Vcv * fill) * lap_T;

            mcpT_new[i] = std::max(0.0, mcpT[i] + dt * dE_dt);
            } else {
                  mcpT_new[i] = mcpT[i];
            }
        }
    }

    phi  = phi_new;
    mcpT = mcpT_new;

    // 4b) mantieni coerente l'energia nodale con T aggiornato (mcp costante)
    for (idx_t i = 0; i < nn; ++i) {
    if (mcp[i] > EPS_MCP) {
        double fL_i = 0.0;
        enthalpy_to_T_fL_node(mcpT[i], mcp[i], T[i], fL_i);
    } else {
        T[i] = T_AMB;
    }
    T[i] = std::clamp(T[i], T_AMB, T_MAX_CLAMP);
    }


    // 5) one su griglia: solo nodi con fill>FILL_MIN (per wpE coerente)
    dist_vec &one = gvars["one"];
    one.assign(nn, 0.0);
    for (idx_t i = 0; i < nn; ++i) {
        const double mcp_full = std::max(RHO_TI_SOLID * CP_TI_SOLID * Vcv, 1e-30);
        const double fill     = mcp[i] / mcp_full;
        if (fill > FILL_MIN) one[i] = 1.0;
    }

    // 6) G2P additivo: mappa phi e T, normalizza con wpE
    dist_vec E_old(pts.num_particles);
    for (idx_t p = 0; p < pts.num_particles; ++p) E_old[p] = pts.dp("mcpT", p);
    
    for (idx_t p = 0; p < pts.num_particles; ++p) {
        pts.dp("phi",  p) = 0.0;
        pts.dp("mcpT", p) = 0.0;   // accumulo energia da griglia
        pts.dp("wpE",  p) = 0.0;   // peso per energia (nuovo campo particellare)
    }

    pts.g2p(gvars, {"phi","mcpT","one"}, {"phi","mcpT","wpE"});

    for (idx_t p = 0; p < pts.num_particles; ++p) {
        const double wp = pts.dp("wpE", p);

        double Ep = pts.dp("mcpT", p);
        if (wp > WP_MIN) Ep /= wp;
        else             Ep = E_old[p];

        pts.dp("mcpT", p) = std::max(0.0, Ep);
        pts.dp("phi",  p) = std::clamp(pts.dp("phi", p), -1.0, 1.0);

        // mcp coerente (cp costante)
        const double mp = pts.dp("mp", p);
        pts.dp("mcp", p) = mp * CP_TI_SOLID;

        // PROIEZIONE: energia -> (T,fL)
        project_enthalpy_solid_liquid(pts, p);
    }

    // 7) BC su superficie (aggiunge energia e riproietta solo sulle particelle colpite)
    apply_plasma_bc_to_surface_points(
        pts, dt, cx, cy, R, grid.hx(), time, plasma_on_time
    );

    // 8) hard constraint globale: psi1=0 sempre, rho metal
    for (idx_t p = 0; p < pts.num_particles; ++p) {
        pts.dp("psi1", p) = 0.0;
        pts.dp("rho",  p) = RHO_TI_SOLID;
        pts.dp("T",    p) = std::clamp(pts.dp("T", p), T_AMB, T_MAX_CLAMP);
    }
}

/* ----------------------------- MAIN ----------------------------- */

int main(int argc, char **argv) {
    MPI_Init(&argc, &argv);

    int rank = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    if (rank == 0) {
        std::cout << "MPM + Phase-field: disco 2D diametro 50um, plasma 15000K (NO gas, solid/liquid con T_M)\n";
        std::cout << "Lx=Ly=200e-6, Nx=Ny=201, Np=50000\n";
        std::cout << "THICKNESS=" << THICKNESS << " m\n";
    }

    // Griglia
    grid_t grid(MPI_COMM_WORLD);
    double Lx = 200e-6, Ly = 200e-6;
    idx_t Nx = 201, Ny = 201;
    grid.set_sizes(Ny, Nx, Lx/(Nx-1), Ly/(Ny-1));

    if (rank == 0) {
        std::cout << "Grid " << Ny << "x" << Nx
                  << " hx=" << grid.hx() << " hy=" << grid.hy()
                  << " nn=" << grid.num_global_nodes() << "\n";
    }

    // Proprietà particelle (psi1 presente ma forzata a 0)
    std::vector<std::string> dprops = {
        "phi","T","psi1","psi2","psi3","rho","w","in_jet",
        "Vp","mp","mcp","mcpT",
        "one","wpE"
    };

    // Disco
    double cx = 100e-6;
    double cy = 0.5 * Ly;
    double R  = 25e-6;

    particles_t pts(50000, {}, dprops, grid);
    init_particle_disk_50um(pts, cx, cy, R);

    // Variabili di griglia
    std::map<std::string, dist_vec> gvars;
    for (auto &s : dprops)
        gvars[s] = dist_vec(grid.num_global_nodes(), 0.0);

    // variabili extra su griglia
    gvars["T"]   = dist_vec(grid.num_global_nodes(), T_AMB);
    gvars["one"] = dist_vec(grid.num_global_nodes(), 0.0);

    // dt diffusivo Ti
    double alpha = K_TI_SOLID / (RHO_TI_SOLID * CP_TI_SOLID);
    double hmin  = std::min(grid.hx(), grid.hy());
    double dt    = CFL_SAFE * hmin * hmin / alpha;

    if (rank == 0)
        std::cout << "dt = " << dt << " s\n";

    double t_end          = 2.0e-3;   // 2 ms
    double plasma_on_time = 1.0e-3;   // plasma ON 1 ms

    int num_steps = (int)std::ceil(t_end / dt);
    double time = 0.0;

    int out_every = 50;

    for (int step = 0; step < num_steps; ++step) {
        if (rank == 0 && (step % out_every == 0))
            std::cout << "Step " << step << "/" << num_steps << " t=" << time << " s\n";

        do_one_mpm_timestep(grid, pts, gvars, dt, time, cx, cy, R, plasma_on_time);

        if (rank == 0 && (step % out_every == 0)) {
            double Tmax = -1e300, Tmin = 1e300;
            int n_s=0, n_l=0;

            for (idx_t p = 0; p < pts.num_particles; ++p) {
                double Tp = pts.dp("T", p);
                Tmin = std::min(Tmin, Tp);
                Tmax = std::max(Tmax, Tp);

                if (pts.dp("psi3", p) > 0.5) n_l++;
                else                         n_s++;
            }

            std::cout << "  [particles] Tmin=" << Tmin << " K, Tmax=" << Tmax
                      << " K | solid=" << n_s << " liquid=" << n_l << "\n";

            const auto &mcp_g = gvars["mcp"];
            double mcp_min = 1e300, mcp_max = -1e300;
            int mcp_n = 0;
            for (double v : mcp_g) {
                if (v > EPS_MCP) { mcp_min = std::min(mcp_min, v); mcp_max = std::max(mcp_max, v); mcp_n++; }
            }
            std::cout << "  [grid] mcp valid nodes=" << mcp_n
                      << " mcp_min=" << (mcp_n?mcp_min:0.0)
                      << " mcp_max=" << (mcp_n?mcp_max:0.0) << "\n";
        }

        if (rank == 0 && (step % out_every == 0)) {
            write_vtk_vts(grid, gvars, step);
            write_particles_vtp(pts, step);
        }

        time += dt;
    }

    MPI_Finalize();
    return 0;
}

