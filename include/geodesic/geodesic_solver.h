// =============================================================================
// geodesic_solver.h — Geodesic Sarim Yolu ODE Cozucu
// =============================================================================
// Phase-2a S4: Tek devre geodesic yol hesabi
//
// Cekirdek ODE (Eq. 5.5):
//   dphi/ds = R_E / (rho(s) * sqrt(rho(s)^2 - R_E^2))
//
// Singularite yonetimi:
//   - rho > R_E + epsilon bolgesi: Dormand-Prince RK45 numerik integrasyon
//   - Son epsilon adimi: birlesik 1. + 2. derece analitik kuyruk
//   - 3 rejim: saf 1. derece, saf 2. derece (isotensoid), genel birlesik
//
// Tam devre: dome1_out -> cyl_fwd -> dome2_in -> dome2_out -> cyl_ret -> dome1_in
// delta_phi = 4 * phi_dome + 2 * phi_cyl
//
// Referans: docs/phase2a_geodesic_math.md Bolum 5, 6
// =============================================================================

#ifndef FILAMENT_GEODESIC_GEODESIC_SOLVER_H
#define FILAMENT_GEODESIC_GEODESIC_SOLVER_H

#include "geometry/filament_types.h"
#include <optional>
#include <vector>

namespace filament {
namespace geodesic {

// ---------------------------------------------------------------------------
// Cikti: tam devre geodesic yolu
// ---------------------------------------------------------------------------
struct GeodesicPath {
    std::vector<double> s;      // Geodesic yay uzunlugu [mm]
    std::vector<double> rho;    // Yaricap [mm]
    std::vector<double> x;      // Eksenel konum [mm]
    std::vector<double> phi;    // Azimuthal aci (kumulatif) [rad]
    std::vector<double> alpha;  // Sarim acisi [rad]

    double delta_phi;           // Tam devre aci ilerlemesi [rad]
    double phi_dome;            // Tek dome katkisi [rad]
    double phi_cyl;             // Tek silindir katkisi [rad]
    double phi_tail;            // Analitik kuyruk katkisi [rad]
    double phi_numerical;       // Numerik integrasyon katkisi [rad]
    double alpha_eq;            // Ekvator sarim acisi [rad]
};

// ---------------------------------------------------------------------------
// Girdi parametreleri
// ---------------------------------------------------------------------------
struct GeodesicParams {
    geometry::DomeType dome_type;
    double R_eq;            // Ekvator yaricapi [mm]
    double r0;              // Polar aciklik yaricapi [mm]
    double L_cyl;           // Silindir uzunlugu [mm]
    double R_E;             // Clairaut sabiti [mm] = r0 + BW_eff/2
    double k = 1.0;         // Elipsoidal en-boy orani (sadece Ellipsoidal)
    double epsilon = 1e-3;  // Singularite tamponu [mm] (Eq. 6.7)
    std::size_t N_dome_points = 500;  // Dome profili interpolasyon noktasi
};

// ---------------------------------------------------------------------------
// GeodesicSolver — statik cozucu
// ---------------------------------------------------------------------------
class GeodesicSolver {
public:
    // Tam devre geodesic yolu hesapla
    // Basarisizlik durumunda nullopt doner (ornegin R_E >= R_eq)
    static std::optional<GeodesicPath> solve(const GeodesicParams& params);
};

// ---------------------------------------------------------------------------
// resample — PCHIP ile esit aralikli yeniden ornekleme
// ---------------------------------------------------------------------------
// Girdi path'in s, rho, x, phi, alpha vektorlerini delta_s adimlarla
// yeniden ornekler. Fritsch-Carlson egim kestirimi kullanir.
// Skaler alanlar (delta_phi, phi_dome, ...) kopyalanir.
GeodesicPath resample(const GeodesicPath& path, double delta_s);

} // namespace geodesic
} // namespace filament

#endif // FILAMENT_GEODESIC_GEODESIC_SOLVER_H
