// =============================================================================
// hemispherical_profile.cpp — Hemispherical Dome Meridyen Profili
// =============================================================================
// Karar-10: HemisphericalProfile : IMeridianProfile
// Karar-5:  Girdi dogrulamasi (r0 < R_eq, r0 > 0)
// Karar-17: Ust katman — gecersiz parametrelerde exception firlatir
//
// Kapali-form analitik hesap:
//   Parametrizasyon: theta = 0 (ekvator) → theta_p (polar aciklik)
//   theta_p = acos(r0 / R_eq)
//   rho(theta) = R_eq * cos(theta)        — yaricap
//   x(theta)   = R_eq * sin(theta)        — aksiyel konum
//   s(theta)   = R_eq * theta             — yay uzunlugu
//   drho/ds    = -sin(theta)              — yaricap turevi
//   dx/ds      = cos(theta)               — aksiyel turev
//   kappa_m    = 1/R_eq                   — sabit meridyen egriligi
//
// Phase-1a MATLAB referansi: hemispherical_dome_profile.m
// =============================================================================

#include "geometry/hemispherical_profile.h"
#include <cmath>
#include <stdexcept>
#include <string>

namespace filament {
namespace geometry {

MeridianLookupTable HemisphericalProfile::generateProfile(
    double R_eq, double r0, std::size_t N_points) const
{
    // -----------------------------------------------------------------------
    // Girdi dogrulama — ust katman (Karar-17): exception firlatir
    // -----------------------------------------------------------------------
    if (R_eq <= 0.0) {
        throw std::invalid_argument(
            "HemisphericalProfile: R_eq pozitif olmali. R_eq = "
            + std::to_string(R_eq));
    }
    if (r0 <= 0.0) {
        throw std::invalid_argument(
            "HemisphericalProfile: r0 pozitif olmali. r0 = "
            + std::to_string(r0));
    }
    if (r0 >= R_eq) {
        throw std::invalid_argument(
            "HemisphericalProfile: r0 < R_eq olmali. r0 = "
            + std::to_string(r0) + ", R_eq = " + std::to_string(R_eq));
    }
    if (N_points < 2) {
        throw std::invalid_argument(
            "HemisphericalProfile: N_points en az 2 olmali. N_points = "
            + std::to_string(N_points));
    }

    // -----------------------------------------------------------------------
    // Analitik profil hesabi
    // -----------------------------------------------------------------------
    const double theta_p = std::acos(r0 / R_eq);
    const double s_total = R_eq * theta_p;
    const double kappa   = 1.0 / R_eq;  // sabit meridyen egriligi

    std::vector<double> s(N_points);
    std::vector<double> rho(N_points);
    std::vector<double> x_local(N_points);
    std::vector<double> drho_ds(N_points);
    std::vector<double> dx_ds(N_points);
    std::vector<double> kappa_m(N_points);

    for (std::size_t i = 0; i < N_points; ++i) {
        // Uniform theta dagitimi: theta_i = theta_p * i / (N-1)
        const double theta = theta_p
            * static_cast<double>(i)
            / static_cast<double>(N_points - 1);

        s[i]       = R_eq * theta;
        rho[i]     = R_eq * std::cos(theta);
        x_local[i] = R_eq * std::sin(theta);
        drho_ds[i] = -std::sin(theta);
        dx_ds[i]   =  std::cos(theta);
        kappa_m[i] = kappa;
    }

    // -----------------------------------------------------------------------
    // Meta-veri hesabi
    // -----------------------------------------------------------------------
    const double h_dome = R_eq * std::sin(theta_p);  // = sqrt(R_eq^2 - r0^2)

    ProfileMetadata meta;
    meta.R_eq      = R_eq;
    meta.r0        = r0;
    meta.s_total   = s_total;
    meta.h_dome    = h_dome;
    meta.A_dome    = 2.0 * constants::PI * R_eq * h_dome;  // kuresel kalotte
    meta.kappa_eq  = kappa;   // ekvator meridyen egriligi
    meta.kappa_pol = kappa;   // polar aciklik meridyen egriligi (kure: sabit)
    meta.alpha_w   = std::asin(r0 / R_eq);  // Clairaut ekvator sarma acisi
    meta.aspect_r  = 1.0;    // hemispherical dome k = 1

    // -----------------------------------------------------------------------
    // Tablo olustur
    // -----------------------------------------------------------------------
    MeridianLookupTable table;
    table.build(s, rho, x_local, drho_ds, dx_ds, kappa_m, meta);

    return table;
}

} // namespace geometry
} // namespace filament
