// =============================================================================
// ellipsoidal_profile.cpp — Elipsoidal Dome Meridyen Profili
// =============================================================================
// Karar-10: EllipsoidalProfile : IMeridianProfile
// S-GEO-03: k_min = 0.15 alt sinir
// S-GEO-04: k=1 → hemispherical ozel hali (ortusme dogrulanmali)
//
// Parametrizasyon:
//   Yari-eksenler: a = R_eq (yaricap), b = k * R_eq (aksiyel)
//   theta: 0 (ekvator) → theta_p = acos(r0/R_eq) (polar aciklik)
//   rho(theta) = R_eq * cos(theta)
//   x(theta)   = k * R_eq * sin(theta)
//   f(theta)   = sqrt(sin^2(theta) + k^2 * cos^2(theta))
//   ds/dtheta  = R_eq * f(theta)
//   drho/ds    = -sin(theta) / f(theta)
//   dx/ds      = k * cos(theta) / f(theta)
//   kappa_m    = k / (R_eq * f^3(theta))
//
// Yay uzunlugu: Eliptik integral — composite Simpson sayisal integrasyon
// Phase-1a MATLAB referansi: ellipsoidal_dome_profile.m
// =============================================================================

#include "geometry/ellipsoidal_profile.h"
#include <cmath>
#include <stdexcept>
#include <string>
#include <vector>

namespace filament {
namespace geometry {

// ---------------------------------------------------------------------------
// Constructor — girdi dogrulama S1 placeholder'dan korunuyor
// ---------------------------------------------------------------------------
EllipsoidalProfile::EllipsoidalProfile(double k)
    : k_(k)
{
    if (k <= 0.0) {
        throw std::invalid_argument(
            "Elipsoidal dome aspect ratio k pozitif olmali. k = "
            + std::to_string(k));
    }
    if (k < limits::ELLIPSOIDAL_K_MIN) {
        throw std::invalid_argument(
            "Elipsoidal dome aspect ratio k = " + std::to_string(k)
            + " alt sinirin altinda (k_min = "
            + std::to_string(limits::ELLIPSOIDAL_K_MIN)
            + "). S-GEO-03: dome dejenere olur.");
    }
}

// ---------------------------------------------------------------------------
// generateProfile — Elipsoidal dome profil uretimi
// ---------------------------------------------------------------------------
MeridianLookupTable EllipsoidalProfile::generateProfile(
    double R_eq, double r0, std::size_t N_points) const
{
    // -----------------------------------------------------------------------
    // Girdi dogrulama — ust katman (Karar-17): exception firlatir
    // -----------------------------------------------------------------------
    if (R_eq <= 0.0) {
        throw std::invalid_argument(
            "EllipsoidalProfile: R_eq pozitif olmali. R_eq = "
            + std::to_string(R_eq));
    }
    if (r0 <= 0.0) {
        throw std::invalid_argument(
            "EllipsoidalProfile: r0 pozitif olmali. r0 = "
            + std::to_string(r0));
    }
    if (r0 >= R_eq) {
        throw std::invalid_argument(
            "EllipsoidalProfile: r0 < R_eq olmali. r0 = "
            + std::to_string(r0) + ", R_eq = " + std::to_string(R_eq));
    }
    if (N_points < 2) {
        throw std::invalid_argument(
            "EllipsoidalProfile: N_points en az 2 olmali. N_points = "
            + std::to_string(N_points));
    }

    // -----------------------------------------------------------------------
    // Elipsoidal profil hesabi
    // -----------------------------------------------------------------------
    const double theta_p = std::acos(r0 / R_eq);
    const double k2 = k_ * k_;  // k^2

    // f(theta) = sqrt(sin^2(theta) + k^2 * cos^2(theta))
    auto f = [k2](double theta) -> double {
        const double st = std::sin(theta);
        const double ct = std::cos(theta);
        return std::sqrt(st * st + k2 * ct * ct);
    };

    // --- Uniform theta dagitimi ---
    std::vector<double> theta(N_points);
    for (std::size_t i = 0; i < N_points; ++i) {
        theta[i] = theta_p
            * static_cast<double>(i)
            / static_cast<double>(N_points - 1);
    }

    // --- rho, x, f degerleri ---
    std::vector<double> rho(N_points);
    std::vector<double> x_local(N_points);
    std::vector<double> f_vals(N_points);

    for (std::size_t i = 0; i < N_points; ++i) {
        rho[i]     = R_eq * std::cos(theta[i]);
        x_local[i] = k_ * R_eq * std::sin(theta[i]);
        f_vals[i]  = f(theta[i]);
    }

    // --- Yay uzunlugu: kumulatif composite Simpson integrasyon ---
    // Her aralik [theta_i, theta_{i+1}] icin 4-noktali composite Simpson
    // ds = R_eq * integral(f(t), theta_i, theta_{i+1})
    constexpr int N_SUB = 4;  // Simpson alt-aralik sayisi (cift olmali)
    std::vector<double> s(N_points);
    s[0] = 0.0;

    for (std::size_t i = 0; i + 1 < N_points; ++i) {
        const double ta = theta[i];
        const double tb = theta[i + 1];
        const double h  = (tb - ta) / N_SUB;

        // Composite Simpson: (h/3) * [f(t0) + 4*f(t1) + 2*f(t2) + 4*f(t3) + f(t4)]
        double sum = f(ta) + f(tb);
        for (int j = 1; j < N_SUB; j += 2) {
            sum += 4.0 * f(ta + j * h);
        }
        for (int j = 2; j < N_SUB; j += 2) {
            sum += 2.0 * f(ta + j * h);
        }
        s[i + 1] = s[i] + R_eq * sum * h / 3.0;
    }

    // --- Turevler ve egrilik ---
    std::vector<double> drho_ds(N_points);
    std::vector<double> dx_ds(N_points);
    std::vector<double> kappa_m(N_points);

    for (std::size_t i = 0; i < N_points; ++i) {
        const double fi  = f_vals[i];
        const double fi3 = fi * fi * fi;

        drho_ds[i] = -std::sin(theta[i]) / fi;
        dx_ds[i]   =  k_ * std::cos(theta[i]) / fi;
        kappa_m[i] =  k_ / (R_eq * fi3);
    }

    // -----------------------------------------------------------------------
    // Meta-veri hesabi
    // -----------------------------------------------------------------------
    const double s_total = s.back();
    const double h_dome  = k_ * R_eq * std::sin(theta_p);

    // Yuzey alani: A = 2*pi * integral(rho * ds) — trapezoidal kural
    double A_dome = 0.0;
    for (std::size_t i = 0; i + 1 < N_points; ++i) {
        A_dome += 0.5 * (rho[i] + rho[i + 1]) * (s[i + 1] - s[i]);
    }
    A_dome *= 2.0 * constants::PI;

    ProfileMetadata meta;
    meta.R_eq      = R_eq;
    meta.r0        = r0;
    meta.s_total   = s_total;
    meta.h_dome    = h_dome;
    meta.A_dome    = A_dome;
    meta.kappa_eq  = kappa_m.front();     // = 1/(R_eq * k^2)
    meta.kappa_pol = kappa_m.back();      // = k/(R_eq * f^3(theta_p))
    meta.alpha_w   = std::asin(r0 / R_eq);
    meta.aspect_r  = k_;

    // -----------------------------------------------------------------------
    // Tablo olustur
    // -----------------------------------------------------------------------
    MeridianLookupTable table;
    table.build(s, rho, x_local, drho_ds, dx_ds, kappa_m, meta);

    return table;
}

} // namespace geometry
} // namespace filament
