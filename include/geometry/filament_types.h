// =============================================================================
// filament_types.h — Ortak Tip Tanimlari ve Sabitler
// =============================================================================
// Karar-8:  Birim sistemi konvansiyonu (mm, radyan dahili)
// Karar-11: Numerik tolerans ve yakinlasma kriterleri
// Karar-15: Singularite katalogu (S-GEO-03 k_min)
// =============================================================================

#ifndef FILAMENT_GEOMETRY_TYPES_H
#define FILAMENT_GEOMETRY_TYPES_H

#include <cstddef>
#include <cmath>

namespace filament {
namespace geometry {

// ---------------------------------------------------------------------------
// Dome tipi enum (Karar-5)
// ---------------------------------------------------------------------------
enum class DomeType {
    Hemispherical,  // Yaricap = R_eq (elipsoidal k=1 ozel hali)
    Ellipsoidal,    // Aspect ratio k = h_dome / R_eq
    Isotensoid      // Netting theory ODE'sinden uretilir
};

// ---------------------------------------------------------------------------
// Numerik Toleranslar (Karar-11)
// ---------------------------------------------------------------------------
namespace tolerances {

    // Katman 1 — ODE Integrasyon Toleransi
    constexpr double ODE_REL_TOL = 1e-8;
    constexpr double ODE_ABS_TOL = 1e-10;

    // Katman 2 — Interpolasyon Hatasi Toleransi
    constexpr double POSITION_ABS_TOL   = 1e-4;  // mm
    constexpr double DERIVATIVE_ABS_TOL = 1e-6;  // boyutsuz
    constexpr double CURVATURE_REL_TOL  = 1e-4;  // bagil (%0.01)

    // Katman 3 — Geometrik Sorgu Toleransi
    // Silindir-dome sinirinda tolerans bandi.
    // Band icindeki noktalar dome tarafina atanir (Karar-11).
    constexpr double JUNCTION_BAND_TOL = 1e-6;   // mm

} // namespace tolerances

// ---------------------------------------------------------------------------
// Geometrik Sinir Degerleri (Karar-15)
// ---------------------------------------------------------------------------
namespace limits {

    // S-GEO-03: Elipsoidal dome minimum aspect ratio
    constexpr double ELLIPSOIDAL_K_MIN = 0.15;

} // namespace limits

// ---------------------------------------------------------------------------
// Matematiksel Sabitler
// ---------------------------------------------------------------------------
namespace constants {

    constexpr double PI     = 3.14159265358979323846;
    constexpr double TWO_PI = 2.0 * PI;

} // namespace constants

} // namespace geometry
} // namespace filament

#endif // FILAMENT_GEOMETRY_TYPES_H
