// =============================================================================
// ellipsoidal_profile.cpp — Placeholder (Phase-1b S4'te implement edilecek)
// =============================================================================

#include "geometry/ellipsoidal_profile.h"
#include <stdexcept>

namespace filament {
namespace geometry {

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

MeridianLookupTable EllipsoidalProfile::generateProfile(
    double /*R_eq*/, double /*r0*/, std::size_t /*N_points*/) const
{
    throw std::runtime_error(
        "EllipsoidalProfile::generateProfile() henuz implement edilmedi. "
        "Phase-1b S4 session'inda yapilacak.");
}

} // namespace geometry
} // namespace filament
