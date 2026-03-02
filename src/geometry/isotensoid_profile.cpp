// =============================================================================
// isotensoid_profile.cpp — Placeholder (Phase-1b S5'te implement edilecek)
// =============================================================================

#include "geometry/isotensoid_profile.h"
#include <stdexcept>

namespace filament {
namespace geometry {

MeridianLookupTable IsotensoidProfile::generateProfile(
    double /*R_eq*/, double /*r0*/, std::size_t /*N_points*/) const
{
    throw std::runtime_error(
        "IsotensoidProfile::generateProfile() henuz implement edilmedi. "
        "Phase-1b S5 session'inda yapilacak.");
}

} // namespace geometry
} // namespace filament
