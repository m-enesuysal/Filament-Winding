// =============================================================================
// hemispherical_profile.cpp — Placeholder (Phase-1b S3'te implement edilecek)
// =============================================================================

#include "geometry/hemispherical_profile.h"
#include <stdexcept>

namespace filament {
namespace geometry {

MeridianLookupTable HemisphericalProfile::generateProfile(
    double /*R_eq*/, double /*r0*/, std::size_t /*N_points*/) const
{
    throw std::runtime_error(
        "HemisphericalProfile::generateProfile() henuz implement edilmedi. "
        "Phase-1b S3 session'inda yapilacak.");
}

} // namespace geometry
} // namespace filament
