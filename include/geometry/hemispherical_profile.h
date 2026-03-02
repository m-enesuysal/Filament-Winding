// =============================================================================
// hemispherical_profile.h — Hemispherical Dome Meridyen Profili
// =============================================================================
// Karar-10: HemisphericalProfile : IMeridianProfile
// Kapali-form analitik hesap. kappa_m = 1/R_eq (sabit).
// Phase-1a MATLAB referansi: hemispherical_dome_profile.m
// Detayli implementasyon: Phase-1b S3
// =============================================================================

#ifndef FILAMENT_GEOMETRY_HEMISPHERICAL_PROFILE_H
#define FILAMENT_GEOMETRY_HEMISPHERICAL_PROFILE_H

#include "geometry/i_meridian_profile.h"

namespace filament {
namespace geometry {

class HemisphericalProfile : public IMeridianProfile {
public:
    HemisphericalProfile() = default;

    MeridianLookupTable generateProfile(
        double R_eq, double r0, std::size_t N_points) const override;

    DomeType domeType() const override { return DomeType::Hemispherical; }
    const char* name() const override { return "Hemispherical"; }
};

} // namespace geometry
} // namespace filament

#endif // FILAMENT_GEOMETRY_HEMISPHERICAL_PROFILE_H
