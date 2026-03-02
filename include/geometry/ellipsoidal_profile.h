// =============================================================================
// ellipsoidal_profile.h — Elipsoidal Dome Meridyen Profili
// =============================================================================
// Karar-10: EllipsoidalProfile : IMeridianProfile
// S-GEO-03: k_min = 0.15 alt sinir
// S-GEO-04: k=1 → hemispherical ortusme
// Phase-1a MATLAB referansi: ellipsoidal_dome_profile.m
// Detayli implementasyon: Phase-1b S4
// =============================================================================

#ifndef FILAMENT_GEOMETRY_ELLIPSOIDAL_PROFILE_H
#define FILAMENT_GEOMETRY_ELLIPSOIDAL_PROFILE_H

#include "geometry/i_meridian_profile.h"

namespace filament {
namespace geometry {

class EllipsoidalProfile : public IMeridianProfile {
public:
    explicit EllipsoidalProfile(double k);

    MeridianLookupTable generateProfile(
        double R_eq, double r0, std::size_t N_points) const override;

    DomeType domeType() const override { return DomeType::Ellipsoidal; }
    const char* name() const override { return "Ellipsoidal"; }
    double k() const { return k_; }

private:
    double k_;
};

} // namespace geometry
} // namespace filament

#endif // FILAMENT_GEOMETRY_ELLIPSOIDAL_PROFILE_H
