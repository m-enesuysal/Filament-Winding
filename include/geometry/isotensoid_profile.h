// =============================================================================
// isotensoid_profile.h — Izotensoid Dome Meridyen Profili
// =============================================================================
// Karar-10: IsotensoidProfile : IMeridianProfile
// Karar-6:  Boost.Odeint Dormand-Prince (RK45) adaptif cozucu
// S-GEO-01: Polar aciklik yakininda stiff davranis
// Phase-1a MATLAB referansi: isotensoid_dome_profile.m
// Detayli implementasyon: Phase-1b S5
// =============================================================================

#ifndef FILAMENT_GEOMETRY_ISOTENSOID_PROFILE_H
#define FILAMENT_GEOMETRY_ISOTENSOID_PROFILE_H

#include "geometry/i_meridian_profile.h"

namespace filament {
namespace geometry {

class IsotensoidProfile : public IMeridianProfile {
public:
    IsotensoidProfile() = default;

    MeridianLookupTable generateProfile(
        double R_eq, double r0, std::size_t N_points) const override;

    DomeType domeType() const override { return DomeType::Isotensoid; }
    const char* name() const override { return "Isotensoid"; }
};

} // namespace geometry
} // namespace filament

#endif // FILAMENT_GEOMETRY_ISOTENSOID_PROFILE_H
