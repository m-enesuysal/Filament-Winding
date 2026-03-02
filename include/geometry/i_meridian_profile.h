// =============================================================================
// i_meridian_profile.h — Meridyen Profil Soyut Arayuzu
// =============================================================================
// Karar-10: IMeridianProfile soyut arayuzu
// Karar-3:  Birlesik parametrik meridyen temsili
//
// Tum dome tipleri (isotensoid, elipsoidal, hemispherical) bu arayuzu
// implement eder. Phase-2 geodesic solver dome tipinden bagimsiz olarak
// bu arayuz uzerinden calisacaktir.
// =============================================================================

#ifndef FILAMENT_GEOMETRY_I_MERIDIAN_PROFILE_H
#define FILAMENT_GEOMETRY_I_MERIDIAN_PROFILE_H

#include "geometry/meridian_lookup_table.h"
#include "geometry/filament_types.h"
#include <memory>

namespace filament {
namespace geometry {

class IMeridianProfile {
public:
    virtual ~IMeridianProfile() = default;

    // Profil uretimi — hesaplanan tum buyuklukleri tabloya yazar
    // Ust katman (Karar-17): gecersiz parametrede exception firlatir
    virtual MeridianLookupTable generateProfile(
        double R_eq,
        double r0,
        std::size_t N_points
    ) const = 0;

    // Dome tipi sorgusu
    virtual DomeType domeType() const = 0;

    // Insan okunur isim
    virtual const char* name() const = 0;

protected:
    IMeridianProfile() = default;
    IMeridianProfile(const IMeridianProfile&) = default;
    IMeridianProfile& operator=(const IMeridianProfile&) = default;
};

// Factory fonksiyonu — dome tipine gore implementasyonu olusturur
std::unique_ptr<IMeridianProfile> createProfile(
    DomeType type,
    double k = 1.0
);

} // namespace geometry
} // namespace filament

#endif // FILAMENT_GEOMETRY_I_MERIDIAN_PROFILE_H
