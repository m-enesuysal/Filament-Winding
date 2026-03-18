// =============================================================================
// mandrel_geometry.h — Tam Mandrel Yuzey Temsili
// =============================================================================
// Karar-10: MandrelGeometry ust katman — silindir + iki dome birlesimi
// Karar-21: Immutable + factory pattern
// Karar-2:  Simetrik mandrel (dome1 = dome2)
//
// Global yay uzunlugu koordinat sistemi:
//   s_global in [0, s_dome + L_cyl + s_dome]
//   Bolge 1: Dome-1       [0, s_dome)
//   Bolge 2: Silindir     [s_dome, s_dome + L_cyl)
//   Bolge 3: Dome-2       [s_dome + L_cyl, 2*s_dome + L_cyl]
//
// Silindirik bolge: analitik tanim (rho = R_eq, kappa_m = 0, beta = 0)
// Detayli implementasyon: Phase-1b S6
// =============================================================================

#ifndef FILAMENT_GEOMETRY_MANDREL_GEOMETRY_H
#define FILAMENT_GEOMETRY_MANDREL_GEOMETRY_H

#include "geometry/meridian_lookup_table.h"
#include "geometry/filament_types.h"
#include <memory>

namespace filament {
namespace geometry {

class IMeridianProfile;

class MandrelGeometry {
public:
    // Factory constructor (Karar-21: immutable, tum hesap constructor'da)
    MandrelGeometry(
        DomeType dome_type,
        double R_eq,
        double r0,
        double L_cyl,
        double k = 1.0,
        std::size_t N_points = 500
    );

    // Sorgu: global yay uzunlugundan nokta
    std::optional<MeridianPoint> point(double s_global) const;

    // Bolge tespiti (Karar-11 Katman 3: +-1e-6 mm tolerans bandi)
    bool isOnDome1(double s_global) const;
    bool isOnCylinder(double s_global) const;
    bool isOnDome2(double s_global) const;

    // Clairaut sabiti hesabi (Karar-7: alfa konvansiyonu)
    double windingAngleToClairaut(double alpha, double rho) const;

    // Meta-veri erisimi
    double R_eq() const;
    double r0() const;
    double L_cyl() const;
    double totalLength() const;
    DomeType domeType() const;
    const ProfileMetadata& domeMetadata() const;
    const MeridianLookupTable& domeTable() const { return dome_table_; }

private:
    MeridianLookupTable dome_table_;
    double R_eq_;
    double r0_;
    double L_cyl_;
    double s_dome_;
    double s_total_;
    DomeType dome_type_;
};

} // namespace geometry
} // namespace filament

#endif // FILAMENT_GEOMETRY_MANDREL_GEOMETRY_H
