// =============================================================================
// mandrel_geometry.cpp — Placeholder (Phase-1b S6'da implement edilecek)
// =============================================================================

#include "geometry/mandrel_geometry.h"
#include "geometry/i_meridian_profile.h"
#include "geometry/hemispherical_profile.h"
#include "geometry/ellipsoidal_profile.h"
#include "geometry/isotensoid_profile.h"
#include <stdexcept>

namespace filament {
namespace geometry {

// Factory fonksiyonu (Karar-10)
std::unique_ptr<IMeridianProfile> createProfile(DomeType type, double k)
{
    switch (type) {
        case DomeType::Hemispherical:
            return std::make_unique<HemisphericalProfile>();
        case DomeType::Ellipsoidal:
            return std::make_unique<EllipsoidalProfile>(k);
        case DomeType::Isotensoid:
            return std::make_unique<IsotensoidProfile>();
        default:
            throw std::invalid_argument("Bilinmeyen dome tipi.");
    }
}

MandrelGeometry::MandrelGeometry(
    DomeType dome_type, double R_eq, double r0,
    double L_cyl, double k, std::size_t N_points)
    : R_eq_(R_eq), r0_(r0), L_cyl_(L_cyl)
    , s_dome_(0.0), s_total_(0.0), dome_type_(dome_type)
{
    (void)k; (void)N_points;
    throw std::runtime_error(
        "MandrelGeometry constructor henuz implement edilmedi. "
        "Phase-1b S6 session'inda yapilacak.");
}

std::optional<MeridianPoint> MandrelGeometry::point(double /*s_global*/) const
{ return std::nullopt; }

bool MandrelGeometry::isOnDome1(double /*s*/) const { return false; }
bool MandrelGeometry::isOnCylinder(double /*s*/) const { return false; }
bool MandrelGeometry::isOnDome2(double /*s*/) const { return false; }

double MandrelGeometry::windingAngleToClairaut(double alpha, double rho) const
{ return rho * std::sin(alpha); }

double MandrelGeometry::R_eq() const { return R_eq_; }
double MandrelGeometry::r0() const { return r0_; }
double MandrelGeometry::L_cyl() const { return L_cyl_; }
double MandrelGeometry::totalLength() const { return s_total_; }
DomeType MandrelGeometry::domeType() const { return dome_type_; }
const ProfileMetadata& MandrelGeometry::domeMetadata() const
{ return dome_table_.metadata(); }

} // namespace geometry
} // namespace filament
