// =============================================================================
// meridian_lookup_table.cpp — Placeholder (Phase-1b S2'de implement edilecek)
// =============================================================================

#include "geometry/meridian_lookup_table.h"

namespace filament {
namespace geometry {

void MeridianLookupTable::build(
    const std::vector<double>& s,
    const std::vector<double>& rho,
    const std::vector<double>& x_local,
    const std::vector<double>& drho_ds,
    const std::vector<double>& dx_ds,
    const std::vector<double>& kappa_m,
    const ProfileMetadata& metadata)
{
    // S2'de implement edilecek: veri kopyalama + spline katsayi hesabi
    s_ = s;
    rho_ = rho;
    x_local_ = x_local;
    drho_ds_ = drho_ds;
    dx_ds_ = dx_ds;
    kappa_m_ = kappa_m;
    metadata_ = metadata;
    valid_ = !s_.empty();
}

std::optional<MeridianPoint> MeridianLookupTable::query(double /*s*/) const
{
    // S2'de implement edilecek: binary search + kubik spline interpolasyon
    return std::nullopt;
}

const ProfileMetadata& MeridianLookupTable::metadata() const { return metadata_; }
std::size_t MeridianLookupTable::size() const { return s_.size(); }
bool MeridianLookupTable::isValid() const { return valid_; }
const std::vector<double>& MeridianLookupTable::rawS() const { return s_; }
const std::vector<double>& MeridianLookupTable::rawRho() const { return rho_; }
const std::vector<double>& MeridianLookupTable::rawX() const { return x_local_; }

} // namespace geometry
} // namespace filament
