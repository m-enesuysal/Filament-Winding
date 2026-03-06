// =============================================================================
// mandrel_geometry.cpp — Tam Mandrel Yuzey Temsili
// =============================================================================
// Karar-10: MandrelGeometry ust katman — silindir + iki dome birlesimi
// Karar-21: Immutable + factory pattern (tum hesap constructor'da)
// Karar-2:  Simetrik mandrel (dome1 = dome2)
// Karar-11: Katman 3 — Junction band tolerance +-1e-6 mm, dome tarafina atama
//
// Global yay uzunlugu koordinat sistemi:
//   s_global in [0, 2*s_dome + L_cyl]
//   Bolge 1: Dome-1       [0, s_dome)         polar(s=0) → ekvator(s=s_dome)
//   Bolge 2: Silindir     [s_dome, s_dome+L)
//   Bolge 3: Dome-2       [s_dome+L, 2*s_dome+L]  ekvator → polar
//
// Dome-1 donusumu (ters parametrizasyon):
//   s_local = s_dome - s_global
//   rho = dome.rho (ayni)
//   x_local = h_dome - dome.x_local (ters)
//   drho/ds_global = -drho/ds_dome (s tersleme)
//   dx/ds_global = dx/ds_dome (cift olumsuzlama: s tersleme + x tersleme)
//   kappa_m_global = -kappa_m_dome (s tersleme)
//
// Dome-2 donusumu (ayni yon):
//   s_local = s_global - s_dome - L_cyl
//   rho = dome.rho, drho/ds = dome.drho/ds, dx/ds = dome.dx/ds
//   x_local = h_dome + L_cyl + dome.x_local
//   kappa_m = dome.kappa_m
//
// Silindirik bolge (analitik):
//   rho = R_eq, drho/ds = 0, dx/ds = 1, kappa_m = 0
//   x_local = h_dome + (s_global - s_dome)
//
// C1 sureklilik: Ekvator noktasinda rho, drho/ds, dx/ds uyumu — runtime assert
// =============================================================================

#include "geometry/mandrel_geometry.h"
#include "geometry/i_meridian_profile.h"
#include "geometry/hemispherical_profile.h"
#include "geometry/ellipsoidal_profile.h"
#include "geometry/isotensoid_profile.h"
#include <stdexcept>
#include <cmath>
#include <cassert>
#include <string>

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

// ===========================================================================
// Constructor (Karar-21: immutable, tum hesap burada)
// ===========================================================================
MandrelGeometry::MandrelGeometry(
    DomeType dome_type, double R_eq, double r0,
    double L_cyl, double k, std::size_t N_points)
    : R_eq_(R_eq), r0_(r0), L_cyl_(L_cyl)
    , s_dome_(0.0), s_total_(0.0), dome_type_(dome_type)
{
    // -------------------------------------------------------------------
    // Girdi dogrulama (Karar-5 + Karar-17)
    // -------------------------------------------------------------------
    if (R_eq <= 0.0) {
        throw std::invalid_argument(
            "MandrelGeometry: R_eq pozitif olmali. R_eq = "
            + std::to_string(R_eq));
    }
    if (r0 <= 0.0) {
        throw std::invalid_argument(
            "MandrelGeometry: r0 pozitif olmali. r0 = "
            + std::to_string(r0));
    }
    if (r0 >= R_eq) {
        throw std::invalid_argument(
            "MandrelGeometry: r0 < R_eq olmali. r0 = "
            + std::to_string(r0) + ", R_eq = " + std::to_string(R_eq));
    }
    if (L_cyl < 0.0) {
        throw std::invalid_argument(
            "MandrelGeometry: L_cyl negatif olamaz. L_cyl = "
            + std::to_string(L_cyl));
    }
    if (N_points < 2) {
        throw std::invalid_argument(
            "MandrelGeometry: N_points en az 2 olmali. N_points = "
            + std::to_string(N_points));
    }

    // -------------------------------------------------------------------
    // Dome profili olustur (factory pattern)
    // -------------------------------------------------------------------
    auto profile = createProfile(dome_type, k);
    dome_table_ = profile->generateProfile(R_eq, r0, N_points);

    // -------------------------------------------------------------------
    // Global koordinat sistemi
    // -------------------------------------------------------------------
    s_dome_  = dome_table_.metadata().s_total;
    s_total_ = 2.0 * s_dome_ + L_cyl_;

    // -------------------------------------------------------------------
    // C1 sureklilik runtime assertion
    // Dome ekvator noktasi (s_local=0) silindir ile uyumlu olmali:
    //   rho(0) = R_eq,  drho/ds(0) = 0,  dx/ds(0) = 1
    // -------------------------------------------------------------------
    auto eq_pt = dome_table_.query(0.0);
    assert(eq_pt.has_value() && "Dome ekvator noktasi sorgulanamadi");

    const double rho_err    = std::abs(eq_pt->rho - R_eq);
    const double drho_err   = std::abs(eq_pt->drho_ds);
    const double dxds_err   = std::abs(eq_pt->dx_ds - 1.0);

    assert(rho_err  < 1e-6 && "C1 sureklilik ihlali: rho(ekvator) != R_eq");
    assert(drho_err < 1e-6 && "C1 sureklilik ihlali: drho/ds(ekvator) != 0");
    assert(dxds_err < 1e-6 && "C1 sureklilik ihlali: dx/ds(ekvator) != 1");

    // Kullanilmayan degiskenleri sustur (Release build'de assert kaldirilir)
    (void)rho_err;
    (void)drho_err;
    (void)dxds_err;
}

// ===========================================================================
// point() — Global yay uzunlugundan MeridianPoint sorgusu
// ===========================================================================
std::optional<MeridianPoint> MandrelGeometry::point(double s_global) const
{
    // Aralik kontrolu
    if (s_global < 0.0 || s_global > s_total_) {
        return std::nullopt;
    }

    const double h_dome = dome_table_.metadata().h_dome;
    const double J      = tolerances::JUNCTION_BAND_TOL;

    // -------------------------------------------------------------------
    // Bolge 1: Dome-1 [0, s_dome + J]
    // Ters parametrizasyon: s_local = s_dome - s_global
    // -------------------------------------------------------------------
    if (s_global <= s_dome_ + J) {
        const double s_local = s_dome_ - s_global;
        // s_local sinirlarini koru
        const double s_clamped = (s_local < 0.0) ? 0.0
                               : (s_local > s_dome_) ? s_dome_
                               : s_local;

        auto dome_pt = dome_table_.query(s_clamped);
        if (!dome_pt.has_value()) return std::nullopt;

        MeridianPoint pt;
        pt.s        = s_global;
        pt.rho      = dome_pt->rho;
        pt.x_local  = h_dome - dome_pt->x_local;
        pt.drho_ds  = -dome_pt->drho_ds;     // s tersleme
        pt.dx_ds    =  dome_pt->dx_ds;        // cift olumsuzlama
        pt.kappa_m  = -dome_pt->kappa_m;      // s tersleme
        return pt;
    }

    // -------------------------------------------------------------------
    // Bolge 3: Dome-2 [s_dome + L_cyl - J, s_total]
    // Ayni yon: s_local = s_global - s_dome - L_cyl
    // -------------------------------------------------------------------
    if (s_global >= s_dome_ + L_cyl_ - J) {
        const double s_local = s_global - s_dome_ - L_cyl_;
        const double s_clamped = (s_local < 0.0) ? 0.0
                               : (s_local > s_dome_) ? s_dome_
                               : s_local;

        auto dome_pt = dome_table_.query(s_clamped);
        if (!dome_pt.has_value()) return std::nullopt;

        MeridianPoint pt;
        pt.s        = s_global;
        pt.rho      = dome_pt->rho;
        pt.x_local  = h_dome + L_cyl_ + dome_pt->x_local;
        pt.drho_ds  = dome_pt->drho_ds;
        pt.dx_ds    = dome_pt->dx_ds;
        pt.kappa_m  = dome_pt->kappa_m;
        return pt;
    }

    // -------------------------------------------------------------------
    // Bolge 2: Silindir (analitik)
    // rho = R_eq, drho/ds = 0, dx/ds = 1, kappa_m = 0
    // -------------------------------------------------------------------
    MeridianPoint pt;
    pt.s        = s_global;
    pt.rho      = R_eq_;
    pt.x_local  = h_dome + (s_global - s_dome_);
    pt.drho_ds  = 0.0;
    pt.dx_ds    = 1.0;
    pt.kappa_m  = 0.0;
    return pt;
}

// ===========================================================================
// Bolge tespiti (Karar-11 Katman 3: +-1e-6 mm tolerans bandi)
// Sinir noktalarinda dome tarafi oncelikli
// ===========================================================================
bool MandrelGeometry::isOnDome1(double s) const
{
    return s >= 0.0
        && s <= s_dome_ + tolerances::JUNCTION_BAND_TOL;
}

bool MandrelGeometry::isOnCylinder(double s) const
{
    return s > s_dome_ + tolerances::JUNCTION_BAND_TOL
        && s < s_dome_ + L_cyl_ - tolerances::JUNCTION_BAND_TOL;
}

bool MandrelGeometry::isOnDome2(double s) const
{
    return s >= s_dome_ + L_cyl_ - tolerances::JUNCTION_BAND_TOL
        && s <= s_total_;
}

// ===========================================================================
// Clairaut sabiti hesabi (Karar-7)
// ===========================================================================
double MandrelGeometry::windingAngleToClairaut(double alpha, double rho) const
{
    return rho * std::sin(alpha);
}

// ===========================================================================
// Meta-veri erisimi
// ===========================================================================
double MandrelGeometry::R_eq() const { return R_eq_; }
double MandrelGeometry::r0() const { return r0_; }
double MandrelGeometry::L_cyl() const { return L_cyl_; }
double MandrelGeometry::totalLength() const { return s_total_; }
DomeType MandrelGeometry::domeType() const { return dome_type_; }
const ProfileMetadata& MandrelGeometry::domeMetadata() const
{ return dome_table_.metadata(); }

} // namespace geometry
} // namespace filament
