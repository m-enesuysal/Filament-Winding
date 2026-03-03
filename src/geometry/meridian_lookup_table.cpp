// =============================================================================
// meridian_lookup_table.cpp — Meridyen Profil Veri Tablosu Implementasyonu
// =============================================================================
// Karar-9:  Onceden hesaplanmis lookup table + interpolasyon
// Karar-11: Tolerans kriterleri (Katman 2)
// Karar-17: Alt katman hata yonetimi — std::optional, exception yok
//
// Interpolasyon yontemi: Clamped kubik spline (GATE-1a secimi)
//   - rho, x_local: endpoint turevleri drho_ds ve dx_ds'den alinir
//   - drho_ds, dx_ds: diferansiyel geometri formulleri ile endpoint turev
//       d(drho_ds)/ds = -dx_ds * kappa_m  (Frenet-Serret, meridyen)
//       d(dx_ds)/ds   = drho_ds * kappa_m
//   - kappa_m: dogal spline (endpoint turev bulunmuyor)
//   - Sinir kosullari: S'(s_0) = f'_0, S'(s_n) = f'_n (clamped)
//   - Thomas algoritmasi ile O(n) katsayi hesabi
//   - O(log n) binary search ile segment tespiti
//
// Karar-19: -ffast-math YASAK (IEEE 754 zorunlu)
// =============================================================================

#include "geometry/meridian_lookup_table.h"

#include <algorithm>   // std::upper_bound
#include <cassert>     // assert
#include <cmath>       // std::isfinite

namespace filament {
namespace geometry {

// =============================================================================
// buildNaturalSpline — Dogal Kubik Spline Katsayi Hesabi
// =============================================================================
// Algoritma: Thomas (tridiagonal matrix) algoritmasi
//
// n+1 dugum (knots[0..n], values[0..n]) → n interval katsayisi (b,c,d)
//
// Her interval [knots[i], knots[i+1]] icin kubik polinom:
//   S_i(s) = values[i] + b[i]*dt + c[i]*dt^2 + d[i]*dt^3
//   dt = s - knots[i]
//
// Moment M icin tridiagonal sistem:
//   h[i-1]*M[i-1] + 2*(h[i-1]+h[i])*M[i] + h[i]*M[i+1] = rhs[i]
//   Sinir: M[0] = M[n] = 0 (dogal spline)
// =============================================================================
MeridianLookupTable::CubicSpline MeridianLookupTable::buildNaturalSpline(
    const std::vector<double>& knots,
    const std::vector<double>& values)
{
    assert(knots.size() == values.size());

    const std::size_t n = knots.size() - 1;  // interval sayisi

    CubicSpline sp;
    sp.b.resize(n, 0.0);
    sp.c.resize(n, 0.0);
    sp.d.resize(n, 0.0);

    if (n == 0) return sp;

    std::vector<double> h(n);
    for (std::size_t i = 0; i < n; ++i) {
        h[i] = knots[i + 1] - knots[i];
        assert(h[i] > 0.0);
    }

    if (n == 1) {
        sp.b[0] = (values[1] - values[0]) / h[0];
        return sp;
    }

    // Birinci dereceden fark bolumler
    std::vector<double> delta(n);
    for (std::size_t i = 0; i < n; ++i) {
        delta[i] = (values[i + 1] - values[i]) / h[i];
    }

    // Thomas forward sweep (ic dugumler: M[1..n-1])
    const std::size_t m = n - 1;
    std::vector<double> cp(m, 0.0), rp(m, 0.0);

    {
        double w = 2.0 * (h[0] + h[1]);
        cp[0] = (m > 1) ? (h[1] / w) : 0.0;
        rp[0] = 6.0 * (delta[1] - delta[0]) / w;
    }
    for (std::size_t k = 1; k < m; ++k) {
        double w = 2.0 * (h[k] + h[k + 1]) - h[k] * cp[k - 1];
        cp[k] = (k < m - 1) ? (h[k + 1] / w) : 0.0;
        rp[k] = (6.0 * (delta[k + 1] - delta[k]) - h[k] * rp[k - 1]) / w;
    }

    // Thomas back substitution
    std::vector<double> M(n + 1, 0.0);
    M[m] = rp[m - 1];
    for (std::size_t k = m - 1; k > 0; --k) {
        M[k] = rp[k - 1] - cp[k - 1] * M[k + 1];
    }

    for (std::size_t i = 0; i < n; ++i) {
        sp.b[i] = delta[i] - h[i] * (2.0 * M[i] + M[i + 1]) / 6.0;
        sp.c[i] = M[i] / 2.0;
        sp.d[i] = (M[i + 1] - M[i]) / (6.0 * h[i]);
    }

    return sp;
}

// =============================================================================
// buildClampedSpline — Clamped Kubik Spline Katsayi Hesabi
// =============================================================================
// Sinir kosullari: S'(x_0) = slope_left, S'(x_n) = slope_right
//
// Dogal spline'dan farki: M[0] ve M[n] bilinmeyen olarak dahil edilir.
// n+1 dugum → n+1 bilinmeyen (M[0]..M[n]) → (n+1)x(n+1) tridiagonal sistem.
//
// Ek denklemler:
//   Sol BC:   2h[0]*M[0] + h[0]*M[1] = 6*(delta[0] - slope_left)
//   Sag BC:   h[n-1]*M[n-1] + 2h[n-1]*M[n] = 6*(slope_right - delta[n-1])
//
// Karar-11: Clamped BC endpoint hatalarini giderir (O(h^4) yerine O(h^3)
//   endpoint davranisi dogal spline'da; clamped O(h^4) her yerde).
// =============================================================================
MeridianLookupTable::CubicSpline MeridianLookupTable::buildClampedSpline(
    const std::vector<double>& knots,
    const std::vector<double>& values,
    double slope_left,
    double slope_right)
{
    assert(knots.size() == values.size());

    const std::size_t n = knots.size() - 1;

    CubicSpline sp;
    sp.b.resize(n, 0.0);
    sp.c.resize(n, 0.0);
    sp.d.resize(n, 0.0);

    if (n == 0) return sp;

    std::vector<double> h(n);
    for (std::size_t i = 0; i < n; ++i) {
        h[i] = knots[i + 1] - knots[i];
        assert(h[i] > 0.0);
    }

    if (n == 1) {
        // Tek interval: clamped Hermite cubici
        // S(t) = values[0] + slope_left*t + c*t^2 + d*t^3
        // S'(0) = slope_left (tatmin)
        // S'(h) = slope_right: slope_left + 2c*h + 3d*h^2
        // S(h) = values[1]: values[0] + slope_left*h + c*h^2 + d*h^3
        double dy = values[1] - values[0];
        double c = (3.0*dy/h[0] - 2.0*slope_left - slope_right) / h[0];
        double d = (slope_left + slope_right - 2.0*dy/h[0]) / (h[0]*h[0]);
        sp.b[0] = slope_left;
        sp.c[0] = c;
        sp.d[0] = d;
        return sp;
    }

    // Birinci dereceden fark bolumler
    std::vector<double> delta(n);
    for (std::size_t i = 0; i < n; ++i) {
        delta[i] = (values[i + 1] - values[i]) / h[i];
    }

    // -----------------------------------------------------------------------
    // (n+1)x(n+1) tridiagonal sistem: M[0..n]
    //
    // Alt diagonal: l[1]=h[0], l[2]=h[1], ..., l[n]=h[n-1]
    // Ana diagonal: d[0]=2h[0], d[i]=2(h[i-1]+h[i]) (i=1..n-1), d[n]=2h[n-1]
    // Ust diagonal: u[0]=h[0], u[1]=h[1], ..., u[n-1]=h[n-1]
    // Sag taraf:
    //   rhs[0]   = 6*(delta[0] - slope_left)
    //   rhs[i]   = 6*(delta[i] - delta[i-1])   (i=1..n-1)
    //   rhs[n]   = 6*(slope_right - delta[n-1])
    // -----------------------------------------------------------------------
    const std::size_t sz = n + 1;  // system size
    std::vector<double> cp(sz, 0.0), rp(sz, 0.0);

    // Forward sweep, satir 0 (sol clamped BC)
    {
        double diag0 = 2.0 * h[0];
        cp[0] = h[0] / diag0;  // u[0] / d[0]
        rp[0] = 6.0 * (delta[0] - slope_left) / diag0;
    }

    // Forward sweep, ic satirlar i=1..n-1
    for (std::size_t i = 1; i < n; ++i) {
        // Ana diagonal: 2*(h[i-1]+h[i])
        // Alt diagonal: h[i-1]
        // Ust diagonal: h[i]
        double w = 2.0 * (h[i - 1] + h[i]) - h[i - 1] * cp[i - 1];
        cp[i] = h[i] / w;
        rp[i] = (6.0 * (delta[i] - delta[i - 1]) - h[i - 1] * rp[i - 1]) / w;
    }

    // Forward sweep, satir n (sag clamped BC)
    {
        // Ana diagonal: 2*h[n-1]
        // Alt diagonal: h[n-1]
        double w = 2.0 * h[n - 1] - h[n - 1] * cp[n - 1];
        cp[n] = 0.0;  // son satirin ustdiyagonal yok
        rp[n] = (6.0 * (slope_right - delta[n - 1]) - h[n - 1] * rp[n - 1]) / w;
    }

    // Back substitution
    std::vector<double> M(sz);
    M[n] = rp[n];
    for (std::size_t i = n; i > 0; --i) {
        M[i - 1] = rp[i - 1] - cp[i - 1] * M[i];
    }

    // Spline katsayilari
    for (std::size_t i = 0; i < n; ++i) {
        sp.b[i] = delta[i] - h[i] * (2.0 * M[i] + M[i + 1]) / 6.0;
        sp.c[i] = M[i] / 2.0;
        sp.d[i] = (M[i + 1] - M[i]) / (6.0 * h[i]);
    }

    return sp;
}

// =============================================================================
// findSegment — Binary Search ile Segment Indeksi
// =============================================================================
std::size_t MeridianLookupTable::findSegment(
    const std::vector<double>& knots,
    double s_query)
{
    auto it = std::upper_bound(knots.begin(), knots.end(), s_query);
    std::size_t idx = static_cast<std::size_t>(it - knots.begin());

    if (idx == 0) idx = 1;
    if (idx >= knots.size()) idx = knots.size() - 1;

    return idx - 1;
}

// =============================================================================
// evalSpline — Spline Degerlendirmesi (Horner yontemi)
// =============================================================================
double MeridianLookupTable::evalSpline(
    const std::vector<double>& values,
    const CubicSpline& sp,
    std::size_t seg,
    double ds)
{
    return values[seg] + ds * (sp.b[seg] + ds * (sp.c[seg] + ds * sp.d[seg]));
}

// =============================================================================
// build — Tabloyu Doldur ve Spline Katsayilarini Hesapla
// =============================================================================
// Clamped spline stratejisi:
//   rho     : sol/sag endpoint turev = drho_ds[0] / drho_ds.back()
//   x_local : sol/sag endpoint turev = dx_ds[0]   / dx_ds.back()
//   drho_ds : sol = -dx_ds[0]*kappa_m[0], sag = -dx_ds.back()*kappa_m.back()
//             (diferansiyel geometri: d(drho_ds)/ds = -dx_ds * kappa_m)
//   dx_ds   : sol = drho_ds[0]*kappa_m[0], sag = drho_ds.back()*kappa_m.back()
//             (diferansiyel geometri: d(dx_ds)/ds = drho_ds * kappa_m)
//   kappa_m : dogal spline (endpoint turev bulunmuyor)
// =============================================================================
void MeridianLookupTable::build(
    const std::vector<double>& s,
    const std::vector<double>& rho,
    const std::vector<double>& x_local,
    const std::vector<double>& drho_ds,
    const std::vector<double>& dx_ds,
    const std::vector<double>& kappa_m,
    const ProfileMetadata& metadata)
{
    if (s.size() < 2) {
        valid_ = false;
        return;
    }

    const std::size_t n = s.size();
    if (rho.size() != n || x_local.size() != n ||
        drho_ds.size() != n || dx_ds.size() != n || kappa_m.size() != n) {
        valid_ = false;
        return;
    }

    s_        = s;
    rho_      = rho;
    x_local_  = x_local;
    drho_ds_  = drho_ds;
    dx_ds_    = dx_ds;
    kappa_m_  = kappa_m;
    metadata_ = metadata;

    // Clamped endpoint turevleri
    const double drho_left  = drho_ds_.front();
    const double drho_right = drho_ds_.back();
    const double dx_left    = dx_ds_.front();
    const double dx_right   = dx_ds_.back();
    const double km_left    = kappa_m_.front();
    const double km_right   = kappa_m_.back();

    // d(drho_ds)/ds = -dx_ds * kappa_m  (Frenet-Serret, genel meridyen)
    const double ddrho_left  = -dx_left  * km_left;
    const double ddrho_right = -dx_right * km_right;

    // d(dx_ds)/ds = drho_ds * kappa_m  (Frenet-Serret, genel meridyen)
    const double ddx_left    =  drho_left  * km_left;
    const double ddx_right   =  drho_right * km_right;

    // Clamped spline: endpoint hatalarini giderir
    sp_rho_   = buildClampedSpline(s_, rho_,     drho_left,  drho_right);
    sp_x_     = buildClampedSpline(s_, x_local_,  dx_left,   dx_right);
    sp_drho_  = buildClampedSpline(s_, drho_ds_, ddrho_left, ddrho_right);
    sp_dx_    = buildClampedSpline(s_, dx_ds_,   ddx_left,   ddx_right);

    // kappa_m: dogal spline (hemisphere icin sabit, genel profil icin yeterli)
    sp_kappa_ = buildNaturalSpline(s_, kappa_m_);

    valid_ = true;
}

// =============================================================================
// query — Interpolasyon Sorgusu (Alt katman: std::optional, exception yok)
// =============================================================================
std::optional<MeridianPoint> MeridianLookupTable::query(double s_val) const
{
    if (!valid_ || s_.empty()) return std::nullopt;

    constexpr double eps = 1e-12;
    if (s_val < s_.front() - eps || s_val > s_.back() + eps) {
        return std::nullopt;
    }

    const double s_clamped = std::min(std::max(s_val, s_.front()), s_.back());
    const std::size_t seg  = findSegment(s_, s_clamped);
    const double ds        = s_clamped - s_[seg];

    MeridianPoint pt;
    pt.s       = s_clamped;
    pt.rho     = evalSpline(rho_,     sp_rho_,   seg, ds);
    pt.x_local = evalSpline(x_local_, sp_x_,     seg, ds);
    pt.drho_ds = evalSpline(drho_ds_, sp_drho_,  seg, ds);
    pt.dx_ds   = evalSpline(dx_ds_,   sp_dx_,    seg, ds);
    pt.kappa_m = evalSpline(kappa_m_, sp_kappa_, seg, ds);

    return pt;
}

// =============================================================================
// Accessor metodlar
// =============================================================================
const ProfileMetadata& MeridianLookupTable::metadata() const { return metadata_; }
std::size_t            MeridianLookupTable::size()     const { return s_.size(); }
bool                   MeridianLookupTable::isValid()  const { return valid_;    }

const std::vector<double>& MeridianLookupTable::rawS()   const { return s_;       }
const std::vector<double>& MeridianLookupTable::rawRho() const { return rho_;     }
const std::vector<double>& MeridianLookupTable::rawX()   const { return x_local_; }

} // namespace geometry
} // namespace filament
