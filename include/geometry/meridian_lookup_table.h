// =============================================================================
// meridian_lookup_table.h — Meridyen Profil Veri Tablosu ve Interpolasyon
// =============================================================================
// Karar-9:  Onceden hesaplanmis lookup table + interpolasyon
// Karar-11: Tolerans kriterleri (Katman 2)
//
// Tablo yapisi: {s_i, rho_i, x_i, drho/ds_i, dx/ds_i, kappa_m_i}
// Interpolasyon: Kubik spline (GATE-1a secimi)
// Sorgulama: O(log n) binary search + interpolasyon
// Hata yonetimi: Alt katman — std::optional (Karar-17)
//
// Detayli implementasyon: Phase-1b S2 session'inda yapilacak.
// =============================================================================

#ifndef FILAMENT_GEOMETRY_MERIDIAN_LOOKUP_TABLE_H
#define FILAMENT_GEOMETRY_MERIDIAN_LOOKUP_TABLE_H

#include "geometry/filament_types.h"
#include <vector>
#include <optional>

namespace filament {
namespace geometry {

// Interpolasyon sorgu sonucu
struct MeridianPoint {
    double s;        // Yay uzunlugu [mm]
    double rho;      // Yaricap [mm]
    double x_local;  // Aksiyel konum [mm]
    double drho_ds;  // drho/ds [-]
    double dx_ds;    // dx/ds [-]
    double kappa_m;  // Meridyen egriligi [1/mm]
};

// Skaler profil meta-verisi (struct arayuz yamasi uyumlu)
struct ProfileMetadata {
    double R_eq      = 0.0;  // Ekvator yaricapi [mm]
    double r0        = 0.0;  // Polar aciklik yaricapi [mm]
    double s_total   = 0.0;  // Toplam yay uzunlugu [mm]
    double h_dome    = 0.0;  // Dome yuksekligi [mm]
    double A_dome    = 0.0;  // Dome yuzey alani [mm^2]
    double kappa_eq  = 0.0;  // Ekvator meridyen egriligi [1/mm]
    double kappa_pol = 0.0;  // Polar aciklik meridyen egriligi [1/mm]
    double alpha_w   = 0.0;  // Ekvator winding acisi [rad]
    double aspect_r  = 0.0;  // Dome aspect ratio [-]
};

class MeridianLookupTable {
public:
    MeridianLookupTable() = default;

    // Tabloyu doldur — profil ureticisi tarafindan cagrilir
    void build(const std::vector<double>& s,
               const std::vector<double>& rho,
               const std::vector<double>& x_local,
               const std::vector<double>& drho_ds,
               const std::vector<double>& dx_ds,
               const std::vector<double>& kappa_m,
               const ProfileMetadata& metadata);

    // s parametresinden interpolasyon sorgusu
    // Alt katman (Karar-17): exception firlatmaz, std::optional doner
    std::optional<MeridianPoint> query(double s) const;

    // Meta-veri erisimi
    const ProfileMetadata& metadata() const;

    // Tablo boyutu
    std::size_t size() const;

    // Tablo dolu mu?
    bool isValid() const;

    // Ham veri erisimi (test amacli)
    const std::vector<double>& rawS() const;
    const std::vector<double>& rawRho() const;
    const std::vector<double>& rawX() const;

private:
    std::vector<double> s_;
    std::vector<double> rho_;
    std::vector<double> x_local_;
    std::vector<double> drho_ds_;
    std::vector<double> dx_ds_;
    std::vector<double> kappa_m_;

    ProfileMetadata metadata_{};
    bool valid_ = false;

    // -----------------------------------------------------------------------
    // Kubik spline katsayi yapisi
    // S(s) = values[i] + b[i]*dt + c[i]*dt^2 + d[i]*dt^3
    // dt = s - s_[i], i: segment indeksi, a[i] = values[i] (harici)
    // n+1 dugum icin n adet interval katsayisi
    // -----------------------------------------------------------------------
    struct CubicSpline {
        std::vector<double> b;  // Lineer katsayi
        std::vector<double> c;  // Karesel katsayi
        std::vector<double> d;  // Kubik katsayi
    };

    CubicSpline sp_rho_;     // rho(s) spline
    CubicSpline sp_x_;       // x_local(s) spline
    CubicSpline sp_drho_;    // drho_ds(s) spline
    CubicSpline sp_dx_;      // dx_ds(s) spline
    CubicSpline sp_kappa_;   // kappa_m(s) spline

    // -----------------------------------------------------------------------
    // Dogal kubik spline katsayi hesabi (Thomas algoritmasi)
    // Sinir kosullari: M_0 = M_n = 0 (ikinci turev sinirda sifir)
    // knots ve values n+1 elemanli olmali
    // -----------------------------------------------------------------------
    static CubicSpline buildNaturalSpline(
        const std::vector<double>& knots,
        const std::vector<double>& values);

    // -----------------------------------------------------------------------
    // Clamped kubik spline katsayi hesabi (Thomas algoritmasi)
    // Sinir kosullari: S'(x_0) = slope_left, S'(x_n) = slope_right
    // knots ve values n+1 elemanli olmali
    // Karar-9: Dogrulugu bilinen endpoint turevleriyle hata azaltimi
    // -----------------------------------------------------------------------
    static CubicSpline buildClampedSpline(
        const std::vector<double>& knots,
        const std::vector<double>& values,
        double slope_left,
        double slope_right);

    // -----------------------------------------------------------------------
    // Binary search ile segment indeksi — O(log n)
    // Dondurur: i oyle ki knots[i] <= s_query <= knots[i+1]
    // Sinir durumu: s_query == knots.back() → son segment indeksi
    // -----------------------------------------------------------------------
    static std::size_t findSegment(
        const std::vector<double>& knots,
        double s_query);

    // -----------------------------------------------------------------------
    // Tek buyukluk icin spline degerlendirmesi
    // seg: findSegment() ciktisi, ds: s - s_[seg]
    // -----------------------------------------------------------------------
    static double evalSpline(
        const std::vector<double>& values,
        const CubicSpline& sp,
        std::size_t seg,
        double ds);
};

} // namespace geometry
} // namespace filament

#endif // FILAMENT_GEOMETRY_MERIDIAN_LOOKUP_TABLE_H
