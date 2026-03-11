// =============================================================================
// meridian_lookup_table.h — Meridyen Profil Veri Tablosu ve Interpolasyon
// =============================================================================
// Karar-9:  Onceden hesaplanmis lookup table + interpolasyon
// Karar-11: Tolerans kriterleri (Katman 2)
//
// Tablo yapisi: {s_i, rho_i, x_i, drho/ds_i, dx/ds_i, kappa_m_i}
// Interpolasyon: Kubik spline (GATE-1a secimi) — dogal sinir kosullari
// Sorgulama: O(log n) binary search + interpolasyon
// Hata yonetimi: Alt katman — std::optional (Karar-17), exception yok
// =============================================================================

#ifndef FILAMENT_GEOMETRY_MERIDIAN_LOOKUP_TABLE_H
#define FILAMENT_GEOMETRY_MERIDIAN_LOOKUP_TABLE_H

#include "geometry/filament_types.h"
#include <vector>
#include <optional>
#include <cstddef>

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

// Ters sorgu sonucu: rho -> (s, x)
struct InversePoint {
    double s;        // Yay uzunlugu [mm]
    double x_local;  // Aksiyel konum [mm]
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
    // Onsart: tum vektorler ayni boyutta, en az 2 eleman, s kesinlikle artan
    // Gecersiz girdi durumunda tablo invalid kalir (isValid() == false)
    void build(const std::vector<double>& s,
               const std::vector<double>& rho,
               const std::vector<double>& x_local,
               const std::vector<double>& drho_ds,
               const std::vector<double>& dx_ds,
               const std::vector<double>& kappa_m,
               const ProfileMetadata& metadata);

    // s parametresinden interpolasyon sorgusu — O(log n) binary search + kubik spline
    // Alt katman (Karar-17): exception firlatmaz, std::optional doner
    // s araligi disinda veya gecersiz tablo icin std::nullopt doner
    std::optional<MeridianPoint> query(double s) const;

    // Ters sorgu: rho -> (s, x_local) — PCHIP interpolasyon
    // Monoton azalan rho(s) profili icin; rho araligi disinda std::nullopt doner
    // Phase-2a S3: geodesic yol hesabinda rho'dan s ve x geri cozumu
    std::optional<InversePoint> inverseLookup(double rho) const;

    // Meta-veri erisimi
    const ProfileMetadata& metadata() const;

    // Tablo boyutu (veri noktasi sayisi)
    std::size_t size() const;

    // Tablo gecerli mi? (build basarili mi tamamlandi?)
    bool isValid() const;

    // Ham veri erisimi (test amacli)
    const std::vector<double>& rawS() const;
    const std::vector<double>& rawRho() const;
    const std::vector<double>& rawX() const;

private:
    // Ham veri noktalari
    std::vector<double> s_;
    std::vector<double> rho_;
    std::vector<double> x_local_;
    std::vector<double> drho_ds_;
    std::vector<double> dx_ds_;
    std::vector<double> kappa_m_;

    ProfileMetadata metadata_{};
    bool valid_ = false;

    // -----------------------------------------------------------------------
    // Kubik spline altyapisi
    // -----------------------------------------------------------------------
    // Dogal kubik spline: S_i(t) = a + b*t + c*t^2 + d*t^3
    // t = s - s_[i], her aralik icin ayri katsayilar
    struct SplineCoeffs {
        std::vector<double> a;  // sabit terim    (aralik sayisi = n-1)
        std::vector<double> b;  // lineer terim
        std::vector<double> c;  // kuadratik terim
        std::vector<double> d;  // kubik terim
    };

    // 5 kanal icin spline katsayilari
    SplineCoeffs spline_rho_;
    SplineCoeffs spline_x_;
    SplineCoeffs spline_drho_ds_;
    SplineCoeffs spline_dx_ds_;
    SplineCoeffs spline_kappa_m_;

    // Aralik genislikleri: h_[i] = s_[i+1] - s_[i]
    std::vector<double> h_;

    // Dogal kubik spline katsayilarini hesapla (Thomas algoritmasi)
    // Sinir kosullari: S''(s_0) = 0, S''(s_{n-1}) = 0
    SplineCoeffs computeNaturalCubicSpline(const std::vector<double>& y) const;

    // Binary search: s_[i] <= s <= s_[i+1] olan aralik indeksini bul
    std::size_t findSegment(double s) const;

    // Kubik polinom degerlendirme (Horner yontemi)
    static double evalPoly(const SplineCoeffs& coeff, std::size_t i, double t);

    // -----------------------------------------------------------------------
    // Ters sorgu altyapisi (forward spline + bisection)
    // -----------------------------------------------------------------------
    // rho monoton azalir: rho_[0] = R_eq > rho_[n-1] = r0
    // inverseLookup: rho -> binary search + bisection on rho(s) spline -> s, x
    // Forward spline hassasiyetini (Karar-11 Katman 2) korur.
};

} // namespace geometry
} // namespace filament

#endif // FILAMENT_GEOMETRY_MERIDIAN_LOOKUP_TABLE_H
