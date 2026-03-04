// =============================================================================
// meridian_lookup_table.cpp — Kubik Spline Interpolasyonlu Meridyen Tablo
// =============================================================================
// Karar-9:  Onceden hesaplanmis lookup table + kubik spline interpolasyon
// Karar-11: Katman 2 toleranslari hedef (pozisyon < 1e-4 mm, turev < 1e-6)
// Karar-17: Alt katman — exception firlatilmaz, std::optional doner
//
// Algoritma:
//   build()  — veri kopyalama + dogal kubik spline katsayi hesabi (Thomas alg.)
//   query(s) — O(log n) binary search + Horner polinom degerlendirme
// =============================================================================

#include "geometry/meridian_lookup_table.h"
#include <algorithm>
#include <cmath>

namespace filament {
namespace geometry {

// ---------------------------------------------------------------------------
// build — Tabloyu doldur ve spline katsayilarini hesapla
// ---------------------------------------------------------------------------
void MeridianLookupTable::build(
    const std::vector<double>& s,
    const std::vector<double>& rho,
    const std::vector<double>& x_local,
    const std::vector<double>& drho_ds,
    const std::vector<double>& dx_ds,
    const std::vector<double>& kappa_m,
    const ProfileMetadata& metadata)
{
    valid_ = false;

    // --- Girdi dogrulama (alt katman: exception yok, sadece erken donus) ---

    const auto n = s.size();

    // Tum vektorler ayni boyutta olmali
    if (rho.size() != n || x_local.size() != n ||
        drho_ds.size() != n || dx_ds.size() != n || kappa_m.size() != n) {
        return;
    }

    // En az 2 veri noktasi gerekli
    if (n < 2) {
        return;
    }

    // s kesinlikle artan olmali
    for (std::size_t i = 1; i < n; ++i) {
        if (s[i] <= s[i - 1]) {
            return;
        }
    }

    // --- Veri kopyalama ---
    s_ = s;
    rho_ = rho;
    x_local_ = x_local;
    drho_ds_ = drho_ds;
    dx_ds_ = dx_ds;
    kappa_m_ = kappa_m;
    metadata_ = metadata;

    // --- Aralik genisliklerini hesapla ---
    const auto m = n - 1;  // aralik sayisi
    h_.resize(m);
    for (std::size_t i = 0; i < m; ++i) {
        h_[i] = s_[i + 1] - s_[i];
    }

    // --- 5 kanal icin kubik spline katsayilarini hesapla ---
    spline_rho_      = computeNaturalCubicSpline(rho_);
    spline_x_        = computeNaturalCubicSpline(x_local_);
    spline_drho_ds_  = computeNaturalCubicSpline(drho_ds_);
    spline_dx_ds_    = computeNaturalCubicSpline(dx_ds_);
    spline_kappa_m_  = computeNaturalCubicSpline(kappa_m_);

    valid_ = true;
}

// ---------------------------------------------------------------------------
// query — s parametresinden interpolasyon sorgusu
// ---------------------------------------------------------------------------
std::optional<MeridianPoint> MeridianLookupTable::query(double s) const
{
    if (!valid_ || s_.size() < 2) {
        return std::nullopt;
    }

    const double s_min = s_.front();
    const double s_max = s_.back();

    // Araligi kucuk toleransla genislet (sinirlardaki kayan nokta hatalari icin)
    constexpr double eps = 1e-10;
    if (s < s_min - eps || s > s_max + eps) {
        return std::nullopt;
    }

    // Gecerli araliga sabitle
    if (s < s_min) s = s_min;
    if (s > s_max) s = s_max;

    // O(log n) binary search ile aralik bul
    const auto seg = findSegment(s);

    // Lokal ofset
    const double t = s - s_[seg];

    // 5 kanalda kubik spline degerlendirme
    MeridianPoint pt;
    pt.s       = s;
    pt.rho     = evalPoly(spline_rho_,     seg, t);
    pt.x_local = evalPoly(spline_x_,       seg, t);
    pt.drho_ds = evalPoly(spline_drho_ds_, seg, t);
    pt.dx_ds   = evalPoly(spline_dx_ds_,   seg, t);
    pt.kappa_m = evalPoly(spline_kappa_m_, seg, t);

    return pt;
}

// ---------------------------------------------------------------------------
// computeNaturalCubicSpline — Thomas algoritmasi ile dogal kubik spline
// ---------------------------------------------------------------------------
// Dogal sinir kosullari: S''(s_0) = 0, S''(s_{n-1}) = 0
// Tridiagonal sistem n-2 bilinmeyenli (c_1, ..., c_{n-2})
// Referans: Numerical Recipes, 3rd Edition, Bölüm 3.3
// ---------------------------------------------------------------------------
MeridianLookupTable::SplineCoeffs
MeridianLookupTable::computeNaturalCubicSpline(const std::vector<double>& y) const
{
    const auto n = s_.size();
    const auto m = n - 1;  // aralik sayisi

    SplineCoeffs coeff;
    coeff.a.resize(m);
    coeff.b.resize(m);
    coeff.c.resize(m);
    coeff.d.resize(m);

    // --- n = 2: Lineer interpolasyon ---
    if (n == 2) {
        coeff.a[0] = y[0];
        coeff.b[0] = (y[1] - y[0]) / h_[0];
        coeff.c[0] = 0.0;
        coeff.d[0] = 0.0;
        return coeff;
    }

    // --- Genel durum: n >= 3 ---
    // c degerlerini hesapla (n adet, c_0 = c_{n-1} = 0 dogal BC)
    std::vector<double> c_all(n, 0.0);

    const auto interior = n - 2;  // ic nokta sayisi

    // Sag taraf vektoru (RHS)
    std::vector<double> rhs(interior);
    for (std::size_t i = 0; i < interior; ++i) {
        const auto j = i + 1;  // veri indeksi
        rhs[i] = 3.0 * ((y[j + 1] - y[j]) / h_[j]
                       - (y[j]     - y[j - 1]) / h_[j - 1]);
    }

    if (interior == 1) {
        // Tek denklem: 2*(h_0 + h_1) * c_1 = rhs[0]
        c_all[1] = rhs[0] / (2.0 * (h_[0] + h_[1]));
    } else {
        // Thomas algoritmasi (ileri tarama + geri yerine koyma)
        // Sistem satirlari (0-tabanli, i = 0..interior-1):
        //   alt kose: h_[i]         (i >= 1)
        //   ana kose: 2*(h_[i] + h_[i+1])
        //   ust kose: h_[i+1]       (i <= interior-2)
        std::vector<double> cp(interior, 0.0);  // degistirilmis ust kosegen
        std::vector<double> dp(interior);        // degistirilmis sag taraf

        // Ilk satir (i=0): alt kosegen yok (c_0 = 0 zaten elimine)
        double diag = 2.0 * (h_[0] + h_[1]);
        cp[0] = h_[1] / diag;
        dp[0] = rhs[0] / diag;

        // Ileri tarama
        for (std::size_t i = 1; i < interior; ++i) {
            const double sub = h_[i];  // alt kosegen elemani
            diag = 2.0 * (h_[i] + h_[i + 1]);
            const double sup = (i + 1 < interior) ? h_[i + 1] : 0.0;

            const double denom = diag - sub * cp[i - 1];
            cp[i] = sup / denom;
            dp[i] = (rhs[i] - sub * dp[i - 1]) / denom;
        }

        // Geri yerine koyma
        c_all[interior] = dp[interior - 1];  // c_{n-2}
        for (int i = static_cast<int>(interior) - 2; i >= 0; --i) {
            c_all[static_cast<std::size_t>(i) + 1] =
                dp[static_cast<std::size_t>(i)]
                - cp[static_cast<std::size_t>(i)]
                  * c_all[static_cast<std::size_t>(i) + 2];
        }
    }

    // --- a, b, c, d katsayilarini hesapla ---
    for (std::size_t i = 0; i < m; ++i) {
        coeff.a[i] = y[i];
        coeff.d[i] = (c_all[i + 1] - c_all[i]) / (3.0 * h_[i]);
        coeff.b[i] = (y[i + 1] - y[i]) / h_[i]
                   - h_[i] * (2.0 * c_all[i] + c_all[i + 1]) / 3.0;
        coeff.c[i] = c_all[i];
    }

    return coeff;
}

// ---------------------------------------------------------------------------
// findSegment — Binary search ile aralik indeksi bul
// ---------------------------------------------------------------------------
// s_[i] <= s <= s_[i+1] olan i'yi doner.
// s == s_[n-1] (son nokta) icin son araligi (n-2) doner.
// ---------------------------------------------------------------------------
std::size_t MeridianLookupTable::findSegment(double s) const
{
    // std::upper_bound: s'den buyuk ilk elemani bul
    auto it = std::upper_bound(s_.begin(), s_.end(), s);

    if (it == s_.begin()) {
        return 0;
    }
    if (it == s_.end()) {
        return s_.size() - 2;
    }

    return static_cast<std::size_t>(std::distance(s_.begin(), it)) - 1;
}

// ---------------------------------------------------------------------------
// evalPoly — Horner yontemiyle kubik polinom degerlendirme
// ---------------------------------------------------------------------------
// S_i(t) = a + t*(b + t*(c + t*d))
// ---------------------------------------------------------------------------
double MeridianLookupTable::evalPoly(
    const SplineCoeffs& coeff, std::size_t i, double t)
{
    return coeff.a[i] + t * (coeff.b[i] + t * (coeff.c[i] + t * coeff.d[i]));
}

// ---------------------------------------------------------------------------
// Accessor'lar
// ---------------------------------------------------------------------------
const ProfileMetadata& MeridianLookupTable::metadata() const { return metadata_; }
std::size_t MeridianLookupTable::size() const { return s_.size(); }
bool MeridianLookupTable::isValid() const { return valid_; }
const std::vector<double>& MeridianLookupTable::rawS() const { return s_; }
const std::vector<double>& MeridianLookupTable::rawRho() const { return rho_; }
const std::vector<double>& MeridianLookupTable::rawX() const { return x_local_; }

} // namespace geometry
} // namespace filament
