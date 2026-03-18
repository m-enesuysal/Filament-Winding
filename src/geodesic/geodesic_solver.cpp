// =============================================================================
// geodesic_solver.cpp — Geodesic Sarim Yolu ODE Cozucu
// =============================================================================
// Phase-2a S4: Tek devre geodesic yol hesabi
//
// ODE: dphi/ds = R_E / (rho(s) * sqrt(rho(s)^2 - R_E^2))
// Singularite: birlesik 1. + 2. derece analitik kuyruk (dome_phi_integration.m)
// Numerik: Dormand-Prince RK45 (isotensoid_profile.cpp ile ayni Butcher tablosu)
//
// Phase-1b API kullanimi:
//   MandrelGeometry -> dome lookup table -> rho(s), drho/ds(s) sorgusu
//   MeridianLookupTable::inverseLookup(rho) -> s_turn, s_epsilon tespiti
// =============================================================================

#include "geodesic/geodesic_solver.h"
#include "geometry/mandrel_geometry.h"
#include "geometry/meridian_lookup_table.h"
#include <cmath>
#include <algorithm>
#include <array>

namespace filament {
namespace geodesic {

namespace {

// ============================================================================
// Dormand-Prince RK45 Butcher Tablosu (isotensoid_profile.cpp ile ayni)
// ============================================================================
namespace dp {
    constexpr double c2 = 1.0 / 5.0;
    constexpr double c3 = 3.0 / 10.0;
    constexpr double c4 = 4.0 / 5.0;
    constexpr double c5 = 8.0 / 9.0;

    constexpr double a21 = 1.0 / 5.0;
    constexpr double a31 = 3.0 / 40.0;
    constexpr double a32 = 9.0 / 40.0;
    constexpr double a41 = 44.0 / 45.0;
    constexpr double a42 = -56.0 / 15.0;
    constexpr double a43 = 32.0 / 9.0;
    constexpr double a51 = 19372.0 / 6561.0;
    constexpr double a52 = -25360.0 / 2187.0;
    constexpr double a53 = 64448.0 / 6561.0;
    constexpr double a54 = -212.0 / 729.0;
    constexpr double a61 = 9017.0 / 3168.0;
    constexpr double a62 = -355.0 / 33.0;
    constexpr double a63 = 46732.0 / 5247.0;
    constexpr double a64 = 49.0 / 176.0;
    constexpr double a65 = -5103.0 / 18656.0;

    constexpr double b1 = 35.0 / 384.0;
    constexpr double b3 = 500.0 / 1113.0;
    constexpr double b4 = 125.0 / 192.0;
    constexpr double b5 = -2187.0 / 6784.0;
    constexpr double b6 = 11.0 / 84.0;

    constexpr double e1 = 71.0 / 57600.0;
    constexpr double e3 = -71.0 / 16695.0;
    constexpr double e4 = 71.0 / 1920.0;
    constexpr double e5 = -17253.0 / 339200.0;
    constexpr double e6 = 22.0 / 525.0;
    constexpr double e7 = -1.0 / 40.0;
} // namespace dp

// ============================================================================
// Skaler RK45 ODE cozucu: dphi/ds = f(s)
// s: dome-lokal yay uzunlugu, phi: azimuthal aci
// ============================================================================
struct PhiOdeSolution {
    std::vector<double> s;
    std::vector<double> phi;
};

PhiOdeSolution solvePhiOde(
    const geometry::MeridianLookupTable& dome_table,
    double R_E,
    double s_end,       // s_epsilon (numerik durdurma noktasi)
    double rel_tol,
    double abs_tol)
{
    constexpr double SAFETY    = 0.9;
    constexpr double H_MIN     = 1e-15;
    constexpr double GROW_MAX  = 5.0;
    constexpr double SHRINK_MAX = 0.2;
    constexpr int    MAX_STEPS = 200000;

    const double h_max = s_end / 100.0;

    // ODE sag taraf: dphi/ds = R_E / (rho * sqrt(rho^2 - R_E^2))
    auto rhs = [&](double s_dome) -> double {
        auto pt = dome_table.query(s_dome);
        if (!pt) return 0.0;
        const double rho = pt->rho;
        const double diff = rho * rho - R_E * R_E;
        if (diff <= 0.0) return 0.0;
        return R_E / (rho * std::sqrt(diff));
    };

    double s = 0.0;
    double phi = 0.0;
    double h = s_end / 200.0;

    PhiOdeSolution sol;
    sol.s.reserve(5000);
    sol.phi.reserve(5000);
    sol.s.push_back(s);
    sol.phi.push_back(phi);

    for (int step = 0; step < MAX_STEPS && s < s_end - H_MIN; ++step) {
        if (s + h > s_end) h = s_end - s;
        if (h < H_MIN) h = H_MIN;

        // 7-stage Dormand-Prince (skaler versiyon)
        // NOT: dphi/ds = f(s) sadece s'ye bagimli, phi'ya degil.
        //      Butcher tablosunda b2=0, dolayisiyla k2 hesaplanmaz.
        const double k1 = rhs(s);

        const double k3 = rhs(s + dp::c3 * h);
        const double k4 = rhs(s + dp::c4 * h);
        const double k5 = rhs(s + dp::c5 * h);
        const double k6 = rhs(s + h);

        // 5. derece cozum
        const double phi_new = phi + h * (dp::b1 * k1 + dp::b3 * k3
                                         + dp::b4 * k4 + dp::b5 * k5
                                         + dp::b6 * k6);

        // k7 (hata tahmini)
        const double k7 = rhs(s + h);

        // Hata tahmini
        const double err_raw = h * (dp::e1 * k1 + dp::e3 * k3
                                   + dp::e4 * k4 + dp::e5 * k5
                                   + dp::e6 * k6 + dp::e7 * k7);
        const double sc = abs_tol + rel_tol * std::max(std::abs(phi), std::abs(phi_new));
        double err = std::abs(err_raw) / sc;

        if (err <= 1.0) {
            s += h;
            phi = phi_new;
            sol.s.push_back(s);
            sol.phi.push_back(phi);
        }

        // Adim boyutu ayarlama
        double factor;
        if (err < 1e-30) {
            factor = GROW_MAX;
        } else {
            factor = SAFETY * std::pow(err, -0.2);
            factor = std::min(factor, GROW_MAX);
            factor = std::max(factor, SHRINK_MAX);
        }
        h *= factor;
        h = std::min(h, h_max);
    }

    return sol;
}

// ============================================================================
// Analitik kuyruk tamamlama — 3 rejim (dome_phi_integration.m birebir)
// ============================================================================
struct TailResult {
    double phi_tail;
    double delta_eps;
    double a;   // |drho/ds| at s_epsilon
    double b;   // |d2rho/ds2| at s_epsilon
};

TailResult computeAnalyticTail(
    const geometry::MeridianLookupTable& dome_table,
    double R_E,
    double epsilon,
    double s_epsilon,
    double s_total)
{
    // 1. turev: a = |drho/ds| at s_epsilon
    auto pt_eps = dome_table.query(s_epsilon);
    const double a = std::abs(pt_eps->drho_ds);

    // 2. turev: b = |d2rho/ds2| at s_epsilon (merkezi sonlu fark)
    double ds_fd = std::min(epsilon / 2.0, (s_total - s_epsilon) / 4.0);
    ds_fd = std::max(ds_fd, 1e-8);

    auto pt_plus  = dome_table.query(s_epsilon + ds_fd);
    auto pt_minus = dome_table.query(s_epsilon - ds_fd);
    const double d2rho = (pt_plus->drho_ds - pt_minus->drho_ds) / (2.0 * ds_fd);
    const double b = std::abs(d2rho);

    TailResult result{};
    result.a = a;
    result.b = b;

    // --- Birlesik kuyruk formulu (3 rejim) ---
    if (b < 1e-12 && a > 1e-12) {
        // Rejim 1: Saf birinci derece (standart dome'lar)
        result.delta_eps = epsilon / a;
        result.phi_tail  = (1.0 / a) * std::sqrt(2.0 * epsilon / R_E);

    } else if (a < 1e-12 && b > 1e-12) {
        // Rejim 2: Saf ikinci derece (isotensoid donus noktasi, a ~ 0)
        result.delta_eps = std::sqrt(2.0 * epsilon / b);
        result.phi_tail  = (1.0 / std::sqrt(R_E * b)) * std::acosh(1.0 + epsilon / R_E);

    } else {
        // Rejim 3: Genel birlesik formul
        const double discriminant = a * a + 2.0 * b * epsilon;
        result.delta_eps = 2.0 * epsilon / (a + std::sqrt(discriminant));

        if (b > 1e-12) {
            if (a > 1e-14) {
                const double arg = std::sqrt(b * result.delta_eps / (2.0 * a));
                result.phi_tail = (2.0 / std::sqrt(R_E * b)) * std::asinh(arg);
            } else {
                // a ~ 0: acosh formulune dusus
                result.phi_tail = (1.0 / std::sqrt(R_E * b)) * std::acosh(1.0 + epsilon / R_E);
            }
        } else {
            // Dejenere: a ~ 0, b ~ 0
            const double a_safe = std::max(a, 1e-8);
            result.delta_eps = epsilon / a_safe;
            result.phi_tail  = (1.0 / a_safe) * std::sqrt(2.0 * epsilon / R_E);
        }
    }

    return result;
}

} // anonymous namespace

// ============================================================================
// GeodesicSolver::solve — Tam devre geodesic yol hesabi
// ============================================================================
std::optional<GeodesicPath> GeodesicSolver::solve(const GeodesicParams& params)
{
    // --- Girdi dogrulama ---
    if (params.R_E >= params.R_eq) return std::nullopt;
    if (params.R_E <= 0.0) return std::nullopt;
    if (params.R_eq <= params.r0) return std::nullopt;
    if (params.r0 <= 0.0) return std::nullopt;

    // --- Mandrel olustur ---
    geometry::MandrelGeometry mandrel(
        params.dome_type,
        params.R_eq,
        params.r0,
        params.L_cyl,
        params.k,
        params.N_dome_points
    );

    const auto& dome_table = mandrel.domeTable();
    const auto& meta = dome_table.metadata();
    const double s_dome = meta.s_total;
    const double R_E = params.R_E;
    const double epsilon = params.epsilon;

    // --- s_epsilon tespiti: rho(s_eps) = R_E + epsilon ---
    auto inv_eps = dome_table.inverseLookup(R_E + epsilon);
    if (!inv_eps) return std::nullopt;
    const double s_epsilon = inv_eps->s;

    // --- Numerik ODE integrasyonu: s in [0, s_epsilon] ---
    constexpr double REL_TOL = 1e-10;
    constexpr double ABS_TOL = 1e-12;

    auto ode_sol = solvePhiOde(dome_table, R_E, s_epsilon, REL_TOL, ABS_TOL);
    const double phi_numerical = ode_sol.phi.back();

    // --- Analitik kuyruk ---
    auto tail = computeAnalyticTail(dome_table, R_E, epsilon, s_epsilon, s_dome);

    // --- Toplam dome phi ---
    const double phi_dome = phi_numerical + tail.phi_tail;

    // --- Silindir katkisi ---
    const double alpha_eq = std::asin(R_E / params.R_eq);
    const double phi_cyl = params.L_cyl * std::tan(alpha_eq) / params.R_eq;

    // --- Tam devre ---
    const double delta_phi = 4.0 * phi_dome + 2.0 * phi_cyl;

    // --- Cikti struct ---
    GeodesicPath path;
    path.delta_phi    = delta_phi;
    path.phi_dome     = phi_dome;
    path.phi_cyl      = phi_cyl;
    path.phi_tail     = tail.phi_tail;
    path.phi_numerical = phi_numerical;
    path.alpha_eq     = alpha_eq;

    // --- Yol vektorlerini olustur (sadece dome parcasi) ---
    // ODE cozum noktalari: dome-lokal s [0, s_epsilon]
    // Clairaut: alpha(s) = asin(R_E / rho(s))
    const std::size_t n_ode = ode_sol.s.size();
    path.s.reserve(n_ode);
    path.rho.reserve(n_ode);
    path.x.reserve(n_ode);
    path.phi.reserve(n_ode);
    path.alpha.reserve(n_ode);

    for (std::size_t i = 0; i < n_ode; ++i) {
        auto pt = dome_table.query(ode_sol.s[i]);
        if (!pt) continue;

        path.s.push_back(ode_sol.s[i]);
        path.rho.push_back(pt->rho);
        path.x.push_back(pt->x_local);
        path.phi.push_back(ode_sol.phi[i]);

        const double sin_alpha = R_E / pt->rho;
        const double clamped = std::min(sin_alpha, 1.0);
        path.alpha.push_back(std::asin(clamped));
    }

    return path;
}

// ============================================================================
// PCHIP interpolasyon altyapisi — Fritsch-Carlson egim kestirimi
// ============================================================================
namespace {

// Fritsch-Carlson monoton egim kestirimi
// Girdi: x (kesinlikle artan), y (fonksiyon degerleri)
// Cikti: d (turevler, y ile ayni boyut)
std::vector<double> fritschCarlsonSlopes(
    const std::vector<double>& x,
    const std::vector<double>& y)
{
    const std::size_t n = x.size();
    std::vector<double> d(n, 0.0);

    if (n < 2) return d;
    if (n == 2) {
        const double slope = (y[1] - y[0]) / (x[1] - x[0]);
        d[0] = slope;
        d[1] = slope;
        return d;
    }

    // Secant egimleri
    std::vector<double> delta(n - 1);
    std::vector<double> h(n - 1);
    for (std::size_t i = 0; i < n - 1; ++i) {
        h[i] = x[i + 1] - x[i];
        delta[i] = (y[i + 1] - y[i]) / h[i];
    }

    // Ic noktalar: harmonik ortalama (Fritsch-Carlson)
    for (std::size_t i = 1; i < n - 1; ++i) {
        if (delta[i - 1] * delta[i] > 0.0) {
            // Ayni isaretli: agirlikli harmonik ortalama
            const double w1 = 2.0 * h[i] + h[i - 1];
            const double w2 = h[i] + 2.0 * h[i - 1];
            d[i] = (w1 + w2) / (w1 / delta[i - 1] + w2 / delta[i]);
        } else {
            // Farkli isaretli veya biri sifir: turev = 0
            d[i] = 0.0;
        }
    }

    // Uc noktalar: tek tarafli formul + monotonluk kontrolu
    // Sol uc
    d[0] = ((2.0 * h[0] + h[1]) * delta[0] - h[0] * delta[1]) / (h[0] + h[1]);
    if (d[0] * delta[0] < 0.0) {
        d[0] = 0.0;
    } else if (delta[0] * delta[1] < 0.0 && std::abs(d[0]) > 3.0 * std::abs(delta[0])) {
        d[0] = 3.0 * delta[0];
    }

    // Sag uc
    const std::size_t m = n - 2;
    d[n - 1] = ((2.0 * h[m] + h[m - 1]) * delta[m] - h[m] * delta[m - 1]) / (h[m] + h[m - 1]);
    if (d[n - 1] * delta[m] < 0.0) {
        d[n - 1] = 0.0;
    } else if (delta[m] * delta[m - 1] < 0.0 && std::abs(d[n - 1]) > 3.0 * std::abs(delta[m])) {
        d[n - 1] = 3.0 * delta[m];
    }

    return d;
}

// Tek nokta PCHIP degerlendirmesi: x'te kesinlikle artan, y ve d eslesmeli
// xi: sorgu noktasi (x[0] <= xi <= x[n-1])
double pchipEval(
    const std::vector<double>& x,
    const std::vector<double>& y,
    const std::vector<double>& d,
    double xi)
{
    const std::size_t n = x.size();

    // Binary search ile aralik bul
    auto it = std::upper_bound(x.begin(), x.end(), xi);
    std::size_t j = (it == x.begin()) ? 0
                  : static_cast<std::size_t>(it - x.begin()) - 1;
    if (j >= n - 1) j = n - 2;

    const double h = x[j + 1] - x[j];
    const double t = (xi - x[j]) / h;
    const double t2 = t * t;
    const double t3 = t2 * t;

    // Hermite baz fonksiyonlari
    const double h00 = 2.0 * t3 - 3.0 * t2 + 1.0;
    const double h10 = t3 - 2.0 * t2 + t;
    const double h01 = -2.0 * t3 + 3.0 * t2;
    const double h11 = t3 - t2;

    return h00 * y[j] + h10 * h * d[j] + h01 * y[j + 1] + h11 * h * d[j + 1];
}

// Vektor boyunca PCHIP yeniden ornekleme
std::vector<double> pchipResampleChannel(
    const std::vector<double>& s_orig,
    const std::vector<double>& y_orig,
    const std::vector<double>& s_new)
{
    auto slopes = fritschCarlsonSlopes(s_orig, y_orig);
    std::vector<double> y_new(s_new.size());
    for (std::size_t i = 0; i < s_new.size(); ++i) {
        y_new[i] = pchipEval(s_orig, y_orig, slopes, s_new[i]);
    }
    return y_new;
}

} // anonymous namespace

// ============================================================================
// resample — PCHIP ile esit aralikli yeniden ornekleme
// ============================================================================
GeodesicPath resample(const GeodesicPath& path, double delta_s)
{
    if (path.s.size() < 2 || delta_s <= 0.0) return path;

    const double s_begin = path.s.front();
    const double s_end   = path.s.back();
    const double range    = s_end - s_begin;

    if (range <= 0.0) return path;

    // Esit aralikli s_new olustur
    const std::size_t n_new = static_cast<std::size_t>(std::floor(range / delta_s)) + 1;
    std::vector<double> s_new(n_new);
    for (std::size_t i = 0; i < n_new; ++i) {
        s_new[i] = s_begin + static_cast<double>(i) * delta_s;
    }
    // Son noktayi ekle (eger delta_s tam bolmuyorsa)
    if (s_new.back() < s_end - 1e-12) {
        s_new.push_back(s_end);
    }

    // Her kanali PCHIP ile yeniden ornekle
    GeodesicPath resampled;
    resampled.s     = s_new;
    resampled.rho   = pchipResampleChannel(path.s, path.rho, s_new);
    resampled.x     = pchipResampleChannel(path.s, path.x, s_new);
    resampled.phi   = pchipResampleChannel(path.s, path.phi, s_new);
    resampled.alpha = pchipResampleChannel(path.s, path.alpha, s_new);

    // Skaler alanlari kopyala
    resampled.delta_phi    = path.delta_phi;
    resampled.phi_dome     = path.phi_dome;
    resampled.phi_cyl      = path.phi_cyl;
    resampled.phi_tail     = path.phi_tail;
    resampled.phi_numerical = path.phi_numerical;
    resampled.alpha_eq     = path.alpha_eq;

    return resampled;
}

} // namespace geodesic
} // namespace filament
