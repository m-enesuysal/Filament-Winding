// =============================================================================
// isotensoid_profile.cpp — Izotensoid Dome Meridyen Profili
// =============================================================================
// Karar-10: IsotensoidProfile : IMeridianProfile
// Karar-6:  Dormand-Prince RK45 adaptif cozucu (standalone, Boost.Odeint API uyumlu)
// Karar-11: Katman 1 toleranslari — RelTol=1e-8, AbsTol=1e-10
// S-GEO-01: Polar aciklik yakininda adaptif adim kontrolu
//
// Matematiksel Model — Koussios Eliptik Integral Cozumu:
//   Y(theta) = sqrt(1 + q*sin^2(theta)),  theta in [0, pi/2]
//   theta=0: polar aciklik,  theta=pi/2: ekvator
//   q = Y_eq^2 - 1 = (R_eq/r0)^2 - 1
//   m = q/(1+2q)  (eliptik integral parametresi)
//
// ODE Sistemi (theta parametresi):
//   dY/dtheta  = q*sin(theta)*cos(theta) / Y
//   dZ/dtheta  = q*(1+cos^2(theta)) / [sqrt(1+2q)*Delta]
//   ds/dtheta  = r0*sqrt((dY/dtheta)^2 + (dZ/dtheta)^2)
//   Delta      = sqrt(1 - m*sin^2(theta))
//
// Turevler (s cinsinden, analitik):
//   drho/ds = -dY_dth / sqrt(dY_dth^2 + dZ_dth^2)
//   dx/ds   =  dZ_dth / sqrt(dY_dth^2 + dZ_dth^2)
//
// Egrilik (parametrik formul):
//   kappa_m = (dY*d2Z - dZ*d2Y) / {r0*[(dY)^2 + (dZ)^2]^(3/2)}
//   Ekvator limiti: kappa_eq  =  (1+q)/(r0*q*Y_eq)   [pozitif, konveks]
//   Polar limiti:   kappa_pol = -(1+2q)/(4*q*r0)      [negatif, konkav]
//
// Phase-1a MATLAB referansi: isotensoid_dome_profile.m (v2)
// =============================================================================

#include "geometry/isotensoid_profile.h"
#include <cmath>
#include <stdexcept>
#include <string>
#include <vector>
#include <array>
#include <algorithm>

namespace filament {
namespace geometry {

namespace {

// ============================================================================
// Dormand-Prince RK45 Standalone Adaptif ODE Cozucu
// ============================================================================
// Boost.Odeint runge_kutta_dopri5 ile ayni Butcher tablosu.
// Karar-6: RelTol=1e-8, AbsTol=1e-10
// ============================================================================

using State3 = std::array<double, 3>;  // [Y, Z, s_ode]

// -----------------------------------------------------------------------
// Butcher tablosu — Dormand & Prince (1980)
// 5. derece cozum (stepping) + 4. derece hata tahmini (FSAL)
// -----------------------------------------------------------------------
namespace dp {
    // Nodes
    constexpr double c2 = 1.0 / 5.0;
    constexpr double c3 = 3.0 / 10.0;
    constexpr double c4 = 4.0 / 5.0;
    constexpr double c5 = 8.0 / 9.0;
    // c6 = 1, c7 = 1

    // Matris katsayilari
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

    // 5. derece agirliklar (= a7x satiri)
    constexpr double b1 = 35.0 / 384.0;
    // b2 = 0
    constexpr double b3 = 500.0 / 1113.0;
    constexpr double b4 = 125.0 / 192.0;
    constexpr double b5 = -2187.0 / 6784.0;
    constexpr double b6 = 11.0 / 84.0;
    // b7 = 0

    // Hata katsayilari (b - b_hat)
    constexpr double e1 = 71.0 / 57600.0;
    // e2 = 0
    constexpr double e3 = -71.0 / 16695.0;
    constexpr double e4 = 71.0 / 1920.0;
    constexpr double e5 = -17253.0 / 339200.0;
    constexpr double e6 = 22.0 / 525.0;
    constexpr double e7 = -1.0 / 40.0;
} // namespace dp

// -----------------------------------------------------------------------
// Izotensoid ODE sag taraf fonksiyonu
// -----------------------------------------------------------------------
struct IsotensoidRHS {
    double q;
    double m_ell;
    double sqrt_1p2q;
    double r0;

    State3 operator()(double theta, const State3& y) const
    {
        const double st  = std::sin(theta);
        const double ct  = std::cos(theta);
        const double st2 = st * st;
        const double ct2 = ct * ct;
        const double Y   = y[0];

        // dY/dtheta
        const double dY = q * st * ct / Y;

        // dZ/dtheta
        const double Delta = std::sqrt(1.0 - m_ell * st2);
        const double dZ = q * (1.0 + ct2) / (sqrt_1p2q * Delta);

        // ds/dtheta
        const double ds = r0 * std::sqrt(dY * dY + dZ * dZ);

        return {dY, dZ, ds};
    }
};

// -----------------------------------------------------------------------
// ODE cozum konteyneri
// -----------------------------------------------------------------------
struct OdeSolution {
    std::vector<double> theta;
    std::vector<double> Y;
    std::vector<double> Z;
    std::vector<double> s_ode;
};

// -----------------------------------------------------------------------
// Dormand-Prince RK45 adaptif integrator
// -----------------------------------------------------------------------
OdeSolution solveDormandPrince(
    const IsotensoidRHS& rhs,
    double theta_end,
    double rel_tol,
    double abs_tol)
{
    constexpr double SAFETY    = 0.9;
    constexpr double H_MIN     = 1e-15;
    constexpr double GROW_MAX  = 5.0;
    constexpr double SHRINK_MAX = 0.2;
    constexpr int    MAX_STEPS = 200000;

    const double h_max = theta_end * 0.005;  // maks adim = %0.5 aralik (ODE interpolasyon dogrulugu)

    double theta = 0.0;
    State3 y = {1.0, 0.0, 0.0};  // Y(0)=1, Z(0)=0, s(0)=0
    double h = theta_end / 200.0;

    OdeSolution sol;
    sol.theta.reserve(5000);
    sol.Y.reserve(5000);
    sol.Z.reserve(5000);
    sol.s_ode.reserve(5000);

    sol.theta.push_back(theta);
    sol.Y.push_back(y[0]);
    sol.Z.push_back(y[1]);
    sol.s_ode.push_back(y[2]);

    for (int step = 0; step < MAX_STEPS && theta < theta_end - H_MIN; ++step) {
        // Son noktayi asma
        if (theta + h > theta_end) {
            h = theta_end - theta;
        }
        if (h < H_MIN) h = H_MIN;

        // --- 7 asama Dormand-Prince ---
        const State3 k1 = rhs(theta, y);

        State3 yt;
        for (int i = 0; i < 3; ++i)
            yt[i] = y[i] + h * dp::a21 * k1[i];
        const State3 k2 = rhs(theta + dp::c2 * h, yt);

        for (int i = 0; i < 3; ++i)
            yt[i] = y[i] + h * (dp::a31 * k1[i] + dp::a32 * k2[i]);
        const State3 k3 = rhs(theta + dp::c3 * h, yt);

        for (int i = 0; i < 3; ++i)
            yt[i] = y[i] + h * (dp::a41 * k1[i] + dp::a42 * k2[i]
                               + dp::a43 * k3[i]);
        const State3 k4 = rhs(theta + dp::c4 * h, yt);

        for (int i = 0; i < 3; ++i)
            yt[i] = y[i] + h * (dp::a51 * k1[i] + dp::a52 * k2[i]
                               + dp::a53 * k3[i] + dp::a54 * k4[i]);
        const State3 k5 = rhs(theta + dp::c5 * h, yt);

        for (int i = 0; i < 3; ++i)
            yt[i] = y[i] + h * (dp::a61 * k1[i] + dp::a62 * k2[i]
                               + dp::a63 * k3[i] + dp::a64 * k4[i]
                               + dp::a65 * k5[i]);
        const State3 k6 = rhs(theta + h, yt);

        // 5. derece cozum
        State3 y_new;
        for (int i = 0; i < 3; ++i)
            y_new[i] = y[i] + h * (dp::b1 * k1[i] + dp::b3 * k3[i]
                                  + dp::b4 * k4[i] + dp::b5 * k5[i]
                                  + dp::b6 * k6[i]);

        // k7 (FSAL — hata tahmini icin)
        const State3 k7 = rhs(theta + h, y_new);

        // --- Hata tahmini ---
        double err = 0.0;
        for (int i = 0; i < 3; ++i) {
            const double ei = h * (dp::e1 * k1[i] + dp::e3 * k3[i]
                                 + dp::e4 * k4[i] + dp::e5 * k5[i]
                                 + dp::e6 * k6[i] + dp::e7 * k7[i]);
            const double sc = abs_tol
                + rel_tol * std::max(std::abs(y[i]), std::abs(y_new[i]));
            const double ratio = std::abs(ei) / sc;
            if (ratio > err) err = ratio;
        }

        if (err <= 1.0) {
            // Adim kabul
            theta += h;
            y = y_new;
            sol.theta.push_back(theta);
            sol.Y.push_back(y[0]);
            sol.Z.push_back(y[1]);
            sol.s_ode.push_back(y[2]);
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

// -----------------------------------------------------------------------
// Kubik Hermite interpolasyon yardimcisi
// x: kesinlikle artan, y: fonksiyon degerleri, dy: turevler (dy/dx)
// O(h^4) dogruluk — lineer interpolasyona gore cok daha iyi
// -----------------------------------------------------------------------
double hermiteLookup(
    const std::vector<double>& x,
    const std::vector<double>& y,
    const std::vector<double>& dy,
    double xi)
{
    auto it = std::upper_bound(x.begin(), x.end(), xi);
    if (it == x.begin()) return y.front();
    if (it == x.end())   return y.back();

    const auto j = static_cast<std::size_t>(it - x.begin());
    const double h = x[j] - x[j - 1];
    const double t = (xi - x[j - 1]) / h;
    const double t2 = t * t;
    const double t3 = t2 * t;

    // Hermite baz fonksiyonlari
    const double h00 =  2.0 * t3 - 3.0 * t2 + 1.0;
    const double h10 =        t3 - 2.0 * t2 + t;
    const double h01 = -2.0 * t3 + 3.0 * t2;
    const double h11 =        t3 -       t2;

    return h00 * y[j - 1] + h10 * h * dy[j - 1]
         + h01 * y[j]     + h11 * h * dy[j];
}

} // anonymous namespace

// ===========================================================================
// IsotensoidProfile::generateProfile
// ===========================================================================
MeridianLookupTable IsotensoidProfile::generateProfile(
    double R_eq, double r0, std::size_t N_points) const
{
    // -------------------------------------------------------------------
    // Girdi dogrulama (Karar-5 + Karar-17: ust katman exception)
    // -------------------------------------------------------------------
    if (R_eq <= 0.0) {
        throw std::invalid_argument(
            "IsotensoidProfile: R_eq pozitif olmali. R_eq = "
            + std::to_string(R_eq));
    }
    if (r0 <= 0.0) {
        throw std::invalid_argument(
            "IsotensoidProfile: r0 pozitif olmali. r0 = "
            + std::to_string(r0));
    }
    if (r0 >= R_eq) {
        throw std::invalid_argument(
            "IsotensoidProfile: r0 < R_eq olmali. r0 = "
            + std::to_string(r0) + ", R_eq = " + std::to_string(R_eq));
    }
    if (N_points < 2) {
        throw std::invalid_argument(
            "IsotensoidProfile: N_points en az 2 olmali. N_points = "
            + std::to_string(N_points));
    }

    // -------------------------------------------------------------------
    // Koussios parametreleri
    // -------------------------------------------------------------------
    const double Y_eq     = R_eq / r0;
    const double q        = Y_eq * Y_eq - 1.0;
    const double m_ell    = q / (1.0 + 2.0 * q);
    const double sqrt_1p2q = std::sqrt(1.0 + 2.0 * q);

    // -------------------------------------------------------------------
    // Faz 1: Dormand-Prince RK45 ODE integrasyonu
    //        theta: 0 (polar) → pi/2 (ekvator)
    //        Durum: [Y, Z, s_ode]
    //        Karar-6: rel_tol=1e-8, abs_tol=1e-10
    // -------------------------------------------------------------------
    IsotensoidRHS rhs_func{q, m_ell, sqrt_1p2q, r0};

    OdeSolution ode = solveDormandPrince(
        rhs_func,
        constants::PI / 2.0,
        tolerances::ODE_REL_TOL,
        tolerances::ODE_ABS_TOL);

    // Son noktayi tam pi/2 olarak sabitle
    ode.theta.back() = constants::PI / 2.0;
    ode.Y.back()     = Y_eq;  // analitik: Y(pi/2) = sqrt(1+q) = Y_eq

    const double s_total = ode.s_ode.back();
    const double h_dome  = r0 * ode.Z.back();

    // -------------------------------------------------------------------
    // Faz 1b: Hermite interpolasyon icin turevleri hesapla
    //         ODE cikti noktalarinda dZ/dtheta ve ds/dtheta
    // -------------------------------------------------------------------
    const std::size_t n_ode = ode.theta.size();
    std::vector<double> ode_dZ_dth(n_ode);
    std::vector<double> ode_ds_dth(n_ode);

    for (std::size_t k = 0; k < n_ode; ++k) {
        State3 state_k = {ode.Y[k], ode.Z[k], ode.s_ode[k]};
        State3 rhs_k = rhs_func(ode.theta[k], state_k);
        ode_dZ_dth[k] = rhs_k[1];  // dZ/dtheta
        ode_ds_dth[k] = rhs_k[2];  // ds_ode/dtheta
    }

    // -------------------------------------------------------------------
    // Faz 2: Uniform theta izgarasi + analitik hesap
    //        theta_i: pi/2 (ekvator, i=0) → 0 (polar, i=N-1)
    //        s_i:     0    (ekvator)       → s_total (polar)
    // -------------------------------------------------------------------
    std::vector<double> s(N_points);
    std::vector<double> rho(N_points);
    std::vector<double> x_local(N_points);
    std::vector<double> drho_ds(N_points);
    std::vector<double> dx_ds(N_points);
    std::vector<double> kappa_m(N_points);

    for (std::size_t i = 0; i < N_points; ++i) {
        // theta: pi/2 → 0 (azalan)
        const double th = (constants::PI / 2.0)
            * static_cast<double>(N_points - 1 - i)
            / static_cast<double>(N_points - 1);

        const double st  = std::sin(th);
        const double ct  = std::cos(th);
        const double st2 = st * st;
        const double ct2 = ct * ct;

        // Y analitik
        const double Y = std::sqrt(1.0 + q * st2);

        // Z ve s_ode: ODE cozumunden kubik Hermite interpolasyon
        // ODE theta kesinlikle artan (0 → pi/2)
        const double Z_val     = hermiteLookup(ode.theta, ode.Z, ode_dZ_dth, th);
        const double s_ode_val = hermiteLookup(ode.theta, ode.s_ode, ode_ds_dth, th);

        // Fiziksel koordinatlar
        rho[i]     = r0 * Y;
        x_local[i] = h_dome - r0 * Z_val;

        // s konvansiyonu: s=0 ekvator (theta=pi/2), s=s_total polar (theta=0)
        s[i] = s_total - s_ode_val;

        // --- Turevler (analitik, theta'dan) ---
        const double dY_dth = q * st * ct / Y;
        const double Delta  = std::sqrt(1.0 - m_ell * st2);
        const double dZ_dth = q * (1.0 + ct2) / (sqrt_1p2q * Delta);
        const double mag2   = dY_dth * dY_dth + dZ_dth * dZ_dth;
        const double mag    = std::sqrt(mag2);

        if (mag > 1e-30) {
            // ds_our/dtheta = -ds_ode/dtheta  (konvansiyon flip)
            // drho/ds = (r0*dY_dth) / (ds_our/dtheta)
            //         = (r0*dY_dth) / (-r0*mag) = -dY_dth/mag
            drho_ds[i] = -dY_dth / mag;
            // dx/dtheta = d(h_dome - r0*Z)/dtheta = -r0*dZ_dth
            // dx/ds = (-r0*dZ_dth) / (-r0*mag) = dZ_dth/mag
            dx_ds[i]   =  dZ_dth / mag;
        } else {
            drho_ds[i] = 0.0;
            dx_ds[i]   = 1.0;
        }

        // --- Egrilik (parametrik formul) ---
        // d2Y/dtheta2
        const double Y3 = Y * Y * Y;
        const double d2Y = q * std::cos(2.0 * th) / Y
                         - q * q * st2 * ct2 / Y3;

        // d2Z/dtheta2 (v2 duzeltilmis — MATLAB referansi)
        const double Delta3 = Delta * Delta * Delta;
        const double pow_1p2q_32 = sqrt_1p2q * (1.0 + 2.0 * q);  // (1+2q)^(3/2)
        const double d2Z = -q * st * ct * (2.0 * (1.0 + q) - q * st2)
                         / (pow_1p2q_32 * Delta3);

        const double numer = dY_dth * d2Z - dZ_dth * d2Y;
        const double denom = mag2 * mag;  // (dY^2 + dZ^2)^(3/2)

        if (denom > 1e-60) {
            kappa_m[i] = numer / (r0 * denom);
        } else {
            kappa_m[i] = 0.0;
        }
    }

    // -------------------------------------------------------------------
    // Sinir degerlerini analitik limitlerle sabitle
    // -------------------------------------------------------------------
    // Ekvator (i=0, theta=pi/2)
    s[0]        = 0.0;
    rho[0]      = R_eq;
    x_local[0]  = 0.0;
    drho_ds[0]  = 0.0;
    dx_ds[0]    = 1.0;
    kappa_m[0]  = (1.0 + q) / (r0 * q * Y_eq);

    // Polar aciklik (i=N-1, theta=0)
    s[N_points - 1]       = s_total;
    rho[N_points - 1]     = r0;
    x_local[N_points - 1] = h_dome;
    drho_ds[N_points - 1] = 0.0;
    dx_ds[N_points - 1]   = 1.0;
    kappa_m[N_points - 1] = -(1.0 + 2.0 * q) / (4.0 * q * r0);

    // -------------------------------------------------------------------
    // s monotoniklik guvenlik kontrolu
    // theta azaldikca s_our artmali (ekvator→polar)
    // -------------------------------------------------------------------
    for (std::size_t i = 1; i < N_points; ++i) {
        if (s[i] <= s[i - 1]) {
            // Numerik yuvarlama durumunda kucuk duzeltme
            s[i] = s[i - 1] + 1e-15;
        }
    }

    // -------------------------------------------------------------------
    // Meta-veri hesabi
    // -------------------------------------------------------------------

    // Yuzey alani (trapezoidal kural)
    double A_dome = 0.0;
    for (std::size_t i = 0; i + 1 < N_points; ++i) {
        A_dome += 0.5 * (rho[i] + rho[i + 1]) * (s[i + 1] - s[i]);
    }
    A_dome *= 2.0 * constants::PI;

    ProfileMetadata meta;
    meta.R_eq      = R_eq;
    meta.r0        = r0;
    meta.s_total   = s_total;
    meta.h_dome    = h_dome;
    meta.A_dome    = A_dome;
    meta.kappa_eq  = kappa_m[0];
    meta.kappa_pol = kappa_m[N_points - 1];
    meta.alpha_w   = std::asin(r0 / R_eq);
    meta.aspect_r  = h_dome / R_eq;

    // -------------------------------------------------------------------
    // Tablo olustur
    // -------------------------------------------------------------------
    MeridianLookupTable table;
    table.build(s, rho, x_local, drho_ds, dx_ds, kappa_m, meta);

    return table;
}

} // namespace geometry
} // namespace filament
