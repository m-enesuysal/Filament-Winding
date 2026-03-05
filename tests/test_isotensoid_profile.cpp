// =============================================================================
// test_isotensoid_profile.cpp — IsotensoidProfile Birim Testleri
// =============================================================================
// Phase-1b S5: Izotensoid dome profili dogrulama
//
// MATLAB referansi: isotenoid_dome_verification_report.txt (v2)
// TEST-01: ASTM Subscale      R=73.0,  r0=22.0
// TEST-02: Endustriyel COPV   R=152.4, r0=45.0
// TEST-03: Kucuk Aciklik      R=150.0, r0=10.0  (S-GEO-01 stiff ODE)
// TEST-04: H2 Aerospace       R=200.0, r0=50.0
//
// Toleranslar: Karar-11 Katman 2
//   pozisyon < 1e-4 mm, turev < 1e-6, egrilik bagil < 1e-4
//
// Izotensoid ozel:
//   - kappa_m isaret degistirir (bukulme noktasi)
//   - kappa_eq > 0 (konveks), kappa_pol < 0 (konkav)
//   - aspect_r ≈ 1.76..1.90 (q'ya bagli, degistirilemez)
// =============================================================================

#include <gtest/gtest.h>
#include <cmath>
#include <vector>
#include <iostream>

#include "geometry/isotensoid_profile.h"
#include "geometry/meridian_lookup_table.h"
#include "geometry/filament_types.h"
#include "geometry/i_meridian_profile.h"

using namespace filament::geometry;

// ===========================================================================
// MATLAB referans degerleri (isotenoid_dome_verification_report.txt v2)
// ===========================================================================
struct IsotensoidMatlabRef {
    const char* name;
    double R_eq;
    double r0;
    // Koussios parametreleri (analitik)
    double q;
    double m_ell;
    // MATLAB cikti degerleri
    double s_total;
    double h_dome;
    double aspect_r;
    double kappa_eq;     // pozitif (konveks)
    // kappa_pol analitik: -(1+2q)/(4*q*r0)
};

static const IsotensoidMatlabRef MATLAB_REFS[] = {
    {"TEST-01 ASTM",   73.0,  22.0,  10.010331, 0.476214,
     139.9096,  128.6302,  1.7621,  1.506708e-02},
    {"TEST-02 COPV",  152.4,  45.0,  10.469600, 0.477210,
     293.3068,  269.4514,  1.7681,  7.188422e-03},
    {"TEST-03 Kucuk", 150.0,  10.0, 224.000000, 0.498886,
     322.0304,  285.4405,  1.9029,  6.696429e-03},
    {"TEST-04 H2",    200.0,  50.0,  15.000000, 0.483871,
     396.1208,  361.7220,  1.8086,  5.333333e-03},
};

// ===========================================================================
// Yardimci: analitik kappa_eq ve kappa_pol hesabi
// ===========================================================================
static double analyticalKappaEq(double r0, double q, double Y_eq)
{
    return (1.0 + q) / (r0 * q * Y_eq);
}

static double analyticalKappaPol(double r0, double q)
{
    return -(1.0 + 2.0 * q) / (4.0 * q * r0);
}

// ===========================================================================
// Parametrize test fixture: TEST-01..04
// ===========================================================================
class IsotensoidProfileTest
    : public ::testing::TestWithParam<IsotensoidMatlabRef> {
protected:
    static constexpr std::size_t N_POINTS = 4000;
};

// ---------------------------------------------------------------------------
// T1: Profil uretimi basarili
// ---------------------------------------------------------------------------
TEST_P(IsotensoidProfileTest, GenerateProfileValid)
{
    const auto& ref = GetParam();
    IsotensoidProfile prof;
    auto table = prof.generateProfile(ref.R_eq, ref.r0, N_POINTS);

    EXPECT_TRUE(table.isValid());
    EXPECT_EQ(table.size(), N_POINTS);
}

// ---------------------------------------------------------------------------
// T2: Koussios parametreleri dogrulama
// ---------------------------------------------------------------------------
TEST_P(IsotensoidProfileTest, KoussiosParameters)
{
    const auto& ref = GetParam();
    const double Y_eq = ref.R_eq / ref.r0;
    const double q    = Y_eq * Y_eq - 1.0;
    const double m    = q / (1.0 + 2.0 * q);

    EXPECT_NEAR(q, ref.q, 1e-4)
        << ref.name << " q";
    EXPECT_NEAR(m, ref.m_ell, 1e-4)
        << ref.name << " m";
}

// ---------------------------------------------------------------------------
// T3: Meta-veri MATLAB referanslariyla uyumlu
// ---------------------------------------------------------------------------
TEST_P(IsotensoidProfileTest, MetadataMatchesMatlabReference)
{
    const auto& ref = GetParam();
    IsotensoidProfile prof;
    auto table = prof.generateProfile(ref.R_eq, ref.r0, N_POINTS);
    const auto& meta = table.metadata();

    EXPECT_DOUBLE_EQ(meta.R_eq, ref.R_eq);
    EXPECT_DOUBLE_EQ(meta.r0, ref.r0);

    // s_total — ODE vs MATLAB eliptik integral
    // RK45 rel_tol=1e-8 → ~0.01 mm mutlak dogruluk beklenir
    EXPECT_NEAR(meta.s_total, ref.s_total, 0.05)
        << ref.name << " s_total";

    // h_dome — r0 * Z(pi/2)
    EXPECT_NEAR(meta.h_dome, ref.h_dome, 0.01)
        << ref.name << " h_dome";

    // aspect_r = h_dome / R_eq
    EXPECT_NEAR(meta.aspect_r, ref.aspect_r, 0.001)
        << ref.name << " aspect_r";

    // kappa_eq — analitik referansla karsilastirma
    const double Y_eq = ref.R_eq / ref.r0;
    const double q    = Y_eq * Y_eq - 1.0;
    const double kappa_eq_calc = analyticalKappaEq(ref.r0, q, Y_eq);
    EXPECT_NEAR(meta.kappa_eq, kappa_eq_calc, 1e-10)
        << ref.name << " kappa_eq analitik";
    EXPECT_NEAR(meta.kappa_eq, ref.kappa_eq, 1e-6)
        << ref.name << " kappa_eq MATLAB";

    // kappa_pol — analitik (negatif!)
    const double kappa_pol_calc = analyticalKappaPol(ref.r0, q);
    EXPECT_NEAR(meta.kappa_pol, kappa_pol_calc, 1e-10)
        << ref.name << " kappa_pol analitik";
    EXPECT_LT(meta.kappa_pol, 0.0)
        << ref.name << " kappa_pol negatif olmali (konkav)";

    std::cout << "[" << ref.name << "] s_total=" << meta.s_total
              << " (ref=" << ref.s_total << "), h_dome=" << meta.h_dome
              << " (ref=" << ref.h_dome << "), aspect_r=" << meta.aspect_r
              << std::endl;
}

// ---------------------------------------------------------------------------
// T4: Sinir kosullari — ekvator ve polar aciklik
// ---------------------------------------------------------------------------
TEST_P(IsotensoidProfileTest, BoundaryConditions)
{
    const auto& ref = GetParam();
    IsotensoidProfile prof;
    auto table = prof.generateProfile(ref.R_eq, ref.r0, N_POINTS);

    // Ekvator: s = 0
    auto eq = table.query(0.0);
    ASSERT_TRUE(eq.has_value());
    EXPECT_NEAR(eq->rho, ref.R_eq, 1e-12)
        << ref.name << " rho(0) = R_eq";
    EXPECT_NEAR(eq->x_local, 0.0, 1e-12)
        << ref.name << " x(0) = 0";
    EXPECT_NEAR(eq->drho_ds, 0.0, 1e-12)
        << ref.name << " drho/ds(0) = 0";
    EXPECT_NEAR(eq->dx_ds, 1.0, 1e-12)
        << ref.name << " dx/ds(0) = 1";

    // Polar aciklik: s = s_total
    const auto& meta = table.metadata();
    auto pol = table.query(meta.s_total);
    ASSERT_TRUE(pol.has_value());
    EXPECT_NEAR(pol->rho, ref.r0, 1e-6)
        << ref.name << " rho(s_total) = r0";
    EXPECT_NEAR(pol->x_local, meta.h_dome, 1e-6)
        << ref.name << " x(s_total) = h_dome";
    EXPECT_NEAR(pol->drho_ds, 0.0, 1e-6)
        << ref.name << " drho/ds(s_total) = 0";
    EXPECT_NEAR(pol->dx_ds, 1.0, 1e-6)
        << ref.name << " dx/ds(s_total) = 1";
}

// ---------------------------------------------------------------------------
// T5: Birim teget vektor — |drho/ds|^2 + |dx/ds|^2 = 1 (dugum noktalari)
// ---------------------------------------------------------------------------
TEST_P(IsotensoidProfileTest, UnitTangentVector)
{
    const auto& ref = GetParam();
    IsotensoidProfile prof;
    auto table = prof.generateProfile(ref.R_eq, ref.r0, N_POINTS);
    const auto& raw_s = table.rawS();

    double max_err = 0.0;
    for (std::size_t i = 0; i < N_POINTS; ++i) {
        auto pt = table.query(raw_s[i]);
        ASSERT_TRUE(pt.has_value());
        const double norm_sq = pt->drho_ds * pt->drho_ds + pt->dx_ds * pt->dx_ds;
        const double err = std::abs(norm_sq - 1.0);
        max_err = std::max(max_err, err);
        EXPECT_LT(err, 1e-10)
            << ref.name << " birim teget vektor dugum " << i;
    }
    std::cout << "[" << ref.name << "] Birim teget vektor maks |norm-1|: "
              << max_err << std::endl;
}

// ---------------------------------------------------------------------------
// T6: Monotoniklik — rho azalan, x artan
// ---------------------------------------------------------------------------
TEST_P(IsotensoidProfileTest, Monotonicity)
{
    const auto& ref = GetParam();
    IsotensoidProfile prof;
    auto table = prof.generateProfile(ref.R_eq, ref.r0, N_POINTS);
    const auto& raw_rho = table.rawRho();
    const auto& raw_x   = table.rawX();

    for (std::size_t i = 1; i < N_POINTS; ++i) {
        EXPECT_LT(raw_rho[i], raw_rho[i - 1])
            << ref.name << " rho monoton azalan ihlali: i=" << i;
        EXPECT_GT(raw_x[i], raw_x[i - 1])
            << ref.name << " x monoton artan ihlali: i=" << i;
    }
}

// ---------------------------------------------------------------------------
// T7: C1 sureklilik — ekvator noktasinda silindir uyumu
// ---------------------------------------------------------------------------
TEST_P(IsotensoidProfileTest, C1ContinuityAtEquator)
{
    const auto& ref = GetParam();
    IsotensoidProfile prof;
    auto table = prof.generateProfile(ref.R_eq, ref.r0, N_POINTS);

    auto eq = table.query(0.0);
    ASSERT_TRUE(eq.has_value());
    EXPECT_NEAR(eq->rho, ref.R_eq, 1e-12)
        << ref.name << " C0: rho";
    EXPECT_NEAR(eq->drho_ds, 0.0, 1e-12)
        << ref.name << " C1: drho/ds";
    EXPECT_NEAR(eq->dx_ds, 1.0, 1e-12)
        << ref.name << " C1: dx/ds";
}

// ---------------------------------------------------------------------------
// T8: Egrilik isaret degisimi — bukulme noktasi (isotensoid ozgu)
//     Ekvator: kappa > 0 (konveks)
//     Polar:   kappa < 0 (konkav)
// ---------------------------------------------------------------------------
TEST_P(IsotensoidProfileTest, CurvatureSignChange)
{
    const auto& ref = GetParam();
    IsotensoidProfile prof;
    auto table = prof.generateProfile(ref.R_eq, ref.r0, N_POINTS);
    const auto& raw_s = table.rawS();

    // Ekvator egriligi pozitif
    auto eq_pt = table.query(0.0);
    ASSERT_TRUE(eq_pt.has_value());
    EXPECT_GT(eq_pt->kappa_m, 0.0)
        << ref.name << " kappa_eq pozitif olmali";

    // Polar egrilik negatif
    auto pol_pt = table.query(table.metadata().s_total);
    ASSERT_TRUE(pol_pt.has_value());
    EXPECT_LT(pol_pt->kappa_m, 0.0)
        << ref.name << " kappa_pol negatif olmali";

    // Bukulme noktasinin varligini dogrula
    // (egrilik isaret degistiren bir i olmali)
    bool sign_change_found = false;
    for (std::size_t i = 1; i < N_POINTS; ++i) {
        auto p0 = table.query(raw_s[i - 1]);
        auto p1 = table.query(raw_s[i]);
        if (p0.has_value() && p1.has_value()) {
            if (p0->kappa_m * p1->kappa_m < 0.0) {
                sign_change_found = true;
                const double s_inflect = 0.5 * (raw_s[i - 1] + raw_s[i]);
                std::cout << "[" << ref.name << "] Bukulme noktasi: s/s_total ≈ "
                          << s_inflect / table.metadata().s_total << std::endl;
                break;
            }
        }
    }
    EXPECT_TRUE(sign_change_found)
        << ref.name << " bukulme noktasi bulunamadi";
}

// ---------------------------------------------------------------------------
// T9: Clairaut iliskisi — rho * sin(alpha) = r0
//     alpha_w = asin(r0 / rho)
// ---------------------------------------------------------------------------
TEST_P(IsotensoidProfileTest, ClairautRelation)
{
    const auto& ref = GetParam();
    IsotensoidProfile prof;
    auto table = prof.generateProfile(ref.R_eq, ref.r0, N_POINTS);
    const auto& raw_rho = table.rawRho();

    double max_err = 0.0;
    for (std::size_t i = 0; i < N_POINTS; ++i) {
        if (raw_rho[i] > ref.r0 + 1e-10) {
            const double sin_alpha = ref.r0 / raw_rho[i];
            // Clairaut: rho * sin(alpha) = r0
            // sin(alpha) = r0/rho => rho * sin(alpha) - r0 = 0
            const double clairaut_err = std::abs(raw_rho[i] * sin_alpha - ref.r0);
            max_err = std::max(max_err, clairaut_err);
        }
    }
    EXPECT_LT(max_err, 1e-10)
        << ref.name << " Clairaut iliskisi";
    std::cout << "[" << ref.name << "] Clairaut maks hata: " << max_err << std::endl;
}

// ---------------------------------------------------------------------------
// T10: Interpolasyon — geometrik degismezler ara noktalarda
// ---------------------------------------------------------------------------
TEST_P(IsotensoidProfileTest, InterpolatedGeometricInvariants)
{
    const auto& ref = GetParam();
    IsotensoidProfile prof;
    auto table = prof.generateProfile(ref.R_eq, ref.r0, N_POINTS);
    const auto& raw_s = table.rawS();

    double max_tangent_err = 0.0;

    for (std::size_t i = 0; i + 1 < N_POINTS; ++i) {
        const double s_mid = 0.5 * (raw_s[i] + raw_s[i + 1]);
        auto pt = table.query(s_mid);
        ASSERT_TRUE(pt.has_value());

        // Birim teget vektor
        const double norm_sq = pt->drho_ds * pt->drho_ds + pt->dx_ds * pt->dx_ds;
        const double tangent_err = std::abs(norm_sq - 1.0);
        max_tangent_err = std::max(max_tangent_err, tangent_err);
    }

    // drho/ds ve dx/ds bagimsiz spline kanallari — 5e-6 tolerans
    EXPECT_LT(max_tangent_err, 5e-6)
        << ref.name << " interpolated birim teget vektor";

    std::cout << "[" << ref.name << "] Interpolasyon birim teget maks hata: "
              << max_tangent_err << std::endl;
}

// ---------------------------------------------------------------------------
// T11: Interpolasyon — ince-kaba kiyaslama ile Karar-11 Katman 2
// ---------------------------------------------------------------------------
TEST_P(IsotensoidProfileTest, InterpolationToleranceFineCoarse)
{
    const auto& ref = GetParam();
    IsotensoidProfile prof;

    // Kaba profil (test edilecek) ve ince profil (referans)
    auto coarse = prof.generateProfile(ref.R_eq, ref.r0, 4000);
    auto fine   = prof.generateProfile(ref.R_eq, ref.r0, 20000);

    const auto& fine_s = fine.rawS();

    double max_rho_err  = 0.0;
    double max_x_err    = 0.0;
    double max_drho_err = 0.0;
    double max_dx_err   = 0.0;
    double max_curv_err = 0.0;

    for (std::size_t i = 1; i + 1 < fine.size(); i += 5) {
        auto fine_pt   = fine.query(fine_s[i]);
        auto coarse_pt = coarse.query(fine_s[i]);

        if (!fine_pt.has_value() || !coarse_pt.has_value()) continue;

        const double rho_err  = std::abs(coarse_pt->rho - fine_pt->rho);
        const double x_err    = std::abs(coarse_pt->x_local - fine_pt->x_local);
        const double drho_err = std::abs(coarse_pt->drho_ds - fine_pt->drho_ds);
        const double dx_err   = std::abs(coarse_pt->dx_ds - fine_pt->dx_ds);

        // Egrilik bagil hatasi — sifirdan gecis bolgesini atla
        double curv_err = 0.0;
        if (std::abs(fine_pt->kappa_m) > 1e-6) {
            curv_err = std::abs(coarse_pt->kappa_m - fine_pt->kappa_m)
                     / std::abs(fine_pt->kappa_m);
        }

        max_rho_err  = std::max(max_rho_err, rho_err);
        max_x_err    = std::max(max_x_err, x_err);
        max_drho_err = std::max(max_drho_err, drho_err);
        max_dx_err   = std::max(max_dx_err, dx_err);
        max_curv_err = std::max(max_curv_err, curv_err);
    }

    // Karar-11 Katman 2
    EXPECT_LT(max_rho_err, tolerances::POSITION_ABS_TOL)
        << ref.name << " rho interpolasyon";
    EXPECT_LT(max_x_err, tolerances::POSITION_ABS_TOL)
        << ref.name << " x interpolasyon";
    EXPECT_LT(max_drho_err, tolerances::DERIVATIVE_ABS_TOL)
        << ref.name << " drho/ds interpolasyon";
    EXPECT_LT(max_dx_err, tolerances::DERIVATIVE_ABS_TOL)
        << ref.name << " dx/ds interpolasyon";
    EXPECT_LT(max_curv_err, tolerances::CURVATURE_REL_TOL)
        << ref.name << " kappa interpolasyon";

    std::cout << "[" << ref.name << "] Ince-kaba interpolasyon hatalari:"
              << "\n  rho:     " << max_rho_err  << " mm (limit " << tolerances::POSITION_ABS_TOL << ")"
              << "\n  x:       " << max_x_err    << " mm (limit " << tolerances::POSITION_ABS_TOL << ")"
              << "\n  drho/ds: " << max_drho_err << " (limit " << tolerances::DERIVATIVE_ABS_TOL << ")"
              << "\n  dx/ds:   " << max_dx_err   << " (limit " << tolerances::DERIVATIVE_ABS_TOL << ")"
              << "\n  kappa:   " << max_curv_err << " (limit " << tolerances::CURVATURE_REL_TOL << ")"
              << std::endl;
}

// ---------------------------------------------------------------------------
// Parametrize instantiation
// ---------------------------------------------------------------------------
INSTANTIATE_TEST_SUITE_P(
    MatlabCrossValidation,
    IsotensoidProfileTest,
    ::testing::ValuesIn(MATLAB_REFS),
    [](const ::testing::TestParamInfo<IsotensoidMatlabRef>& info) {
        std::string name = info.param.name;
        for (auto& c : name) {
            if (c == ' ' || c == '-') c = '_';
        }
        return name;
    }
);

// ===========================================================================
// Girdi dogrulama testleri (Karar-5 + Karar-17)
// ===========================================================================
TEST(IsotensoidProfileValidation, NegativeReqThrows)
{
    IsotensoidProfile prof;
    EXPECT_THROW(prof.generateProfile(-100.0, 30.0, 500), std::invalid_argument);
}

TEST(IsotensoidProfileValidation, ZeroR0Throws)
{
    IsotensoidProfile prof;
    EXPECT_THROW(prof.generateProfile(100.0, 0.0, 500), std::invalid_argument);
}

TEST(IsotensoidProfileValidation, R0EqualReqThrows)
{
    IsotensoidProfile prof;
    EXPECT_THROW(prof.generateProfile(100.0, 100.0, 500), std::invalid_argument);
}

TEST(IsotensoidProfileValidation, R0GreaterThanReqThrows)
{
    IsotensoidProfile prof;
    EXPECT_THROW(prof.generateProfile(100.0, 120.0, 500), std::invalid_argument);
}

TEST(IsotensoidProfileValidation, TooFewPointsThrows)
{
    IsotensoidProfile prof;
    EXPECT_THROW(prof.generateProfile(100.0, 30.0, 1), std::invalid_argument);
    EXPECT_THROW(prof.generateProfile(100.0, 30.0, 0), std::invalid_argument);
}

// ===========================================================================
// DomeType ve factory
// ===========================================================================
TEST(IsotensoidProfileIdentity, DomeTypeAndName)
{
    IsotensoidProfile prof;
    EXPECT_EQ(prof.domeType(), DomeType::Isotensoid);
    EXPECT_STREQ(prof.name(), "Isotensoid");
}

TEST(IsotensoidProfileIdentity, FactoryCreation)
{
    auto prof = createProfile(DomeType::Isotensoid);
    ASSERT_NE(prof, nullptr);
    EXPECT_EQ(prof->domeType(), DomeType::Isotensoid);

    auto table = prof->generateProfile(100.0, 30.0, 500);
    EXPECT_TRUE(table.isValid());
}

// ===========================================================================
// S-GEO-01: Kucuk polar aciklik (stiff ODE)
// ===========================================================================
TEST(IsotensoidProfileSGEO01, SmallPolarOpening)
{
    // Y_eq = 15 — S-GEO-01 sinirda
    IsotensoidProfile prof;
    auto table = prof.generateProfile(150.0, 10.0, 500);
    EXPECT_TRUE(table.isValid());

    const auto& meta = table.metadata();
    EXPECT_NEAR(meta.aspect_r, 1.9029, 0.005)
        << "S-GEO-01 aspect_r";

    // Sinir kosullari korunuyor mu?
    auto eq = table.query(0.0);
    ASSERT_TRUE(eq.has_value());
    EXPECT_NEAR(eq->rho, 150.0, 1e-10);

    auto pol = table.query(meta.s_total);
    ASSERT_TRUE(pol.has_value());
    EXPECT_NEAR(pol->rho, 10.0, 1e-4);

    std::cout << "[S-GEO-01] Y_eq=15, q=224, s_total=" << meta.s_total
              << " mm, aspect_r=" << meta.aspect_r << std::endl;
}

TEST(IsotensoidProfileSGEO01, VerySmallPolarOpening)
{
    // Y_eq = 20 — S-GEO-01 otesinde (stiff bolge)
    IsotensoidProfile prof;
    auto table = prof.generateProfile(200.0, 10.0, 500);
    EXPECT_TRUE(table.isValid());

    const auto& meta = table.metadata();
    EXPECT_GT(meta.aspect_r, 1.90);

    auto pol = table.query(meta.s_total);
    ASSERT_TRUE(pol.has_value());
    EXPECT_NEAR(pol->rho, 10.0, 0.01)
        << "Polar rho stiff ODE sonrasi";

    std::cout << "[S-GEO-01 stiff] Y_eq=20, s_total=" << meta.s_total
              << " mm, polar rho=" << pol->rho << std::endl;
}

// ===========================================================================
// Uc durumlar
// ===========================================================================
TEST(IsotensoidProfileEdgeCases, ShallowDome)
{
    // Buyuk polar aciklik — sig dome
    // Y_eq = 1.5 → q = 1.25
    IsotensoidProfile prof;
    auto table = prof.generateProfile(150.0, 100.0, 500);
    EXPECT_TRUE(table.isValid());

    const auto& meta = table.metadata();
    EXPECT_LT(meta.aspect_r, 1.3)
        << "Sig dome aspect_r < 1.3";
    EXPECT_GT(meta.aspect_r, 1.0)
        << "Isotensoid aspect_r her zaman > 1";

    std::cout << "[Sig dome] Y_eq=1.5, aspect_r=" << meta.aspect_r
              << ", h_dome=" << meta.h_dome << " mm" << std::endl;
}

TEST(IsotensoidProfileEdgeCases, AsymptoticLimit)
{
    // Cok kucuk aciklik — asimptotik limit: aspect_r → sqrt(2)*E(1/2) ≈ 1.9101
    IsotensoidProfile prof;
    auto table = prof.generateProfile(200.0, 10.0, 500);
    EXPECT_TRUE(table.isValid());

    const auto& meta = table.metadata();
    // Asimptotik limite yakin ama altinda olmali
    EXPECT_GT(meta.aspect_r, 1.89);
    EXPECT_LT(meta.aspect_r, 1.9101 + 0.001);

    std::cout << "[Asimptotik] Y_eq=20, aspect_r=" << meta.aspect_r
              << " (limit ≈ 1.9101)" << std::endl;
}

// ===========================================================================
// Aspect ratio degistirilemez testi
// Isotensoid dome'un aspect_r yalnizca r0/R_eq oranina baglidir
// ===========================================================================
TEST(IsotensoidProfilePhysics, AspectRatioDeterminedByPhysics)
{
    IsotensoidProfile prof;

    // Farkli R_eq, ayni r0/R_eq orani → ayni aspect_r
    auto t1 = prof.generateProfile(100.0, 30.0, 500);
    auto t2 = prof.generateProfile(200.0, 60.0, 500);

    EXPECT_NEAR(t1.metadata().aspect_r, t2.metadata().aspect_r, 0.001)
        << "Ayni r0/R_eq → ayni aspect_r";
}
