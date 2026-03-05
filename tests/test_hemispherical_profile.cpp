// =============================================================================
// test_hemispherical_profile.cpp — HemisphericalProfile Birim Testleri
// =============================================================================
// Phase-1b S3: Kapali-form hemispherical dome profili dogrulama
//
// MATLAB referansi: hemispherical_dome_verification_report.txt
// TEST-01: ASTM Subscale      R_eq=73.0,  r0=22.0
// TEST-02: Endustriyel COPV   R_eq=152.4, r0=45.0
// TEST-03: Kucuk Aciklik      R_eq=150.0, r0=10.0
// TEST-04: H2 Aerospace       R_eq=200.0, r0=50.0
//
// Toleranslar: Karar-11 Katman 2
//   Pozisyon  |eps| < 1e-4 mm   (mutlak)
//   Turev     |eps| < 1e-6      (mutlak, boyutsuz)
//   Egrilik   |eps_rel| < 1e-4  (bagil, %0.01)
// =============================================================================

#include <gtest/gtest.h>
#include <cmath>
#include <vector>
#include <iostream>

#include "geometry/hemispherical_profile.h"
#include "geometry/meridian_lookup_table.h"
#include "geometry/filament_types.h"
#include "geometry/i_meridian_profile.h"

using namespace filament::geometry;

// ===========================================================================
// MATLAB referans degerleri (hemispherical_dome_verification_report.txt)
// ===========================================================================
struct MatlabReference {
    const char* name;
    double R_eq;
    double r0;

    // MATLAB cikti degerleri
    double r0_over_R;   // r0/R_eq
    double theta_p_deg; // polar aci [derece]
    double s_total;     // toplam yay uzunlugu [mm]
    double h_dome;      // dome yuksekligi [mm]
    double A_dome;      // dome yuzey alani [mm^2]
};

// MATLAB raporundan alinmis tam referans degerleri
static const MatlabReference MATLAB_REFS[] = {
    {"TEST-01 ASTM Subscale",     73.0,  22.0,  0.3014,  72.46,   92.3207,  69.6060,  31926.38},
    {"TEST-02 Endustriyel COPV", 152.4,  45.0,  0.2953,  72.83,  193.7084, 145.6048, 139424.97},
    {"TEST-03 Kucuk Aciklik",    150.0,  10.0,  0.0667,  86.18,  225.6120, 149.6663, 141057.16},
    {"TEST-04 H2 Aerospace",     200.0,  50.0,  0.2500,  75.52,  263.6232, 193.6492, 243346.72},
};

// ===========================================================================
// Analitik yardimci fonksiyon
// ===========================================================================
static MeridianPoint hemisphericalAnalytic(double s, double R_eq)
{
    const double theta = s / R_eq;
    MeridianPoint pt;
    pt.s       = s;
    pt.rho     = R_eq * std::cos(theta);
    pt.x_local = R_eq * std::sin(theta);
    pt.drho_ds = -std::sin(theta);
    pt.dx_ds   =  std::cos(theta);
    pt.kappa_m = 1.0 / R_eq;
    return pt;
}

// ===========================================================================
// Parametrize test fixture: TEST-01..04
// ===========================================================================
class HemisphericalProfileTest
    : public ::testing::TestWithParam<MatlabReference> {
protected:
    static constexpr std::size_t N_POINTS = 500;
    HemisphericalProfile profile_;
};

// ---------------------------------------------------------------------------
// T1: Profil uretimi basarili + tablo gecerli
// ---------------------------------------------------------------------------
TEST_P(HemisphericalProfileTest, GenerateProfileValid)
{
    const auto& ref = GetParam();
    auto table = profile_.generateProfile(ref.R_eq, ref.r0, N_POINTS);

    EXPECT_TRUE(table.isValid());
    EXPECT_EQ(table.size(), N_POINTS);
}

// ---------------------------------------------------------------------------
// T2: Meta-veri MATLAB referanslariyla uyumlu
// ---------------------------------------------------------------------------
TEST_P(HemisphericalProfileTest, MetadataMatchesMatlabReference)
{
    const auto& ref = GetParam();
    auto table = profile_.generateProfile(ref.R_eq, ref.r0, N_POINTS);
    const auto& meta = table.metadata();

    // Temel parametreler
    EXPECT_DOUBLE_EQ(meta.R_eq, ref.R_eq);
    EXPECT_DOUBLE_EQ(meta.r0, ref.r0);

    // MATLAB referans karsilastirmasi (4 ondalik hassasiyet, MATLAB rapor formati)
    EXPECT_NEAR(meta.s_total, ref.s_total, 0.001)
        << ref.name << " s_total";
    EXPECT_NEAR(meta.h_dome, ref.h_dome, 0.001)
        << ref.name << " h_dome";
    EXPECT_NEAR(meta.A_dome, ref.A_dome, 0.1)
        << ref.name << " A_dome";

    // Egrilik: hemispherical dome icin sabit 1/R_eq
    EXPECT_DOUBLE_EQ(meta.kappa_eq, 1.0 / ref.R_eq);
    EXPECT_DOUBLE_EQ(meta.kappa_pol, 1.0 / ref.R_eq);

    // Aspect ratio: hemispherical = 1.0
    EXPECT_DOUBLE_EQ(meta.aspect_r, 1.0);

    // Winding angle: alpha_w = asin(r0/R_eq)
    const double expected_alpha = std::asin(ref.r0 / ref.R_eq);
    EXPECT_NEAR(meta.alpha_w, expected_alpha, 1e-14);

    // r0/R_eq orani
    EXPECT_NEAR(ref.r0 / ref.R_eq, ref.r0_over_R, 0.001);

    // theta_p derece
    const double theta_p_deg = std::acos(ref.r0 / ref.R_eq) * 180.0 / constants::PI;
    EXPECT_NEAR(theta_p_deg, ref.theta_p_deg, 0.01)
        << ref.name << " theta_p derece";
}

// ---------------------------------------------------------------------------
// T3: Sinir kosullari — ekvator ve polar aciklik
// ---------------------------------------------------------------------------
TEST_P(HemisphericalProfileTest, BoundaryConditions)
{
    const auto& ref = GetParam();
    auto table = profile_.generateProfile(ref.R_eq, ref.r0, N_POINTS);

    // Ekvator sinir kosullari: s = 0
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

    // Polar aciklik sinir kosullari: s = s_total
    const auto& meta = table.metadata();
    auto pol = table.query(meta.s_total);
    ASSERT_TRUE(pol.has_value());
    EXPECT_NEAR(pol->rho, ref.r0, 1e-8)
        << ref.name << " rho(s_total) = r0";
    EXPECT_NEAR(pol->x_local, ref.h_dome, 0.001)
        << ref.name << " x(s_total) = h_dome";
}

// ---------------------------------------------------------------------------
// T4: Kure denklemi — rho^2 + x^2 = R_eq^2 (tum dugum noktalarinda)
// ---------------------------------------------------------------------------
TEST_P(HemisphericalProfileTest, SphereEquation)
{
    const auto& ref = GetParam();
    auto table = profile_.generateProfile(ref.R_eq, ref.r0, N_POINTS);
    const auto& raw_rho = table.rawRho();
    const auto& raw_x   = table.rawX();
    const double R_sq = ref.R_eq * ref.R_eq;

    double max_err = 0.0;
    for (std::size_t i = 0; i < N_POINTS; ++i) {
        const double sphere_err = std::abs(
            raw_rho[i] * raw_rho[i] + raw_x[i] * raw_x[i] - R_sq);
        max_err = std::max(max_err, sphere_err);
        EXPECT_LT(sphere_err, 1e-8)
            << ref.name << " kure denklemi hatasi dugum " << i;
    }
    std::cout << "[" << ref.name << "] Kure denklemi maks hata: "
              << max_err << " mm^2" << std::endl;
}

// ---------------------------------------------------------------------------
// T5: Egrilik sabitligi — kappa_m = 1/R_eq (sabit)
// ---------------------------------------------------------------------------
TEST_P(HemisphericalProfileTest, ConstantCurvature)
{
    const auto& ref = GetParam();
    auto table = profile_.generateProfile(ref.R_eq, ref.r0, N_POINTS);
    const auto& raw_s = table.rawS();
    const double kappa_ref = 1.0 / ref.R_eq;

    for (std::size_t i = 0; i < N_POINTS; ++i) {
        auto pt = table.query(raw_s[i]);
        ASSERT_TRUE(pt.has_value());
        // Sabit egrilik: bagil hata = 0 (tam eslesme beklenir)
        EXPECT_NEAR(pt->kappa_m, kappa_ref, 1e-14)
            << ref.name << " kappa_m dugum " << i;
    }
}

// ---------------------------------------------------------------------------
// T6: Birim teget vektor — |drho/ds|^2 + |dx/ds|^2 = 1
// ---------------------------------------------------------------------------
TEST_P(HemisphericalProfileTest, UnitTangentVector)
{
    const auto& ref = GetParam();
    auto table = profile_.generateProfile(ref.R_eq, ref.r0, N_POINTS);
    const auto& raw_s = table.rawS();

    double max_err = 0.0;
    for (std::size_t i = 0; i < N_POINTS; ++i) {
        auto pt = table.query(raw_s[i]);
        ASSERT_TRUE(pt.has_value());
        const double norm_sq = pt->drho_ds * pt->drho_ds + pt->dx_ds * pt->dx_ds;
        const double err = std::abs(norm_sq - 1.0);
        max_err = std::max(max_err, err);
        EXPECT_LT(err, 1e-14)
            << ref.name << " birim teget vektor hatasi dugum " << i;
    }
    std::cout << "[" << ref.name << "] Birim teget vektor maks |norm-1|: "
              << max_err << std::endl;
}

// ---------------------------------------------------------------------------
// T7: Monotoniklik — rho azalan, x artan
// ---------------------------------------------------------------------------
TEST_P(HemisphericalProfileTest, Monotonicity)
{
    const auto& ref = GetParam();
    auto table = profile_.generateProfile(ref.R_eq, ref.r0, N_POINTS);
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
// T8: C1 sureklilik — ekvator noktasinda silindir uyumu
// ---------------------------------------------------------------------------
TEST_P(HemisphericalProfileTest, C1ContinuityAtEquator)
{
    const auto& ref = GetParam();
    auto table = profile_.generateProfile(ref.R_eq, ref.r0, N_POINTS);

    // s=0 noktasinda dome degerleri (silindir tarafindan bakis)
    auto eq = table.query(0.0);
    ASSERT_TRUE(eq.has_value());

    // Silindir: rho = R_eq, drho/ds = 0, dx/ds = 1, kappa_m = 0
    // Dome:     rho = R_eq, drho/ds = 0, dx/ds = 1, kappa_m = 1/R_eq
    EXPECT_NEAR(eq->rho, ref.R_eq, 1e-12)
        << ref.name << " C0: rho surekli";
    EXPECT_NEAR(eq->drho_ds, 0.0, 1e-12)
        << ref.name << " C1: drho/ds surekli";
    EXPECT_NEAR(eq->dx_ds, 1.0, 1e-12)
        << ref.name << " C1: dx/ds surekli";
    // C2 sureksizligi beklenir: kappa 0 -> 1/R_eq (Karar-10: C1 yeterli)
}

// ---------------------------------------------------------------------------
// T9: Ara noktalarda Karar-11 Katman 2 toleranslari (spline interpolasyon)
// ---------------------------------------------------------------------------
TEST_P(HemisphericalProfileTest, InterpolationTolerances)
{
    const auto& ref = GetParam();
    auto table = profile_.generateProfile(ref.R_eq, ref.r0, N_POINTS);
    const auto& raw_s = table.rawS();

    double max_rho_err  = 0.0;
    double max_x_err    = 0.0;
    double max_drho_err = 0.0;
    double max_dx_err   = 0.0;
    double max_curv_err = 0.0;
    const double kappa_ref = 1.0 / ref.R_eq;

    // Her araligin ortasinda sorgu
    for (std::size_t i = 0; i + 1 < N_POINTS; ++i) {
        const double s_mid = 0.5 * (raw_s[i] + raw_s[i + 1]);
        auto result = table.query(s_mid);
        ASSERT_TRUE(result.has_value());

        const auto analytic = hemisphericalAnalytic(s_mid, ref.R_eq);

        const double rho_err  = std::abs(result->rho - analytic.rho);
        const double x_err    = std::abs(result->x_local - analytic.x_local);
        const double drho_err = std::abs(result->drho_ds - analytic.drho_ds);
        const double dx_err   = std::abs(result->dx_ds - analytic.dx_ds);
        const double curv_err = std::abs(result->kappa_m - kappa_ref)
                              / std::abs(kappa_ref);

        max_rho_err  = std::max(max_rho_err, rho_err);
        max_x_err    = std::max(max_x_err, x_err);
        max_drho_err = std::max(max_drho_err, drho_err);
        max_dx_err   = std::max(max_dx_err, dx_err);
        max_curv_err = std::max(max_curv_err, curv_err);

        // Karar-11 Katman 2
        EXPECT_LT(rho_err, tolerances::POSITION_ABS_TOL)
            << ref.name << " rho aralik " << i;
        EXPECT_LT(x_err, tolerances::POSITION_ABS_TOL)
            << ref.name << " x aralik " << i;
        EXPECT_LT(drho_err, tolerances::DERIVATIVE_ABS_TOL)
            << ref.name << " drho/ds aralik " << i;
        EXPECT_LT(dx_err, tolerances::DERIVATIVE_ABS_TOL)
            << ref.name << " dx/ds aralik " << i;
        EXPECT_LT(curv_err, tolerances::CURVATURE_REL_TOL)
            << ref.name << " kappa aralik " << i;
    }

    std::cout << "[" << ref.name << "] Interpolasyon maks hatalari:"
              << "\n  rho:     " << max_rho_err  << " mm (limit " << tolerances::POSITION_ABS_TOL << ")"
              << "\n  x:       " << max_x_err    << " mm (limit " << tolerances::POSITION_ABS_TOL << ")"
              << "\n  drho/ds: " << max_drho_err << " (limit " << tolerances::DERIVATIVE_ABS_TOL << ")"
              << "\n  dx/ds:   " << max_dx_err   << " (limit " << tolerances::DERIVATIVE_ABS_TOL << ")"
              << "\n  kappa:   " << max_curv_err << " (limit " << tolerances::CURVATURE_REL_TOL << ")"
              << std::endl;
}

// ---------------------------------------------------------------------------
// Parametrize: 4 test senaryosu
// ---------------------------------------------------------------------------
INSTANTIATE_TEST_SUITE_P(
    MatlabCrossValidation,
    HemisphericalProfileTest,
    ::testing::ValuesIn(MATLAB_REFS),
    [](const ::testing::TestParamInfo<MatlabReference>& info) {
        // Test ismi icin _ kullan (Google Test sinirlamasi)
        std::string name = info.param.name;
        // bosluk ve tire yerine _
        for (auto& c : name) {
            if (c == ' ' || c == '-') c = '_';
        }
        return name;
    }
);

// ===========================================================================
// Girdi dogrulama testleri (Karar-5 + Karar-17 ust katman)
// ===========================================================================
TEST(HemisphericalProfileValidation, NegativeReqThrows)
{
    HemisphericalProfile prof;
    EXPECT_THROW(prof.generateProfile(-100.0, 30.0, 500), std::invalid_argument);
}

TEST(HemisphericalProfileValidation, ZeroReqThrows)
{
    HemisphericalProfile prof;
    EXPECT_THROW(prof.generateProfile(0.0, 30.0, 500), std::invalid_argument);
}

TEST(HemisphericalProfileValidation, NegativeR0Throws)
{
    HemisphericalProfile prof;
    EXPECT_THROW(prof.generateProfile(100.0, -5.0, 500), std::invalid_argument);
}

TEST(HemisphericalProfileValidation, ZeroR0Throws)
{
    HemisphericalProfile prof;
    EXPECT_THROW(prof.generateProfile(100.0, 0.0, 500), std::invalid_argument);
}

TEST(HemisphericalProfileValidation, R0EqualReqThrows)
{
    HemisphericalProfile prof;
    EXPECT_THROW(prof.generateProfile(100.0, 100.0, 500), std::invalid_argument);
}

TEST(HemisphericalProfileValidation, R0GreaterThanReqThrows)
{
    HemisphericalProfile prof;
    EXPECT_THROW(prof.generateProfile(100.0, 120.0, 500), std::invalid_argument);
}

TEST(HemisphericalProfileValidation, TooFewPointsThrows)
{
    HemisphericalProfile prof;
    EXPECT_THROW(prof.generateProfile(100.0, 30.0, 1), std::invalid_argument);
    EXPECT_THROW(prof.generateProfile(100.0, 30.0, 0), std::invalid_argument);
}

TEST(HemisphericalProfileValidation, MinimalPointsSucceeds)
{
    HemisphericalProfile prof;
    auto table = prof.generateProfile(100.0, 30.0, 2);
    EXPECT_TRUE(table.isValid());
    EXPECT_EQ(table.size(), 2u);
}

// ===========================================================================
// Factory fonksiyonu testi
// ===========================================================================
TEST(HemisphericalProfileFactory, CreateViaFactory)
{
    auto prof = createProfile(DomeType::Hemispherical);
    ASSERT_NE(prof, nullptr);
    EXPECT_EQ(prof->domeType(), DomeType::Hemispherical);
    EXPECT_STREQ(prof->name(), "Hemispherical");

    auto table = prof->generateProfile(100.0, 30.0, 500);
    EXPECT_TRUE(table.isValid());
    EXPECT_EQ(table.size(), 500u);
}

// ===========================================================================
// DomeType ve isim dogrulama
// ===========================================================================
TEST(HemisphericalProfileIdentity, DomeTypeAndName)
{
    HemisphericalProfile prof;
    EXPECT_EQ(prof.domeType(), DomeType::Hemispherical);
    EXPECT_STREQ(prof.name(), "Hemispherical");
}

// ===========================================================================
// Uc durum: sik dome (r0/R_eq ~ 0.95) ve derin dome (r0 kucuk)
// ===========================================================================
TEST(HemisphericalProfileEdgeCases, ShallowDome)
{
    // r0/R_eq = 0.95 — cok sik dome (theta_p kucuk)
    HemisphericalProfile prof;
    auto table = prof.generateProfile(152.4, 0.95 * 152.4, 500);
    EXPECT_TRUE(table.isValid());

    const auto& meta = table.metadata();
    EXPECT_GT(meta.s_total, 0.0);
    EXPECT_GT(meta.h_dome, 0.0);
    EXPECT_LT(meta.h_dome, meta.R_eq);  // sik dome: h < R

    // Ekvator sinirlari
    auto eq = table.query(0.0);
    ASSERT_TRUE(eq.has_value());
    EXPECT_NEAR(eq->rho, 152.4, 1e-10);
}

TEST(HemisphericalProfileEdgeCases, DeepDome)
{
    // r0 = 5mm, R_eq = 152.4mm — derin dome (theta_p ~ pi/2)
    HemisphericalProfile prof;
    auto table = prof.generateProfile(152.4, 5.0, 500);
    EXPECT_TRUE(table.isValid());

    const auto& meta = table.metadata();
    // Derin dome: h ~ R_eq
    EXPECT_GT(meta.h_dome, 0.99 * meta.R_eq);

    // Polar aciklik noktasi
    auto pol = table.query(meta.s_total);
    ASSERT_TRUE(pol.has_value());
    EXPECT_NEAR(pol->rho, 5.0, 1e-8);
}
