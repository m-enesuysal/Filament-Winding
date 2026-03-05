// =============================================================================
// test_ellipsoidal_profile.cpp — EllipsoidalProfile Birim Testleri
// =============================================================================
// Phase-1b S4: Elipsoidal dome profili dogrulama
//
// MATLAB referansi: ellipsoidal_dome_verification_report.txt
// TEST-01: ASTM Subscale      R=73.0,  r0=22.0, k=0.60
// TEST-02: Endustriyel COPV   R=152.4, r0=45.0, k=0.70
// TEST-03: Kucuk Aciklik      R=150.0, r0=10.0, k=0.50
// TEST-04: H2 Aerospace       R=200.0, r0=50.0, k=0.70
//
// S-GEO-04: k=1 → hemispherical ortusme dogrulamasi
// S-GEO-03: k_min = 0.15 sinir dogrulamasi
//
// Toleranslar: Karar-11 Katman 2
// =============================================================================

#include <gtest/gtest.h>
#include <cmath>
#include <vector>
#include <iostream>

#include "geometry/ellipsoidal_profile.h"
#include "geometry/hemispherical_profile.h"
#include "geometry/meridian_lookup_table.h"
#include "geometry/filament_types.h"
#include "geometry/i_meridian_profile.h"

using namespace filament::geometry;

// ===========================================================================
// MATLAB referans degerleri (ellipsoidal_dome_verification_report.txt)
// ===========================================================================
struct EllipsoidalMatlabRef {
    const char* name;
    double R_eq;
    double r0;
    double k;
    double theta_p_deg;
    double s_total;
    double h_dome;
    double kappa_eq;
    double kappa_pol;
};

static const EllipsoidalMatlabRef MATLAB_REFS[] = {
    {"TEST-01 ASTM",    73.0,  22.0, 0.60, 72.46,  71.0473,  41.7636, 3.805175e-02, 8.991663e-03},
    {"TEST-02 COPV",   152.4,  45.0, 0.70, 72.83, 159.7323, 101.9234, 1.339118e-02, 4.917492e-03},
    {"TEST-03 Kucuk",  150.0,  10.0, 0.50, 86.18, 171.6565,  74.8331, 2.666667e-02, 3.350070e-03},
    {"TEST-04 H2",     200.0,  50.0, 0.70, 75.52, 218.8545, 135.5544, 1.020408e-02, 3.674269e-03},
};

// ===========================================================================
// Parametrize test fixture: TEST-01..04
// ===========================================================================
class EllipsoidalProfileTest
    : public ::testing::TestWithParam<EllipsoidalMatlabRef> {
protected:
    static constexpr std::size_t N_POINTS = 500;
};

// ---------------------------------------------------------------------------
// T1: Profil uretimi basarili
// ---------------------------------------------------------------------------
TEST_P(EllipsoidalProfileTest, GenerateProfileValid)
{
    const auto& ref = GetParam();
    EllipsoidalProfile prof(ref.k);
    auto table = prof.generateProfile(ref.R_eq, ref.r0, N_POINTS);

    EXPECT_TRUE(table.isValid());
    EXPECT_EQ(table.size(), N_POINTS);
}

// ---------------------------------------------------------------------------
// T2: Meta-veri MATLAB referanslariyla uyumlu
// ---------------------------------------------------------------------------
TEST_P(EllipsoidalProfileTest, MetadataMatchesMatlabReference)
{
    const auto& ref = GetParam();
    EllipsoidalProfile prof(ref.k);
    auto table = prof.generateProfile(ref.R_eq, ref.r0, N_POINTS);
    const auto& meta = table.metadata();

    EXPECT_DOUBLE_EQ(meta.R_eq, ref.R_eq);
    EXPECT_DOUBLE_EQ(meta.r0, ref.r0);
    EXPECT_DOUBLE_EQ(meta.aspect_r, ref.k);

    // MATLAB rapor degerleriyle karsilastirma
    EXPECT_NEAR(meta.s_total, ref.s_total, 0.01)
        << ref.name << " s_total";
    EXPECT_NEAR(meta.h_dome, ref.h_dome, 0.001)
        << ref.name << " h_dome";

    // Egrilik: kappa_eq = 1/(R_eq * k^2)
    const double kappa_eq_expected = 1.0 / (ref.R_eq * ref.k * ref.k);
    EXPECT_NEAR(meta.kappa_eq, kappa_eq_expected, 1e-10)
        << ref.name << " kappa_eq analitik";
    EXPECT_NEAR(meta.kappa_eq, ref.kappa_eq, 1e-6)
        << ref.name << " kappa_eq MATLAB";
    EXPECT_NEAR(meta.kappa_pol, ref.kappa_pol, 1e-6)
        << ref.name << " kappa_pol MATLAB";

    // theta_p derece
    const double theta_p_deg = std::acos(ref.r0 / ref.R_eq) * 180.0 / constants::PI;
    EXPECT_NEAR(theta_p_deg, ref.theta_p_deg, 0.01)
        << ref.name << " theta_p";

    std::cout << "[" << ref.name << "] s_total=" << meta.s_total
              << " (ref=" << ref.s_total << "), h_dome=" << meta.h_dome
              << " (ref=" << ref.h_dome << ")" << std::endl;
}

// ---------------------------------------------------------------------------
// T3: Sinir kosullari — ekvator ve polar aciklik
// ---------------------------------------------------------------------------
TEST_P(EllipsoidalProfileTest, BoundaryConditions)
{
    const auto& ref = GetParam();
    EllipsoidalProfile prof(ref.k);
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
    EXPECT_NEAR(pol->rho, ref.r0, 1e-8)
        << ref.name << " rho(s_total) = r0";
    EXPECT_NEAR(pol->x_local, ref.h_dome, 0.001)
        << ref.name << " x(s_total) = h_dome";
}

// ---------------------------------------------------------------------------
// T4: Elips denklemi — (rho/a)^2 + (x/b)^2 = 1 (tum dugum noktalarinda)
// ---------------------------------------------------------------------------
TEST_P(EllipsoidalProfileTest, EllipseEquation)
{
    const auto& ref = GetParam();
    EllipsoidalProfile prof(ref.k);
    auto table = prof.generateProfile(ref.R_eq, ref.r0, N_POINTS);
    const auto& raw_rho = table.rawRho();
    const auto& raw_x   = table.rawX();
    const double a = ref.R_eq;
    const double b = ref.k * ref.R_eq;

    double max_err = 0.0;
    for (std::size_t i = 0; i < N_POINTS; ++i) {
        const double ratio = (raw_rho[i] / a) * (raw_rho[i] / a)
                           + (raw_x[i] / b) * (raw_x[i] / b);
        const double err = std::abs(ratio - 1.0);
        max_err = std::max(max_err, err);
        EXPECT_LT(err, 1e-13)
            << ref.name << " elips denklemi dugum " << i;
    }
    std::cout << "[" << ref.name << "] Elips denklemi maks hata: "
              << max_err << std::endl;
}

// ---------------------------------------------------------------------------
// T5: Birim teget vektor — |drho/ds|^2 + |dx/ds|^2 = 1
// ---------------------------------------------------------------------------
TEST_P(EllipsoidalProfileTest, UnitTangentVector)
{
    const auto& ref = GetParam();
    EllipsoidalProfile prof(ref.k);
    auto table = prof.generateProfile(ref.R_eq, ref.r0, N_POINTS);
    const auto& raw_s = table.rawS();

    double max_err = 0.0;
    for (std::size_t i = 0; i < N_POINTS; ++i) {
        auto pt = table.query(raw_s[i]);
        ASSERT_TRUE(pt.has_value());
        const double norm_sq = pt->drho_ds * pt->drho_ds + pt->dx_ds * pt->dx_ds;
        const double err = std::abs(norm_sq - 1.0);
        max_err = std::max(max_err, err);
        EXPECT_LT(err, 1e-13)
            << ref.name << " birim teget vektoru dugum " << i;
    }
    std::cout << "[" << ref.name << "] Birim teget vektor maks |norm-1|: "
              << max_err << std::endl;
}

// ---------------------------------------------------------------------------
// T6: Monotoniklik — rho azalan, x artan
// ---------------------------------------------------------------------------
TEST_P(EllipsoidalProfileTest, Monotonicity)
{
    const auto& ref = GetParam();
    EllipsoidalProfile prof(ref.k);
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
TEST_P(EllipsoidalProfileTest, C1ContinuityAtEquator)
{
    const auto& ref = GetParam();
    EllipsoidalProfile prof(ref.k);
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
// T8: Interpolasyon — geometrik degismezler ara noktalarda
//     Elips denklemi ve birim teget vektor korunuyor mu?
// ---------------------------------------------------------------------------
TEST_P(EllipsoidalProfileTest, InterpolatedGeometricInvariants)
{
    const auto& ref = GetParam();
    EllipsoidalProfile prof(ref.k);
    auto table = prof.generateProfile(ref.R_eq, ref.r0, N_POINTS);
    const auto& raw_s = table.rawS();
    const double a = ref.R_eq;
    const double b = ref.k * ref.R_eq;

    double max_ellipse_err = 0.0;
    double max_tangent_err = 0.0;

    for (std::size_t i = 0; i + 1 < N_POINTS; ++i) {
        const double s_mid = 0.5 * (raw_s[i] + raw_s[i + 1]);
        auto pt = table.query(s_mid);
        ASSERT_TRUE(pt.has_value());

        // Elips denklemi
        const double ratio = (pt->rho / a) * (pt->rho / a)
                           + (pt->x_local / b) * (pt->x_local / b);
        const double ellipse_err = std::abs(ratio - 1.0);
        max_ellipse_err = std::max(max_ellipse_err, ellipse_err);

        // Birim teget vektor
        const double norm_sq = pt->drho_ds * pt->drho_ds + pt->dx_ds * pt->dx_ds;
        const double tangent_err = std::abs(norm_sq - 1.0);
        max_tangent_err = std::max(max_tangent_err, tangent_err);
    }

    // Elips hatasi < 1e-6 (rho ve x interpolasyon dogrulugunun bir olcusu)
    EXPECT_LT(max_ellipse_err, 1e-6)
        << ref.name << " interpolated elips denklemi";

    // Teget vektor hatasi < 5e-6
    // drho/ds ve dx/ds bagimsiz spline kanallari olarak interpolasyon yapilir,
    // bu nedenle turetilmis kisit |t|^2=1 ara noktalarda tam korunmaz.
    EXPECT_LT(max_tangent_err, 5e-6)
        << ref.name << " interpolated birim teget vektor";

    std::cout << "[" << ref.name << "] Interpolasyon geometrik degismezler:"
              << "\n  Elips denklemi maks hata:  " << max_ellipse_err
              << "\n  Birim teget maks hata:     " << max_tangent_err
              << std::endl;
}

// ---------------------------------------------------------------------------
// T9: Interpolasyon — ince-kaba kiyaslama ile Karar-11 Katman 2
// ---------------------------------------------------------------------------
TEST_P(EllipsoidalProfileTest, InterpolationToleranceFineCoarse)
{
    const auto& ref = GetParam();
    EllipsoidalProfile prof(ref.k);

    // Kaba profil (test edilecek) ve ince profil (referans)
    // 1000 nokta: k=0.5 gibi yassi dome'larda turev toleransini saglar
    auto coarse = prof.generateProfile(ref.R_eq, ref.r0, 1000);
    auto fine   = prof.generateProfile(ref.R_eq, ref.r0, 5000);

    const auto& fine_s   = fine.rawS();

    double max_rho_err  = 0.0;
    double max_x_err    = 0.0;
    double max_drho_err = 0.0;
    double max_dx_err   = 0.0;
    double max_curv_err = 0.0;

    // Ince gridin ic noktalarinda kaba tabloyu sorgula
    for (std::size_t i = 1; i + 1 < fine.size(); i += 5) {
        auto fine_pt   = fine.query(fine_s[i]);   // dugum noktasi: tam eslesme
        auto coarse_pt = coarse.query(fine_s[i]); // interpolasyon

        if (!fine_pt.has_value() || !coarse_pt.has_value()) continue;

        const double rho_err = std::abs(coarse_pt->rho - fine_pt->rho);
        const double x_err   = std::abs(coarse_pt->x_local - fine_pt->x_local);
        const double drho_err = std::abs(coarse_pt->drho_ds - fine_pt->drho_ds);
        const double dx_err   = std::abs(coarse_pt->dx_ds - fine_pt->dx_ds);

        // Egrilik bagil hatasi (sifir bolunmesinden kacinma)
        double curv_err = 0.0;
        if (std::abs(fine_pt->kappa_m) > 1e-15) {
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
    EllipsoidalProfileTest,
    ::testing::ValuesIn(MATLAB_REFS),
    [](const ::testing::TestParamInfo<EllipsoidalMatlabRef>& info) {
        std::string name = info.param.name;
        for (auto& c : name) {
            if (c == ' ' || c == '-') c = '_';
        }
        return name;
    }
);

// ===========================================================================
// S-GEO-04: k=1 Hemispherical Ortusme Dogrulamasi
// ===========================================================================
class SGEO04Test : public ::testing::Test {
protected:
    static constexpr std::size_t N = 500;
    static constexpr double R_eq = 152.4;
    static constexpr double r0   = 45.0;
};

TEST_F(SGEO04Test, ScalarMetadataMatch)
{
    EllipsoidalProfile  ellip(1.0);
    HemisphericalProfile hemi;

    auto table_e = ellip.generateProfile(R_eq, r0, N);
    auto table_h = hemi.generateProfile(R_eq, r0, N);

    const auto& me = table_e.metadata();
    const auto& mh = table_h.metadata();

    // s_total tam eslesmeli (k=1 → f(theta)=1 → integral analitige esit)
    EXPECT_NEAR(me.s_total, mh.s_total, 1e-10)
        << "s_total ortusme";
    EXPECT_NEAR(me.h_dome, mh.h_dome, 1e-10)
        << "h_dome ortusme";
    // A_dome: elipsoidal trapezoidal kural vs hemispherical analitik formula
    // N=500 icin trapezoidal hata ~O(h^2) ≈ 0.1 mm²; 139000 mm² yuzey icin <1e-6 bagil
    EXPECT_NEAR(me.A_dome, mh.A_dome, 1.0)
        << "A_dome ortusme";
    EXPECT_NEAR(me.kappa_eq, mh.kappa_eq, 1e-14)
        << "kappa_eq ortusme";
    EXPECT_NEAR(me.kappa_pol, mh.kappa_pol, 1e-14)
        << "kappa_pol ortusme";
    EXPECT_NEAR(me.alpha_w, mh.alpha_w, 1e-14)
        << "alpha_w ortusme";
}

TEST_F(SGEO04Test, PointByPointOverlap)
{
    EllipsoidalProfile  ellip(1.0);
    HemisphericalProfile hemi;

    auto table_e = ellip.generateProfile(R_eq, r0, N);
    auto table_h = hemi.generateProfile(R_eq, r0, N);

    const auto& s_e = table_e.rawS();
    const auto& s_h = table_h.rawS();

    double max_rho_err = 0.0, max_x_err = 0.0;
    double max_drho_err = 0.0, max_dx_err = 0.0;
    double max_kappa_err = 0.0;

    for (std::size_t i = 0; i < N; ++i) {
        // Arclength eslesmeli
        EXPECT_NEAR(s_e[i], s_h[i], 1e-10)
            << "s ortusme nokta " << i;

        auto pe = table_e.query(s_e[i]);
        auto ph = table_h.query(s_h[i]);
        ASSERT_TRUE(pe.has_value() && ph.has_value());

        const double rho_err = std::abs(pe->rho - ph->rho);
        const double x_err   = std::abs(pe->x_local - ph->x_local);
        const double drho_err = std::abs(pe->drho_ds - ph->drho_ds);
        const double dx_err   = std::abs(pe->dx_ds - ph->dx_ds);
        const double kappa_rel = std::abs(pe->kappa_m - ph->kappa_m)
                               / std::abs(ph->kappa_m);

        max_rho_err   = std::max(max_rho_err, rho_err);
        max_x_err     = std::max(max_x_err, x_err);
        max_drho_err  = std::max(max_drho_err, drho_err);
        max_dx_err    = std::max(max_dx_err, dx_err);
        max_kappa_err = std::max(max_kappa_err, kappa_rel);
    }

    // k=1 durumunda cok siki toleranslar
    EXPECT_LT(max_rho_err, 1e-10)   << "rho S-GEO-04";
    EXPECT_LT(max_x_err, 1e-10)     << "x S-GEO-04";
    EXPECT_LT(max_drho_err, 1e-12)  << "drho/ds S-GEO-04";
    EXPECT_LT(max_dx_err, 1e-12)    << "dx/ds S-GEO-04";
    EXPECT_LT(max_kappa_err, 1e-12) << "kappa S-GEO-04";

    std::cout << "[S-GEO-04] k=1 ortusme hatalari:"
              << "\n  |Delta rho|      = " << max_rho_err << " mm"
              << "\n  |Delta x|        = " << max_x_err << " mm"
              << "\n  |Delta drho/ds|  = " << max_drho_err
              << "\n  |Delta dx/ds|    = " << max_dx_err
              << "\n  |Delta kappa|/kappa = " << max_kappa_err
              << std::endl;
}

// Tum TEST senaryolarinda S-GEO-04 dogrulamasi
TEST_F(SGEO04Test, AllScenariosOverlap)
{
    struct Scenario { double R_eq; double r0; };
    const Scenario scenarios[] = {
        {73.0, 22.0}, {152.4, 45.0}, {150.0, 10.0}, {200.0, 50.0}
    };

    for (const auto& sc : scenarios) {
        EllipsoidalProfile  ellip(1.0);
        HemisphericalProfile hemi;

        auto te = ellip.generateProfile(sc.R_eq, sc.r0, 500);
        auto th = hemi.generateProfile(sc.R_eq, sc.r0, 500);

        EXPECT_NEAR(te.metadata().s_total, th.metadata().s_total, 1e-10)
            << "R=" << sc.R_eq << " r0=" << sc.r0;
        EXPECT_NEAR(te.metadata().h_dome, th.metadata().h_dome, 1e-10)
            << "R=" << sc.R_eq << " r0=" << sc.r0;
    }
}

// ===========================================================================
// S-GEO-03: k_min sinir dogrulamasi
// ===========================================================================
TEST(EllipsoidalProfileSGEO03, KBelowMinThrows)
{
    EXPECT_THROW(EllipsoidalProfile(0.10), std::invalid_argument);
    EXPECT_THROW(EllipsoidalProfile(0.14), std::invalid_argument);
    EXPECT_THROW(EllipsoidalProfile(0.0),  std::invalid_argument);
    EXPECT_THROW(EllipsoidalProfile(-0.5), std::invalid_argument);
}

TEST(EllipsoidalProfileSGEO03, KAtMinSucceeds)
{
    EXPECT_NO_THROW(EllipsoidalProfile(0.15));
    EllipsoidalProfile prof(0.15);
    auto table = prof.generateProfile(152.4, 45.0, 500);
    EXPECT_TRUE(table.isValid());
    EXPECT_DOUBLE_EQ(prof.k(), 0.15);

    // Meta-veri: kappa_eq = 1/(R_eq * k^2) = 1/(152.4 * 0.0225) = 291.58...
    // Asiri yuksek egrilik — dome cok yassi
    EXPECT_NEAR(table.metadata().kappa_eq,
                1.0 / (152.4 * 0.15 * 0.15), 1e-8);
}

TEST(EllipsoidalProfileSGEO03, KAboveMinSucceeds)
{
    EXPECT_NO_THROW(EllipsoidalProfile(0.30));
    EXPECT_NO_THROW(EllipsoidalProfile(1.00));
    EXPECT_NO_THROW(EllipsoidalProfile(2.00));
}

// ===========================================================================
// Girdi dogrulama testleri (Karar-5 + Karar-17)
// ===========================================================================
TEST(EllipsoidalProfileValidation, NegativeReqThrows)
{
    EllipsoidalProfile prof(0.7);
    EXPECT_THROW(prof.generateProfile(-100.0, 30.0, 500), std::invalid_argument);
}

TEST(EllipsoidalProfileValidation, ZeroR0Throws)
{
    EllipsoidalProfile prof(0.7);
    EXPECT_THROW(prof.generateProfile(100.0, 0.0, 500), std::invalid_argument);
}

TEST(EllipsoidalProfileValidation, R0EqualReqThrows)
{
    EllipsoidalProfile prof(0.7);
    EXPECT_THROW(prof.generateProfile(100.0, 100.0, 500), std::invalid_argument);
}

TEST(EllipsoidalProfileValidation, R0GreaterThanReqThrows)
{
    EllipsoidalProfile prof(0.7);
    EXPECT_THROW(prof.generateProfile(100.0, 120.0, 500), std::invalid_argument);
}

TEST(EllipsoidalProfileValidation, TooFewPointsThrows)
{
    EllipsoidalProfile prof(0.7);
    EXPECT_THROW(prof.generateProfile(100.0, 30.0, 1), std::invalid_argument);
    EXPECT_THROW(prof.generateProfile(100.0, 30.0, 0), std::invalid_argument);
}

// ===========================================================================
// DomeType ve factory
// ===========================================================================
TEST(EllipsoidalProfileIdentity, DomeTypeAndName)
{
    EllipsoidalProfile prof(0.7);
    EXPECT_EQ(prof.domeType(), DomeType::Ellipsoidal);
    EXPECT_STREQ(prof.name(), "Ellipsoidal");
    EXPECT_DOUBLE_EQ(prof.k(), 0.7);
}

TEST(EllipsoidalProfileIdentity, FactoryCreation)
{
    auto prof = createProfile(DomeType::Ellipsoidal, 0.6);
    ASSERT_NE(prof, nullptr);
    EXPECT_EQ(prof->domeType(), DomeType::Ellipsoidal);

    auto table = prof->generateProfile(100.0, 30.0, 500);
    EXPECT_TRUE(table.isValid());
}

// ===========================================================================
// Uc durumlar
// ===========================================================================
TEST(EllipsoidalProfileEdgeCases, TallDome)
{
    // k=2.0 — uzun dome
    EllipsoidalProfile prof(2.0);
    auto table = prof.generateProfile(152.4, 45.0, 500);
    EXPECT_TRUE(table.isValid());

    const auto& meta = table.metadata();
    // h_dome = k * R_eq * sin(theta_p) > R_eq (uzun dome)
    EXPECT_GT(meta.h_dome, meta.R_eq);
    EXPECT_DOUBLE_EQ(meta.aspect_r, 2.0);

    std::cout << "[Tall dome k=2.0] h_dome=" << meta.h_dome
              << " mm, s_total=" << meta.s_total << " mm" << std::endl;
}

TEST(EllipsoidalProfileEdgeCases, FlatDome)
{
    // k=0.15 — en yassi izin verilen dome (S-GEO-03 sinirda)
    EllipsoidalProfile prof(0.15);
    auto table = prof.generateProfile(152.4, 45.0, 500);
    EXPECT_TRUE(table.isValid());

    const auto& meta = table.metadata();
    EXPECT_LT(meta.h_dome, meta.R_eq * 0.2);  // cok sik

    std::cout << "[Flat dome k=0.15] h_dome=" << meta.h_dome
              << " mm, kappa_eq=" << meta.kappa_eq << " 1/mm" << std::endl;
}
