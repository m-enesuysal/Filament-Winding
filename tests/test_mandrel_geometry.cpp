// =============================================================================
// test_mandrel_geometry.cpp — MandrelGeometry Birim Testleri
// =============================================================================
// Phase-1b S6: Silindir + dome birlesimi dogrulama
//
// Tum dome tipleriyle entegrasyon testi:
//   - Hemispherical (R=73, r0=22, L_cyl=200)
//   - Ellipsoidal   (R=152.4, r0=45, L_cyl=300, k=0.7)
//   - Isotensoid    (R=150, r0=10, L_cyl=250)
//
// Karar-11 Katman 3: Junction band tolerance +-1e-6 mm
// Karar-21: Immutable design — tum hesap constructor'da
// =============================================================================

#include <gtest/gtest.h>
#include <cmath>
#include <iostream>

#include "geometry/mandrel_geometry.h"
#include "geometry/filament_types.h"

using namespace filament::geometry;

// ===========================================================================
// Parametrize edilmis test verisi — 3 dome tipi
// ===========================================================================
struct MandrelTestParams {
    const char* name;
    DomeType dome_type;
    double R_eq;
    double r0;
    double L_cyl;
    double k;
    std::size_t N_points;
};

static const MandrelTestParams MANDREL_PARAMS[] = {
    {"Hemispherical", DomeType::Hemispherical, 73.0,  22.0,  200.0, 1.0, 500},
    {"Ellipsoidal",   DomeType::Ellipsoidal,  152.4,  45.0,  300.0, 0.7, 500},
    {"Isotensoid",    DomeType::Isotensoid,   150.0,  10.0,  250.0, 1.0, 500},
};

// ===========================================================================
// Parametrize test fixture
// ===========================================================================
class MandrelGeometryTest
    : public ::testing::TestWithParam<MandrelTestParams> {
};

// ---------------------------------------------------------------------------
// T1: Construction — tum dome tipleriyle basarili olusturma
// ---------------------------------------------------------------------------
TEST_P(MandrelGeometryTest, ConstructionSucceeds)
{
    const auto& p = GetParam();
    EXPECT_NO_THROW(
        MandrelGeometry(p.dome_type, p.R_eq, p.r0, p.L_cyl, p.k, p.N_points)
    ) << p.name;
}

// ---------------------------------------------------------------------------
// T2: Metadata — R_eq, r0, L_cyl, domeType, totalLength
// ---------------------------------------------------------------------------
TEST_P(MandrelGeometryTest, MetadataAccessors)
{
    const auto& p = GetParam();
    MandrelGeometry mg(p.dome_type, p.R_eq, p.r0, p.L_cyl, p.k, p.N_points);

    EXPECT_DOUBLE_EQ(mg.R_eq(), p.R_eq) << p.name << " R_eq";
    EXPECT_DOUBLE_EQ(mg.r0(), p.r0) << p.name << " r0";
    EXPECT_DOUBLE_EQ(mg.L_cyl(), p.L_cyl) << p.name << " L_cyl";
    EXPECT_EQ(mg.domeType(), p.dome_type) << p.name << " domeType";

    // totalLength = 2 * s_dome + L_cyl
    const double s_dome = mg.domeMetadata().s_total;
    EXPECT_NEAR(mg.totalLength(), 2.0 * s_dome + p.L_cyl, 1e-12)
        << p.name << " totalLength";

    std::cout << "[" << p.name << "] s_dome=" << s_dome
              << " mm, s_total=" << mg.totalLength()
              << " mm, h_dome=" << mg.domeMetadata().h_dome << " mm"
              << std::endl;
}

// ---------------------------------------------------------------------------
// T3: Region detection — isOnDome1, isOnCylinder, isOnDome2
// ---------------------------------------------------------------------------
TEST_P(MandrelGeometryTest, RegionDetectionInterior)
{
    const auto& p = GetParam();
    MandrelGeometry mg(p.dome_type, p.R_eq, p.r0, p.L_cyl, p.k, p.N_points);

    const double s_dome = mg.domeMetadata().s_total;
    const double s_total = mg.totalLength();

    // Dome-1 ic noktasi
    const double s_d1 = s_dome * 0.5;
    EXPECT_TRUE(mg.isOnDome1(s_d1)) << p.name << " dome-1 ic";
    EXPECT_FALSE(mg.isOnCylinder(s_d1)) << p.name;
    EXPECT_FALSE(mg.isOnDome2(s_d1)) << p.name;

    // Silindir ic noktasi
    const double s_cyl = s_dome + p.L_cyl * 0.5;
    EXPECT_FALSE(mg.isOnDome1(s_cyl)) << p.name << " cylinder ic";
    EXPECT_TRUE(mg.isOnCylinder(s_cyl)) << p.name;
    EXPECT_FALSE(mg.isOnDome2(s_cyl)) << p.name;

    // Dome-2 ic noktasi
    const double s_d2 = s_dome + p.L_cyl + s_dome * 0.5;
    EXPECT_FALSE(mg.isOnDome1(s_d2)) << p.name << " dome-2 ic";
    EXPECT_FALSE(mg.isOnCylinder(s_d2)) << p.name;
    EXPECT_TRUE(mg.isOnDome2(s_d2)) << p.name;

    // Uc noktalar
    EXPECT_TRUE(mg.isOnDome1(0.0)) << p.name << " s=0";
    EXPECT_TRUE(mg.isOnDome2(s_total)) << p.name << " s=s_total";
}

// ---------------------------------------------------------------------------
// T3b: Region detection — junction band (Karar-11 Katman 3)
// ---------------------------------------------------------------------------
TEST_P(MandrelGeometryTest, RegionDetectionJunctionBand)
{
    const auto& p = GetParam();
    MandrelGeometry mg(p.dome_type, p.R_eq, p.r0, p.L_cyl, p.k, p.N_points);

    const double s_dome = mg.domeMetadata().s_total;
    const double J = tolerances::JUNCTION_BAND_TOL;

    // Dome-1/Silindir siniri (s_dome): band icinde dome tarafi
    EXPECT_TRUE(mg.isOnDome1(s_dome)) << p.name << " s_dome tam sinir → dome-1";
    EXPECT_TRUE(mg.isOnDome1(s_dome + J * 0.5))
        << p.name << " s_dome + J/2 → dome-1 (band icinde)";

    // Silindir/Dome-2 siniri (s_dome + L_cyl): band icinde dome tarafi
    const double s_junc2 = s_dome + p.L_cyl;
    EXPECT_TRUE(mg.isOnDome2(s_junc2))
        << p.name << " s_dome+L_cyl tam sinir → dome-2";
    EXPECT_TRUE(mg.isOnDome2(s_junc2 - J * 0.5))
        << p.name << " s_dome+L_cyl - J/2 → dome-2 (band icinde)";
}

// ---------------------------------------------------------------------------
// T4: C1 sureklilik — dome-silindir kavsaklarinda
// ---------------------------------------------------------------------------
TEST_P(MandrelGeometryTest, C1ContinuityAtJunctions)
{
    const auto& p = GetParam();
    MandrelGeometry mg(p.dome_type, p.R_eq, p.r0, p.L_cyl, p.k, p.N_points);

    const double s_dome = mg.domeMetadata().s_total;
    constexpr double TOL = 1e-6;

    // --- Kavsak 1: Dome-1 ekvator / Silindir sol (s = s_dome) ---
    {
        // Dome-1 tarafinda (s_dome - eps)
        auto dome_pt = mg.point(s_dome - 1e-10);
        // Silindir tarafinda (s_dome + 2*J)
        auto cyl_pt = mg.point(s_dome + 2.0 * tolerances::JUNCTION_BAND_TOL);

        ASSERT_TRUE(dome_pt.has_value()) << p.name << " dome-1 ekvator";
        ASSERT_TRUE(cyl_pt.has_value()) << p.name << " silindir sol";

        EXPECT_NEAR(dome_pt->rho, p.R_eq, TOL)
            << p.name << " kavsak-1 rho";
        EXPECT_NEAR(dome_pt->drho_ds, 0.0, TOL)
            << p.name << " kavsak-1 drho/ds";
        EXPECT_NEAR(dome_pt->dx_ds, 1.0, TOL)
            << p.name << " kavsak-1 dx/ds";

        // Silindir noktasiyla uyum
        EXPECT_NEAR(dome_pt->rho, cyl_pt->rho, TOL)
            << p.name << " kavsak-1 rho uyumu";
        EXPECT_NEAR(dome_pt->drho_ds, cyl_pt->drho_ds, TOL)
            << p.name << " kavsak-1 drho/ds uyumu";
        EXPECT_NEAR(dome_pt->dx_ds, cyl_pt->dx_ds, TOL)
            << p.name << " kavsak-1 dx/ds uyumu";
    }

    // --- Kavsak 2: Silindir sag / Dome-2 ekvator (s = s_dome + L_cyl) ---
    {
        const double s_junc2 = s_dome + p.L_cyl;

        // Silindir tarafinda (s_junc2 - 2*J)
        auto cyl_pt = mg.point(s_junc2 - 2.0 * tolerances::JUNCTION_BAND_TOL);
        // Dome-2 tarafinda (s_junc2 + eps)
        auto dome_pt = mg.point(s_junc2 + 1e-10);

        ASSERT_TRUE(cyl_pt.has_value()) << p.name << " silindir sag";
        ASSERT_TRUE(dome_pt.has_value()) << p.name << " dome-2 ekvator";

        EXPECT_NEAR(dome_pt->rho, p.R_eq, TOL)
            << p.name << " kavsak-2 rho";
        EXPECT_NEAR(dome_pt->drho_ds, 0.0, TOL)
            << p.name << " kavsak-2 drho/ds";
        EXPECT_NEAR(dome_pt->dx_ds, 1.0, TOL)
            << p.name << " kavsak-2 dx/ds";

        EXPECT_NEAR(dome_pt->rho, cyl_pt->rho, TOL)
            << p.name << " kavsak-2 rho uyumu";
        EXPECT_NEAR(dome_pt->drho_ds, cyl_pt->drho_ds, TOL)
            << p.name << " kavsak-2 drho/ds uyumu";
        EXPECT_NEAR(dome_pt->dx_ds, cyl_pt->dx_ds, TOL)
            << p.name << " kavsak-2 dx/ds uyumu";
    }
}

// ---------------------------------------------------------------------------
// T5: Dome-1 sinir kosullari — s=0 (polar aciklik)
// ---------------------------------------------------------------------------
TEST_P(MandrelGeometryTest, Dome1PolarBoundary)
{
    const auto& p = GetParam();
    MandrelGeometry mg(p.dome_type, p.R_eq, p.r0, p.L_cyl, p.k, p.N_points);

    auto pt = mg.point(0.0);
    ASSERT_TRUE(pt.has_value()) << p.name;

    EXPECT_NEAR(pt->rho, p.r0, 1e-4)
        << p.name << " rho(0) = r0";
    EXPECT_NEAR(pt->x_local, 0.0, 1e-4)
        << p.name << " x(0) = 0";

    // Polar noktasinda drho/ds > 0 (dome-1 ters: ekvator yonune dogru)
    // ve birim teget vektor kisitlamasi
    EXPECT_GT(pt->drho_ds, -1e-6)
        << p.name << " drho/ds(0) >= 0 (ekvator yonune)";
    const double norm = pt->drho_ds * pt->drho_ds + pt->dx_ds * pt->dx_ds;
    EXPECT_NEAR(norm, 1.0, 1e-6)
        << p.name << " birim teget polar";

    std::cout << "[" << p.name << "] Dome-1 polar: rho=" << pt->rho
              << ", x=" << pt->x_local
              << ", drho/ds=" << pt->drho_ds
              << ", dx/ds=" << pt->dx_ds << std::endl;
}

// ---------------------------------------------------------------------------
// T6: Dome-2 sinir kosullari — s=s_total (polar aciklik)
// ---------------------------------------------------------------------------
TEST_P(MandrelGeometryTest, Dome2PolarBoundary)
{
    const auto& p = GetParam();
    MandrelGeometry mg(p.dome_type, p.R_eq, p.r0, p.L_cyl, p.k, p.N_points);

    const double s_total = mg.totalLength();
    const double h_dome = mg.domeMetadata().h_dome;

    auto pt = mg.point(s_total);
    ASSERT_TRUE(pt.has_value()) << p.name;

    EXPECT_NEAR(pt->rho, p.r0, 1e-4)
        << p.name << " rho(s_total) = r0";
    EXPECT_NEAR(pt->x_local, 2.0 * h_dome + p.L_cyl, 1e-4)
        << p.name << " x(s_total) = 2*h_dome + L_cyl";

    // Polar noktasinda drho/ds <= 0 (dome-2 ayni yon: ekvatordan uzaklasir)
    // ve birim teget vektor kisitlamasi
    EXPECT_LT(pt->drho_ds, 1e-6)
        << p.name << " drho/ds(s_total) <= 0 (polar yonune)";
    const double norm = pt->drho_ds * pt->drho_ds + pt->dx_ds * pt->dx_ds;
    EXPECT_NEAR(norm, 1.0, 1e-6)
        << p.name << " birim teget polar dome-2";

    std::cout << "[" << p.name << "] Dome-2 polar: rho=" << pt->rho
              << ", x=" << pt->x_local
              << " (expected " << 2.0 * h_dome + p.L_cyl << ")"
              << ", drho/ds=" << pt->drho_ds
              << ", dx/ds=" << pt->dx_ds << std::endl;
}

// ---------------------------------------------------------------------------
// T7: Silindir ozellikleri — rho=R_eq, kappa_m=0, dx/ds=1
// ---------------------------------------------------------------------------
TEST_P(MandrelGeometryTest, CylinderProperties)
{
    const auto& p = GetParam();
    MandrelGeometry mg(p.dome_type, p.R_eq, p.r0, p.L_cyl, p.k, p.N_points);

    const double s_dome = mg.domeMetadata().s_total;
    const double h_dome = mg.domeMetadata().h_dome;

    // Silindir orta noktasi
    const double s_mid = s_dome + p.L_cyl * 0.5;
    auto pt = mg.point(s_mid);
    ASSERT_TRUE(pt.has_value()) << p.name;

    EXPECT_DOUBLE_EQ(pt->rho, p.R_eq) << p.name << " rho = R_eq";
    EXPECT_DOUBLE_EQ(pt->drho_ds, 0.0) << p.name << " drho/ds = 0";
    EXPECT_DOUBLE_EQ(pt->dx_ds, 1.0) << p.name << " dx/ds = 1";
    EXPECT_DOUBLE_EQ(pt->kappa_m, 0.0) << p.name << " kappa_m = 0";

    // x_local = h_dome + L_cyl/2
    EXPECT_NEAR(pt->x_local, h_dome + p.L_cyl * 0.5, 1e-12)
        << p.name << " x_local silindir orta";
}

// ---------------------------------------------------------------------------
// T8: Simetri — rho(s) = rho(s_total - s)
// ---------------------------------------------------------------------------
TEST_P(MandrelGeometryTest, RadialSymmetry)
{
    const auto& p = GetParam();
    MandrelGeometry mg(p.dome_type, p.R_eq, p.r0, p.L_cyl, p.k, p.N_points);

    const double s_total = mg.totalLength();

    // Dome bolgesinde 10 nokta test et
    const double s_dome = mg.domeMetadata().s_total;
    for (int i = 0; i <= 10; ++i) {
        const double s_left = s_dome * i / 10.0;
        const double s_right = s_total - s_left;

        auto pt_left = mg.point(s_left);
        auto pt_right = mg.point(s_right);

        ASSERT_TRUE(pt_left.has_value()) << p.name << " s=" << s_left;
        ASSERT_TRUE(pt_right.has_value()) << p.name << " s=" << s_right;

        EXPECT_NEAR(pt_left->rho, pt_right->rho, 1e-10)
            << p.name << " simetri i=" << i
            << " rho_left=" << pt_left->rho
            << " rho_right=" << pt_right->rho;
    }
}

// ---------------------------------------------------------------------------
// T9: Monotoniklik — rho artar (dome-1), sabit (cyl), azalir (dome-2)
// ---------------------------------------------------------------------------
TEST_P(MandrelGeometryTest, RhoMonotonicity)
{
    const auto& p = GetParam();
    MandrelGeometry mg(p.dome_type, p.R_eq, p.r0, p.L_cyl, p.k, p.N_points);

    const double s_dome = mg.domeMetadata().s_total;
    const double s_total = mg.totalLength();
    constexpr int N_SAMPLES = 50;

    // Dome-1: rho artmali (polar → ekvator)
    for (int i = 1; i <= N_SAMPLES; ++i) {
        const double s_prev = s_dome * (i - 1) / N_SAMPLES;
        const double s_curr = s_dome * i / N_SAMPLES;
        auto p0 = mg.point(s_prev);
        auto p1 = mg.point(s_curr);
        ASSERT_TRUE(p0.has_value() && p1.has_value());
        EXPECT_LE(p0->rho, p1->rho + 1e-6)
            << p.name << " dome-1 rho artmali i=" << i;
    }

    // Silindir: rho sabit
    for (int i = 1; i <= N_SAMPLES; ++i) {
        const double s_curr = s_dome + p.L_cyl * i / N_SAMPLES;
        auto pt = mg.point(s_curr);
        ASSERT_TRUE(pt.has_value());
        EXPECT_NEAR(pt->rho, p.R_eq, 1e-12)
            << p.name << " silindir rho sabit i=" << i;
    }

    // Dome-2: rho azalmali (ekvator → polar)
    // Not: s_dome + L_cyl + s_dome hesabi floating-point'ta s_total'den
    // farkli olabilir, bu yuzden std::min ile sinirla
    for (int i = 1; i <= N_SAMPLES; ++i) {
        const double s_prev = std::min(
            s_dome + p.L_cyl + s_dome * (i - 1.0) / N_SAMPLES, s_total);
        const double s_curr = std::min(
            s_dome + p.L_cyl + s_dome * static_cast<double>(i) / N_SAMPLES,
            s_total);
        auto p0 = mg.point(s_prev);
        auto p1 = mg.point(s_curr);
        ASSERT_TRUE(p0.has_value() && p1.has_value())
            << p.name << " dome-2 nullopt i=" << i
            << " s_prev=" << s_prev << " s_curr=" << s_curr;
        EXPECT_GE(p0->rho, p1->rho - 1e-6)
            << p.name << " dome-2 rho azalmali i=" << i;
    }
}

// ---------------------------------------------------------------------------
// T10: Aralik disi sorgular — nullopt donmeli
// ---------------------------------------------------------------------------
TEST_P(MandrelGeometryTest, OutOfRangeReturnsNullopt)
{
    const auto& p = GetParam();
    MandrelGeometry mg(p.dome_type, p.R_eq, p.r0, p.L_cyl, p.k, p.N_points);

    EXPECT_FALSE(mg.point(-1.0).has_value()) << p.name << " s < 0";
    EXPECT_FALSE(mg.point(-1e-10).has_value()) << p.name << " s = -eps";
    EXPECT_FALSE(mg.point(mg.totalLength() + 1.0).has_value())
        << p.name << " s > s_total";
}

// ---------------------------------------------------------------------------
// T11: Clairaut sabiti
// ---------------------------------------------------------------------------
TEST_P(MandrelGeometryTest, ClairautConstant)
{
    const auto& p = GetParam();
    MandrelGeometry mg(p.dome_type, p.R_eq, p.r0, p.L_cyl, p.k, p.N_points);

    const double alpha = 0.5;  // rad
    const double rho_test = 80.0;
    const double clairaut = mg.windingAngleToClairaut(alpha, rho_test);
    EXPECT_NEAR(clairaut, rho_test * std::sin(alpha), 1e-12)
        << p.name << " Clairaut = rho * sin(alpha)";
}

// ---------------------------------------------------------------------------
// T12: x_local monotoniklik — tum mandrel boyunca x artmali
// ---------------------------------------------------------------------------
TEST_P(MandrelGeometryTest, XLocalMonotonicallyIncreasing)
{
    const auto& p = GetParam();
    MandrelGeometry mg(p.dome_type, p.R_eq, p.r0, p.L_cyl, p.k, p.N_points);

    const double s_total = mg.totalLength();
    constexpr int N_SAMPLES = 200;

    double prev_x = -1.0;
    for (int i = 0; i <= N_SAMPLES; ++i) {
        const double s = s_total * i / N_SAMPLES;
        auto pt = mg.point(s);
        ASSERT_TRUE(pt.has_value()) << p.name << " i=" << i;
        EXPECT_GT(pt->x_local, prev_x - 1e-12)
            << p.name << " x_local monoton artan i=" << i
            << " prev=" << prev_x << " curr=" << pt->x_local;
        prev_x = pt->x_local;
    }
}

// ---------------------------------------------------------------------------
// T13: Birim teget vektor — tum mandrel boyunca |drho/ds|^2 + |dx/ds|^2 = 1
// ---------------------------------------------------------------------------
TEST_P(MandrelGeometryTest, UnitTangentVector)
{
    const auto& p = GetParam();
    MandrelGeometry mg(p.dome_type, p.R_eq, p.r0, p.L_cyl, p.k, p.N_points);

    const double s_total = mg.totalLength();
    constexpr int N_SAMPLES = 200;

    double max_err = 0.0;
    for (int i = 0; i <= N_SAMPLES; ++i) {
        const double s = s_total * i / N_SAMPLES;
        auto pt = mg.point(s);
        ASSERT_TRUE(pt.has_value());

        const double norm_sq = pt->drho_ds * pt->drho_ds + pt->dx_ds * pt->dx_ds;
        const double err = std::abs(norm_sq - 1.0);
        max_err = std::max(max_err, err);

        EXPECT_LT(err, 1e-6)
            << p.name << " birim teget i=" << i << " s=" << s;
    }

    std::cout << "[" << p.name << "] Birim teget maks hata: "
              << max_err << std::endl;
}

// ---------------------------------------------------------------------------
// T14: domeMetadata() tutarliligi
// ---------------------------------------------------------------------------
TEST_P(MandrelGeometryTest, DomeMetadataConsistency)
{
    const auto& p = GetParam();
    MandrelGeometry mg(p.dome_type, p.R_eq, p.r0, p.L_cyl, p.k, p.N_points);

    const auto& meta = mg.domeMetadata();
    EXPECT_DOUBLE_EQ(meta.R_eq, p.R_eq) << p.name;
    EXPECT_DOUBLE_EQ(meta.r0, p.r0) << p.name;
    EXPECT_GT(meta.s_total, 0.0) << p.name << " s_total > 0";
    EXPECT_GT(meta.h_dome, 0.0) << p.name << " h_dome > 0";
    EXPECT_GT(meta.A_dome, 0.0) << p.name << " A_dome > 0";
}

// ---------------------------------------------------------------------------
// Parametrize instantiation
// ---------------------------------------------------------------------------
INSTANTIATE_TEST_SUITE_P(
    AllDomeTypes,
    MandrelGeometryTest,
    ::testing::ValuesIn(MANDREL_PARAMS),
    [](const ::testing::TestParamInfo<MandrelTestParams>& info) {
        return std::string(info.param.name);
    }
);

// ===========================================================================
// Girdi dogrulama testleri (Karar-5 + Karar-17)
// ===========================================================================
TEST(MandrelGeometryValidation, NegativeReqThrows)
{
    EXPECT_THROW(
        MandrelGeometry(DomeType::Hemispherical, -100.0, 30.0, 200.0),
        std::invalid_argument);
}

TEST(MandrelGeometryValidation, ZeroR0Throws)
{
    EXPECT_THROW(
        MandrelGeometry(DomeType::Hemispherical, 100.0, 0.0, 200.0),
        std::invalid_argument);
}

TEST(MandrelGeometryValidation, R0GreaterEqualReqThrows)
{
    EXPECT_THROW(
        MandrelGeometry(DomeType::Hemispherical, 100.0, 100.0, 200.0),
        std::invalid_argument);
    EXPECT_THROW(
        MandrelGeometry(DomeType::Hemispherical, 100.0, 120.0, 200.0),
        std::invalid_argument);
}

TEST(MandrelGeometryValidation, NegativeLcylThrows)
{
    EXPECT_THROW(
        MandrelGeometry(DomeType::Hemispherical, 100.0, 30.0, -10.0),
        std::invalid_argument);
}

TEST(MandrelGeometryValidation, TooFewPointsThrows)
{
    EXPECT_THROW(
        MandrelGeometry(DomeType::Hemispherical, 100.0, 30.0, 200.0, 1.0, 1),
        std::invalid_argument);
    EXPECT_THROW(
        MandrelGeometry(DomeType::Hemispherical, 100.0, 30.0, 200.0, 1.0, 0),
        std::invalid_argument);
}

// ===========================================================================
// Uc durum: L_cyl = 0 (dome'dan dome'a direkt gecis)
// ===========================================================================
TEST(MandrelGeometryEdgeCases, ZeroCylinderLength)
{
    MandrelGeometry mg(DomeType::Hemispherical, 100.0, 30.0, 0.0);

    const double s_dome = mg.domeMetadata().s_total;

    // totalLength = 2 * s_dome (silindir yok)
    EXPECT_NEAR(mg.totalLength(), 2.0 * s_dome, 1e-12);

    // Kavsak noktasinda (s = s_dome): dome-1 ekvator = dome-2 ekvator
    auto pt = mg.point(s_dome);
    ASSERT_TRUE(pt.has_value());
    EXPECT_NEAR(pt->rho, 100.0, 1e-6) << "L_cyl=0 kavsak rho";

    // Simetri korunuyor mu?
    auto pt_left = mg.point(s_dome * 0.3);
    auto pt_right = mg.point(2.0 * s_dome - s_dome * 0.3);
    ASSERT_TRUE(pt_left.has_value() && pt_right.has_value());
    EXPECT_NEAR(pt_left->rho, pt_right->rho, 1e-10)
        << "L_cyl=0 simetri";

    // isOnCylinder false olmali (silindir yok)
    EXPECT_FALSE(mg.isOnCylinder(s_dome))
        << "L_cyl=0: silindir bolge yok";
}

// ===========================================================================
// Uc durum: Ellipsoidal k=1 (hemispherical ile ozmanlasma)
// ===========================================================================
TEST(MandrelGeometryEdgeCases, EllipsoidalK1MatchesHemispherical)
{
    const double R = 100.0, r0 = 30.0, L = 200.0;

    MandrelGeometry mg_hemi(DomeType::Hemispherical, R, r0, L);
    MandrelGeometry mg_elli(DomeType::Ellipsoidal, R, r0, L, 1.0);

    // Total uzunluk uyumlu olmali
    EXPECT_NEAR(mg_hemi.totalLength(), mg_elli.totalLength(), 0.1)
        << "k=1 total length";

    // Dome ortasinda rho uyumlu
    const double s_mid = mg_hemi.domeMetadata().s_total * 0.5;
    auto pt_h = mg_hemi.point(s_mid);
    auto pt_e = mg_elli.point(s_mid);
    ASSERT_TRUE(pt_h.has_value() && pt_e.has_value());
    EXPECT_NEAR(pt_h->rho, pt_e->rho, 0.1)
        << "k=1 rho overlap";
}

// ===========================================================================
// Tum dome tipleriyle mandrel olusturma (entegrasyon testi)
// ===========================================================================
TEST(MandrelGeometryIntegration, AllDomeTypesQueryable)
{
    const double R = 100.0, r0 = 30.0, L = 200.0;

    // Her dome tipi icin basarili olusturma + 3 bolge sorgusu
    for (auto dt : {DomeType::Hemispherical,
                    DomeType::Ellipsoidal,
                    DomeType::Isotensoid}) {
        double k = (dt == DomeType::Ellipsoidal) ? 0.8 : 1.0;
        MandrelGeometry mg(dt, R, r0, L, k);

        const double s_dome = mg.domeMetadata().s_total;

        // 3 bolge sorgulanabilir
        auto p1 = mg.point(s_dome * 0.5);           // Dome-1
        auto p2 = mg.point(s_dome + L * 0.5);        // Silindir
        auto p3 = mg.point(s_dome + L + s_dome * 0.5); // Dome-2

        ASSERT_TRUE(p1.has_value());
        ASSERT_TRUE(p2.has_value());
        ASSERT_TRUE(p3.has_value());

        // Temel fizik kontrolleri
        EXPECT_GT(p1->rho, r0);
        EXPECT_LT(p1->rho, R);
        EXPECT_DOUBLE_EQ(p2->rho, R);
        EXPECT_GT(p3->rho, r0);
        EXPECT_LT(p3->rho, R);
    }
}
