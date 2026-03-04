// =============================================================================
// test_meridian_lookup_table.cpp — MeridianLookupTable Birim Testleri
// =============================================================================
// Phase-1b S2: Kubik spline interpolasyon ve binary search dogrulama
//
// Referans profil: Hemispherical dome (R=100, r0=30)
//   - Analitik formuller (kapali form) ile karsilastirma
//   - Karar-11 Katman 2 toleranslari:
//       Pozisyon  |eps| < 1e-4 mm  (mutlak)
//       Turev     |eps| < 1e-6     (mutlak, boyutsuz)
//       Egrilik   |eps_rel| < 1e-4 (bagil, %0.01)
//
// Test senaryolari:
//   1. Tablo build + gecerlilik
//   2. Dugum noktalarinda tam eslesme
//   3. Ara noktalarda Karar-11 toleranslari
//   4. Araligi disinda nullopt
//   5. Bos/gecersiz tablo
//   6. Meta-veri koruma
//   7. Lineer interpolasyon (2 nokta)
//   8. Gecersiz girdi dogrulama
// =============================================================================

#include <gtest/gtest.h>
#include <cmath>
#include <vector>
#include <optional>

#include "geometry/meridian_lookup_table.h"
#include "geometry/filament_types.h"

using namespace filament::geometry;

// ===========================================================================
// Test yardimlari: Hemispherical dome analitik profil uretici
// ===========================================================================
// Hemispherical dome icin:
//   theta: 0 (ekvator) -> theta_max (polar aciklik)
//   rho(theta) = R * cos(theta)
//   x(theta)   = R * sin(theta)
//   s(theta)    = R * theta         (yay uzunlugu)
//   drho/ds    = -sin(theta)
//   dx/ds      = cos(theta)
//   kappa_m    = 1/R               (sabit)
// ===========================================================================

struct HemisphericalTestData {
    double R_eq;
    double r0;
    double theta_max;
    double s_total;
    std::size_t N;

    std::vector<double> s;
    std::vector<double> rho;
    std::vector<double> x;
    std::vector<double> drho_ds;
    std::vector<double> dx_ds;
    std::vector<double> kappa_m;
    ProfileMetadata meta;
};

static HemisphericalTestData generateHemisphericalData(
    double R_eq, double r0, std::size_t N)
{
    HemisphericalTestData data;
    data.R_eq = R_eq;
    data.r0 = r0;
    data.theta_max = std::acos(r0 / R_eq);
    data.s_total = R_eq * data.theta_max;
    data.N = N;

    data.s.resize(N);
    data.rho.resize(N);
    data.x.resize(N);
    data.drho_ds.resize(N);
    data.dx_ds.resize(N);
    data.kappa_m.resize(N);

    for (std::size_t i = 0; i < N; ++i) {
        const double theta = data.theta_max
            * static_cast<double>(i) / static_cast<double>(N - 1);

        data.s[i]       = R_eq * theta;
        data.rho[i]     = R_eq * std::cos(theta);
        data.x[i]       = R_eq * std::sin(theta);
        data.drho_ds[i] = -std::sin(theta);
        data.dx_ds[i]   = std::cos(theta);
        data.kappa_m[i] = 1.0 / R_eq;
    }

    // Meta-veri
    data.meta.R_eq      = R_eq;
    data.meta.r0        = r0;
    data.meta.s_total   = data.s_total;
    data.meta.h_dome    = R_eq * std::sin(data.theta_max);
    data.meta.A_dome    = 2.0 * constants::PI * R_eq * R_eq
                        * (1.0 - std::cos(data.theta_max));
    data.meta.kappa_eq  = 1.0 / R_eq;
    data.meta.kappa_pol = 1.0 / R_eq;
    data.meta.alpha_w   = std::asin(r0 / R_eq);
    data.meta.aspect_r  = 1.0;

    return data;
}

// Analitik sorgu: verilen s'de hemispherical dome degerleri
static MeridianPoint hemisphericalAnalytic(double s, double R_eq)
{
    const double theta = s / R_eq;
    MeridianPoint pt;
    pt.s       = s;
    pt.rho     = R_eq * std::cos(theta);
    pt.x_local = R_eq * std::sin(theta);
    pt.drho_ds = -std::sin(theta);
    pt.dx_ds   = std::cos(theta);
    pt.kappa_m = 1.0 / R_eq;
    return pt;
}

// ===========================================================================
// Test fixture
// ===========================================================================
class MeridianLookupTableTest : public ::testing::Test {
protected:
    // Standart test parametreleri: R=100, r0=30, N=500
    static constexpr double R_eq = 100.0;
    static constexpr double r0   = 30.0;
    static constexpr std::size_t N_points = 500;

    // Karar-11 Katman 2 toleranslari
    static constexpr double POS_TOL  = tolerances::POSITION_ABS_TOL;   // 1e-4 mm
    static constexpr double DER_TOL  = tolerances::DERIVATIVE_ABS_TOL; // 1e-6
    static constexpr double CURV_TOL = tolerances::CURVATURE_REL_TOL;  // 1e-4

    HemisphericalTestData testData_;
    MeridianLookupTable table_;

    void SetUp() override {
        testData_ = generateHemisphericalData(R_eq, r0, N_points);
        table_.build(testData_.s, testData_.rho, testData_.x,
                     testData_.drho_ds, testData_.dx_ds, testData_.kappa_m,
                     testData_.meta);
    }
};

// ===========================================================================
// 1. Tablo build ve gecerlilik testleri
// ===========================================================================
TEST_F(MeridianLookupTableTest, BuildSucceeds)
{
    EXPECT_TRUE(table_.isValid());
    EXPECT_EQ(table_.size(), N_points);
}

TEST_F(MeridianLookupTableTest, RawDataPreserved)
{
    const auto& raw_s   = table_.rawS();
    const auto& raw_rho = table_.rawRho();
    const auto& raw_x   = table_.rawX();

    ASSERT_EQ(raw_s.size(), N_points);
    ASSERT_EQ(raw_rho.size(), N_points);
    ASSERT_EQ(raw_x.size(), N_points);

    for (std::size_t i = 0; i < N_points; ++i) {
        EXPECT_DOUBLE_EQ(raw_s[i], testData_.s[i]);
        EXPECT_DOUBLE_EQ(raw_rho[i], testData_.rho[i]);
        EXPECT_DOUBLE_EQ(raw_x[i], testData_.x[i]);
    }
}

// ===========================================================================
// 2. Dugum noktalarinda tam eslesme
// ===========================================================================
TEST_F(MeridianLookupTableTest, ExactAtKnots)
{
    // Spline dugum noktalarinda a_i = y_i, t = 0 oldugu icin tam eslesme
    for (std::size_t i = 0; i < N_points; ++i) {
        auto result = table_.query(testData_.s[i]);
        ASSERT_TRUE(result.has_value()) << "Dugum noktasi " << i << " nullopt dondurdu";

        // Kayan nokta hassasiyetinde eslesme (a_i = y_i, t ≈ 0)
        EXPECT_NEAR(result->rho,     testData_.rho[i],     1e-12)
            << "rho dugum " << i;
        EXPECT_NEAR(result->x_local, testData_.x[i],       1e-12)
            << "x dugum " << i;
        EXPECT_NEAR(result->drho_ds, testData_.drho_ds[i], 1e-12)
            << "drho_ds dugum " << i;
        EXPECT_NEAR(result->dx_ds,   testData_.dx_ds[i],   1e-12)
            << "dx_ds dugum " << i;
        EXPECT_NEAR(result->kappa_m, testData_.kappa_m[i], 1e-12)
            << "kappa_m dugum " << i;
    }
}

// ===========================================================================
// 3. Ara noktalarda Karar-11 Katman 2 toleranslari
// ===========================================================================
TEST_F(MeridianLookupTableTest, IntermediatePositionTolerance)
{
    // Her araligin orta noktasinda sorgula, analitik referansla karsilastir
    double max_rho_err = 0.0;
    double max_x_err   = 0.0;

    for (std::size_t i = 0; i + 1 < N_points; ++i) {
        const double s_mid = 0.5 * (testData_.s[i] + testData_.s[i + 1]);
        auto result = table_.query(s_mid);
        ASSERT_TRUE(result.has_value());

        const auto ref = hemisphericalAnalytic(s_mid, R_eq);

        const double rho_err = std::abs(result->rho - ref.rho);
        const double x_err   = std::abs(result->x_local - ref.x_local);

        max_rho_err = std::max(max_rho_err, rho_err);
        max_x_err   = std::max(max_x_err, x_err);

        EXPECT_LT(rho_err, POS_TOL)
            << "rho pozisyon hatasi aralik " << i << ": " << rho_err;
        EXPECT_LT(x_err, POS_TOL)
            << "x pozisyon hatasi aralik " << i << ": " << x_err;
    }

    // Bilgilendirme: maksimum hatalar
    std::cout << "[INFO] Maks rho pozisyon hatasi: " << max_rho_err
              << " mm (limit: " << POS_TOL << ")" << std::endl;
    std::cout << "[INFO] Maks x pozisyon hatasi:   " << max_x_err
              << " mm (limit: " << POS_TOL << ")" << std::endl;
}

TEST_F(MeridianLookupTableTest, IntermediateDerivativeTolerance)
{
    double max_drho_err = 0.0;
    double max_dx_err   = 0.0;

    for (std::size_t i = 0; i + 1 < N_points; ++i) {
        const double s_mid = 0.5 * (testData_.s[i] + testData_.s[i + 1]);
        auto result = table_.query(s_mid);
        ASSERT_TRUE(result.has_value());

        const auto ref = hemisphericalAnalytic(s_mid, R_eq);

        const double drho_err = std::abs(result->drho_ds - ref.drho_ds);
        const double dx_err   = std::abs(result->dx_ds   - ref.dx_ds);

        max_drho_err = std::max(max_drho_err, drho_err);
        max_dx_err   = std::max(max_dx_err, dx_err);

        EXPECT_LT(drho_err, DER_TOL)
            << "drho/ds hatasi aralik " << i << ": " << drho_err;
        EXPECT_LT(dx_err, DER_TOL)
            << "dx/ds hatasi aralik " << i << ": " << dx_err;
    }

    std::cout << "[INFO] Maks drho/ds hatasi: " << max_drho_err
              << " (limit: " << DER_TOL << ")" << std::endl;
    std::cout << "[INFO] Maks dx/ds hatasi:   " << max_dx_err
              << " (limit: " << DER_TOL << ")" << std::endl;
}

TEST_F(MeridianLookupTableTest, IntermediateCurvatureTolerance)
{
    double max_curv_rel_err = 0.0;
    const double kappa_ref  = 1.0 / R_eq;  // hemispherical: sabit

    for (std::size_t i = 0; i + 1 < N_points; ++i) {
        const double s_mid = 0.5 * (testData_.s[i] + testData_.s[i + 1]);
        auto result = table_.query(s_mid);
        ASSERT_TRUE(result.has_value());

        const double curv_rel_err =
            std::abs(result->kappa_m - kappa_ref) / std::abs(kappa_ref);

        max_curv_rel_err = std::max(max_curv_rel_err, curv_rel_err);

        EXPECT_LT(curv_rel_err, CURV_TOL)
            << "egrilik bagil hatasi aralik " << i << ": " << curv_rel_err;
    }

    std::cout << "[INFO] Maks egrilik bagil hatasi: " << max_curv_rel_err
              << " (limit: " << CURV_TOL << ")" << std::endl;
}

// ===========================================================================
// 3b. Rastgele s degerlerinde tolerans (aralik ortasi disinda)
// ===========================================================================
TEST_F(MeridianLookupTableTest, RandomQueryPositionTolerance)
{
    // Aralik icinde %25 ve %75 noktalarinda da test et
    for (std::size_t i = 0; i + 1 < N_points; i += 10) {
        const double h = testData_.s[i + 1] - testData_.s[i];

        for (double frac : {0.25, 0.75}) {
            const double s_q = testData_.s[i] + frac * h;
            auto result = table_.query(s_q);
            ASSERT_TRUE(result.has_value());

            const auto ref = hemisphericalAnalytic(s_q, R_eq);

            EXPECT_LT(std::abs(result->rho - ref.rho), POS_TOL);
            EXPECT_LT(std::abs(result->x_local - ref.x_local), POS_TOL);
            EXPECT_LT(std::abs(result->drho_ds - ref.drho_ds), DER_TOL);
            EXPECT_LT(std::abs(result->dx_ds - ref.dx_ds), DER_TOL);
        }
    }
}

// ===========================================================================
// 4. Aralik disinda nullopt
// ===========================================================================
TEST_F(MeridianLookupTableTest, OutOfRangeReturnsNullopt)
{
    // Asagida
    auto below = table_.query(-1.0);
    EXPECT_FALSE(below.has_value());

    // Yukarda
    auto above = table_.query(testData_.s_total + 1.0);
    EXPECT_FALSE(above.has_value());
}

TEST_F(MeridianLookupTableTest, BoundaryQueriesSucceed)
{
    // s = 0 (ekvator)
    auto at_start = table_.query(testData_.s.front());
    ASSERT_TRUE(at_start.has_value());
    EXPECT_NEAR(at_start->rho, R_eq, 1e-12);

    // s = s_total (polar aciklik)
    auto at_end = table_.query(testData_.s.back());
    ASSERT_TRUE(at_end.has_value());
    EXPECT_NEAR(at_end->rho, r0, 1e-10);
}

// ===========================================================================
// 5. Bos ve gecersiz tablo
// ===========================================================================
TEST_F(MeridianLookupTableTest, EmptyTableReturnsNullopt)
{
    MeridianLookupTable empty;
    EXPECT_FALSE(empty.isValid());
    EXPECT_EQ(empty.size(), 0u);

    auto result = empty.query(0.0);
    EXPECT_FALSE(result.has_value());
}

// ===========================================================================
// 6. Meta-veri koruma
// ===========================================================================
TEST_F(MeridianLookupTableTest, MetadataPreserved)
{
    const auto& meta = table_.metadata();

    EXPECT_DOUBLE_EQ(meta.R_eq, R_eq);
    EXPECT_DOUBLE_EQ(meta.r0, r0);
    EXPECT_NEAR(meta.s_total, testData_.s_total, 1e-10);
    EXPECT_NEAR(meta.kappa_eq, 1.0 / R_eq, 1e-14);
    EXPECT_NEAR(meta.kappa_pol, 1.0 / R_eq, 1e-14);
    EXPECT_NEAR(meta.aspect_r, 1.0, 1e-14);
}

// ===========================================================================
// 7. Lineer interpolasyon (2 nokta)
// ===========================================================================
TEST(MeridianLookupTableLinear, TwoPointLinear)
{
    MeridianLookupTable table;

    std::vector<double> s       = {0.0, 10.0};
    std::vector<double> rho     = {100.0, 90.0};
    std::vector<double> x       = {0.0, 5.0};
    std::vector<double> drho_ds = {-1.0, -1.0};
    std::vector<double> dx_ds   = {0.5, 0.5};
    std::vector<double> kappa   = {0.01, 0.01};

    ProfileMetadata meta;
    meta.R_eq = 100.0;
    meta.r0   = 30.0;

    table.build(s, rho, x, drho_ds, dx_ds, kappa, meta);
    ASSERT_TRUE(table.isValid());
    EXPECT_EQ(table.size(), 2u);

    // Orta noktada lineer interpolasyon
    auto mid = table.query(5.0);
    ASSERT_TRUE(mid.has_value());
    EXPECT_NEAR(mid->rho,     95.0, 1e-12);  // (100+90)/2
    EXPECT_NEAR(mid->x_local,  2.5, 1e-12);  // (0+5)/2
    EXPECT_NEAR(mid->drho_ds, -1.0, 1e-12);
    EXPECT_NEAR(mid->dx_ds,    0.5, 1e-12);
    EXPECT_NEAR(mid->kappa_m, 0.01, 1e-12);
}

// ===========================================================================
// 8. Gecersiz girdi dogrulama
// ===========================================================================
TEST(MeridianLookupTableValidation, DifferentSizes)
{
    MeridianLookupTable table;
    ProfileMetadata meta;

    std::vector<double> s3 = {0, 1, 2};
    std::vector<double> v2 = {0, 1};
    std::vector<double> v3 = {0, 1, 2};

    // rho boyutu farkli
    table.build(s3, v2, v3, v3, v3, v3, meta);
    EXPECT_FALSE(table.isValid());
}

TEST(MeridianLookupTableValidation, NonMonotonicS)
{
    MeridianLookupTable table;
    ProfileMetadata meta;

    // s monoton degil
    std::vector<double> s   = {0.0, 2.0, 1.0, 3.0};
    std::vector<double> val = {1.0, 2.0, 3.0, 4.0};

    table.build(s, val, val, val, val, val, meta);
    EXPECT_FALSE(table.isValid());
}

TEST(MeridianLookupTableValidation, DuplicateS)
{
    MeridianLookupTable table;
    ProfileMetadata meta;

    // s tekrarlayan deger iceriyor (kesinlikle artan degil)
    std::vector<double> s   = {0.0, 1.0, 1.0, 2.0};
    std::vector<double> val = {1.0, 2.0, 3.0, 4.0};

    table.build(s, val, val, val, val, val, meta);
    EXPECT_FALSE(table.isValid());
}

TEST(MeridianLookupTableValidation, SinglePoint)
{
    MeridianLookupTable table;
    ProfileMetadata meta;

    std::vector<double> s   = {0.0};
    std::vector<double> val = {1.0};

    table.build(s, val, val, val, val, val, meta);
    EXPECT_FALSE(table.isValid());
}

TEST(MeridianLookupTableValidation, EmptyVectors)
{
    MeridianLookupTable table;
    ProfileMetadata meta;

    std::vector<double> empty;
    table.build(empty, empty, empty, empty, empty, empty, meta);
    EXPECT_FALSE(table.isValid());
}

// ===========================================================================
// 9. Standart test senaryolari (Karar-16 TEST-01..04 parametreleri)
// ===========================================================================
struct TestScenario {
    const char* name;
    double R_eq;
    double r0;
};

class MeridianLookupTableScenarios
    : public ::testing::TestWithParam<TestScenario> {};

TEST_P(MeridianLookupTableScenarios, TolerancesMet)
{
    const auto& scenario = GetParam();
    const std::size_t N = 500;

    auto data = generateHemisphericalData(scenario.R_eq, scenario.r0, N);

    MeridianLookupTable table;
    table.build(data.s, data.rho, data.x,
                data.drho_ds, data.dx_ds, data.kappa_m, data.meta);
    ASSERT_TRUE(table.isValid());

    // Her araligin ortasinda Karar-11 Katman 2 toleranslarini kontrol et
    for (std::size_t i = 0; i + 1 < N; ++i) {
        const double s_mid = 0.5 * (data.s[i] + data.s[i + 1]);
        auto result = table.query(s_mid);
        ASSERT_TRUE(result.has_value());

        const auto ref = hemisphericalAnalytic(s_mid, scenario.R_eq);
        const double kappa_ref = 1.0 / scenario.R_eq;

        EXPECT_LT(std::abs(result->rho - ref.rho),
                  tolerances::POSITION_ABS_TOL)
            << scenario.name << " rho aralik " << i;
        EXPECT_LT(std::abs(result->x_local - ref.x_local),
                  tolerances::POSITION_ABS_TOL)
            << scenario.name << " x aralik " << i;
        EXPECT_LT(std::abs(result->drho_ds - ref.drho_ds),
                  tolerances::DERIVATIVE_ABS_TOL)
            << scenario.name << " drho/ds aralik " << i;
        EXPECT_LT(std::abs(result->dx_ds - ref.dx_ds),
                  tolerances::DERIVATIVE_ABS_TOL)
            << scenario.name << " dx/ds aralik " << i;
        EXPECT_LT(std::abs(result->kappa_m - kappa_ref) / std::abs(kappa_ref),
                  tolerances::CURVATURE_REL_TOL)
            << scenario.name << " kappa aralik " << i;
    }
}

INSTANTIATE_TEST_SUITE_P(
    StandardScenarios,
    MeridianLookupTableScenarios,
    ::testing::Values(
        TestScenario{"TEST-01_ASTM_Subscale",       73.0,  22.0},
        TestScenario{"TEST-02_Endustriyel_COPV",   152.4,  45.0},
        TestScenario{"TEST-03_Kucuk_Aciklik",      150.0,  10.0},
        TestScenario{"TEST-04_H2_Aerospace",       200.0,  50.0}
    )
);
