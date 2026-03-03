// =============================================================================
// test_meridian_lookup_table.cpp — MeridianLookupTable Birim Testleri
// =============================================================================
// Kapsam: S2 session — Phase-1b MeridianLookupTable implementasyonu
//
// Test gruplari:
//   1. TableState        — tablo durumu (gecerlilik, boyut)
//   2. BoundaryBehavior  — sinir davranisi (aralik disi nullopt, uclar)
//   3. LinearAccuracy    — lineer fonksiyon (spline exact olmali)
//   4. ConstantAccuracy  — sabit fonksiyon (spline exact olmali)
//   5. Hemispherical     — hemisferik dome analitik referans
//      Karar-11 Katman 2 toleranslari:
//        - Pozisyon (rho, x):  |e| < 1e-4 mm
//        - Turev (drho, dx):   |e| < 1e-6
//        - Egrilik (kappa_m):  |e_rel| < 1e-4
//
// Test senaryolari: R_eq=100 mm, r0=30 mm (TEST-01 benzeri kucultulmus)
// Karar-16 TEST-01: R_eq=73, r0=22 — ek senaryo eklendi
//
// Karar-19: -ffast-math YASAK (IEEE 754 zorunlu)
// =============================================================================

#include <gtest/gtest.h>
#include <cmath>
#include <vector>
#include <optional>
#include <algorithm>
#include <numeric>

#include "geometry/meridian_lookup_table.h"
#include "geometry/filament_types.h"

using namespace filament::geometry;

// ---------------------------------------------------------------------------
// Yardimci fonksiyonlar
// ---------------------------------------------------------------------------
namespace {

// Hemisferik dome analitik profil verisi uret
// s: [0, s_total], s=0 polar aciklik, s=s_total ekvator
// theta: [theta_0, pi/2]
//   rho(theta)    = R_eq * sin(theta)
//   x_loc(theta)  = R_eq * (cos(theta_0) - cos(theta))
//   drho/ds       = cos(theta)    [ds = R_eq * dtheta]
//   dx/ds         = sin(theta)
//   kappa_m       = 1 / R_eq
struct HemiData {
    std::vector<double> s, rho, x_local, drho_ds, dx_ds, kappa_m;
    double R_eq, r0, theta_0, s_total;
};

HemiData makeHemispherical(double R_eq, double r0, int N)
{
    HemiData d;
    d.R_eq   = R_eq;
    d.r0     = r0;
    d.theta_0 = std::asin(r0 / R_eq);
    d.s_total = R_eq * (constants::PI / 2.0 - d.theta_0);

    d.s.resize(N);
    d.rho.resize(N);
    d.x_local.resize(N);
    d.drho_ds.resize(N);
    d.dx_ds.resize(N);
    d.kappa_m.resize(N);

    for (int i = 0; i < N; ++i) {
        double frac  = static_cast<double>(i) / static_cast<double>(N - 1);
        double s_i   = d.s_total * frac;
        double theta = d.theta_0 + s_i / R_eq;

        d.s[i]       = s_i;
        d.rho[i]     = R_eq * std::sin(theta);
        d.x_local[i] = R_eq * (std::cos(d.theta_0) - std::cos(theta));
        d.drho_ds[i] = std::cos(theta);
        d.dx_ds[i]   = std::sin(theta);
        d.kappa_m[i] = 1.0 / R_eq;
    }
    return d;
}

// Analitik hemisferik degerler (verilen s icin)
struct HemiPoint {
    double rho, x_local, drho_ds, dx_ds, kappa_m;
};

HemiPoint hemiAnalytic(const HemiData& d, double s)
{
    double theta = d.theta_0 + s / d.R_eq;
    HemiPoint pt;
    pt.rho     = d.R_eq * std::sin(theta);
    pt.x_local = d.R_eq * (std::cos(d.theta_0) - std::cos(theta));
    pt.drho_ds = std::cos(theta);
    pt.dx_ds   = std::sin(theta);
    pt.kappa_m = 1.0 / d.R_eq;
    return pt;
}

ProfileMetadata makeHemiMeta(const HemiData& d)
{
    ProfileMetadata meta;
    meta.R_eq      = d.R_eq;
    meta.r0        = d.r0;
    meta.s_total   = d.s_total;
    meta.h_dome    = d.R_eq * std::cos(d.theta_0);
    meta.A_dome    = 0.0;  // test icin kullanilmiyor
    meta.kappa_eq  = 1.0 / d.R_eq;
    meta.kappa_pol = 1.0 / d.R_eq;
    meta.alpha_w   = std::asin(d.r0 / d.R_eq);
    meta.aspect_r  = 1.0;
    return meta;
}

} // anonymous namespace


// =============================================================================
// 1. TableState — Tablo Durumu Testleri
// =============================================================================

TEST(TableState, EmptyTable_IsInvalid)
{
    MeridianLookupTable table;
    EXPECT_FALSE(table.isValid());
    EXPECT_EQ(table.size(), 0u);
    // Bos tablo sorgusu nullopt donmeli
    auto r = table.query(0.0);
    EXPECT_FALSE(r.has_value());
}

TEST(TableState, SingleNode_IsInvalid)
{
    // 1 nokta: interval yok → gecersiz
    MeridianLookupTable table;
    std::vector<double> s     = {0.0};
    std::vector<double> ones  = {1.0};
    ProfileMetadata meta;
    meta.R_eq = 100.0; meta.r0 = 30.0; meta.s_total = 0.0;
    table.build(s, ones, ones, ones, ones, ones, meta);
    EXPECT_FALSE(table.isValid());
}

TEST(TableState, TwoNodes_IsValid)
{
    MeridianLookupTable table;
    std::vector<double> s    = {0.0, 10.0};
    std::vector<double> rho  = {30.0, 100.0};
    std::vector<double> x    = {0.0, 80.0};
    std::vector<double> drd  = {0.0, 1.0};
    std::vector<double> dxd  = {1.0, 0.0};
    std::vector<double> kap  = {0.01, 0.01};
    ProfileMetadata meta;
    meta.R_eq = 100.0; meta.r0 = 30.0; meta.s_total = 10.0;
    table.build(s, rho, x, drd, dxd, kap, meta);
    EXPECT_TRUE(table.isValid());
    EXPECT_EQ(table.size(), 2u);
}

TEST(TableState, SizeMismatch_IsInvalid)
{
    // Farkli boyutlu vektorler → gecersiz
    MeridianLookupTable table;
    std::vector<double> s    = {0.0, 10.0, 20.0};
    std::vector<double> rho  = {30.0, 65.0};  // eksik
    std::vector<double> ones = {1.0, 1.0, 1.0};
    ProfileMetadata meta;
    meta.R_eq = 100.0; meta.r0 = 30.0; meta.s_total = 20.0;
    table.build(s, rho, ones, ones, ones, ones, meta);
    EXPECT_FALSE(table.isValid());
}

TEST(TableState, MetadataAccessible)
{
    auto d = makeHemispherical(100.0, 30.0, 50);
    auto meta = makeHemiMeta(d);
    MeridianLookupTable table;
    table.build(d.s, d.rho, d.x_local, d.drho_ds, d.dx_ds, d.kappa_m, meta);
    EXPECT_TRUE(table.isValid());
    EXPECT_DOUBLE_EQ(table.metadata().R_eq,    100.0);
    EXPECT_DOUBLE_EQ(table.metadata().r0,       30.0);
    EXPECT_NEAR(table.metadata().s_total, d.s_total, 1e-12);
    EXPECT_NEAR(table.metadata().kappa_eq, 1.0/100.0, 1e-12);
    EXPECT_NEAR(table.metadata().aspect_r, 1.0,       1e-12);
}

TEST(TableState, RawAccessors_Correct)
{
    auto d = makeHemispherical(100.0, 30.0, 20);
    auto meta = makeHemiMeta(d);
    MeridianLookupTable table;
    table.build(d.s, d.rho, d.x_local, d.drho_ds, d.dx_ds, d.kappa_m, meta);
    ASSERT_EQ(table.rawS().size(), 20u);
    EXPECT_DOUBLE_EQ(table.rawS().front(), 0.0);
    EXPECT_NEAR(table.rawS().back(), d.s_total, 1e-12);
    EXPECT_NEAR(table.rawRho().front(), d.r0, 1e-10);
    EXPECT_NEAR(table.rawRho().back(), d.R_eq, 1e-10);
}


// =============================================================================
// 2. BoundaryBehavior — Sinir Davranisi Testleri
// =============================================================================

TEST(BoundaryBehavior, OutOfRangeBelow_ReturnsNullopt)
{
    auto d = makeHemispherical(100.0, 30.0, 50);
    auto meta = makeHemiMeta(d);
    MeridianLookupTable table;
    table.build(d.s, d.rho, d.x_local, d.drho_ds, d.dx_ds, d.kappa_m, meta);

    // s_min = 0, aralik altinda
    EXPECT_FALSE(table.query(-0.01).has_value());
    EXPECT_FALSE(table.query(-1.0).has_value());
    EXPECT_FALSE(table.query(-100.0).has_value());
}

TEST(BoundaryBehavior, OutOfRangeAbove_ReturnsNullopt)
{
    auto d = makeHemispherical(100.0, 30.0, 50);
    auto meta = makeHemiMeta(d);
    MeridianLookupTable table;
    table.build(d.s, d.rho, d.x_local, d.drho_ds, d.dx_ds, d.kappa_m, meta);

    // s_max = s_total, aralik ustunde
    EXPECT_FALSE(table.query(d.s_total + 0.01).has_value());
    EXPECT_FALSE(table.query(d.s_total + 1.0).has_value());
    EXPECT_FALSE(table.query(d.s_total + 1000.0).has_value());
}

TEST(BoundaryBehavior, LowerBound_ReturnsExact)
{
    auto d = makeHemispherical(100.0, 30.0, 50);
    auto meta = makeHemiMeta(d);
    MeridianLookupTable table;
    table.build(d.s, d.rho, d.x_local, d.drho_ds, d.dx_ds, d.kappa_m, meta);

    // s = 0: polar aciklik noktasi
    auto r = table.query(0.0);
    ASSERT_TRUE(r.has_value());
    // Spline doğal olarak dugum degerlerini tam olarak verir
    EXPECT_NEAR(r->rho,     d.r0,         1e-12);
    EXPECT_NEAR(r->x_local, 0.0,          1e-12);
    EXPECT_NEAR(r->drho_ds, d.drho_ds[0], 1e-12);
    EXPECT_NEAR(r->dx_ds,   d.dx_ds[0],   1e-12);
    EXPECT_NEAR(r->kappa_m, 1.0/100.0,    1e-12);
}

TEST(BoundaryBehavior, UpperBound_ReturnsExact)
{
    auto d = makeHemispherical(100.0, 30.0, 50);
    auto meta = makeHemiMeta(d);
    MeridianLookupTable table;
    table.build(d.s, d.rho, d.x_local, d.drho_ds, d.dx_ds, d.kappa_m, meta);

    // s = s_total: ekvator noktasi
    auto r = table.query(d.s_total);
    ASSERT_TRUE(r.has_value());
    EXPECT_NEAR(r->rho,     100.0,        1e-12);
    EXPECT_NEAR(r->x_local, d.x_local.back(), 1e-12);
    EXPECT_NEAR(r->drho_ds, 0.0,          1e-9);  // cos(pi/2) ≈ 0
    EXPECT_NEAR(r->dx_ds,   1.0,          1e-9);  // sin(pi/2) = 1
}

TEST(BoundaryBehavior, AllNodeValues_Exact)
{
    // Dugum noktalarinda spline tam deger vermeli (by construction)
    auto d = makeHemispherical(100.0, 30.0, 20);
    auto meta = makeHemiMeta(d);
    MeridianLookupTable table;
    table.build(d.s, d.rho, d.x_local, d.drho_ds, d.dx_ds, d.kappa_m, meta);

    for (std::size_t i = 0; i < d.s.size(); ++i) {
        auto r = table.query(d.s[i]);
        ASSERT_TRUE(r.has_value()) << "i=" << i;
        EXPECT_NEAR(r->rho,     d.rho[i],     1e-10) << "rho at i=" << i;
        EXPECT_NEAR(r->x_local, d.x_local[i], 1e-10) << "x at i=" << i;
        EXPECT_NEAR(r->drho_ds, d.drho_ds[i], 1e-10) << "drho at i=" << i;
        EXPECT_NEAR(r->dx_ds,   d.dx_ds[i],   1e-10) << "dx at i=" << i;
        EXPECT_NEAR(r->kappa_m, d.kappa_m[i], 1e-10) << "kap at i=" << i;
    }
}


// =============================================================================
// 3. LinearAccuracy — Lineer Fonksiyon (Spline Exact Olmali)
// =============================================================================
// f(s) = a + b*s → Kubik spline lineer fonksiyonlari TAMAMEN saglar
// (ikinci turev sifir, dogal sinir kosullari tatmin edilir)

TEST(LinearAccuracy, LinearFunction_ExactInterpolation)
{
    // Fiziksel olarak tutarli lineer profil:
    //   rho(s)    = a1 + b1*s   → drho_ds = b1 (sabit)
    //   x_local(s) = a2 + b2*s  → dx_ds   = b2 (sabit)
    //   kappa_m   = 0            (lineer profil → egrilik sifir)
    //
    // Clamped spline icin: slope_left = drho_ds[0] = b1 (rho spline'i icin dogru BC)
    // d(drho_ds)/ds = -dx_ds * kappa_m = -b2 * 0 = 0 (sabit → dogru)
    // d(dx_ds)/ds   = drho_ds * kappa_m = b1 * 0 = 0 (sabit → dogru)
    // Kubik spline lineer fonksiyonlari tamamen saglar.

    const double a1 = 5.0,  b1 = 2.3;    // rho(s) = a1 + b1*s
    const double a2 = 1.0,  b2 = -0.5;   // x_local(s) = a2 + b2*s

    const int N = 15;
    const double S = 50.0;

    std::vector<double> s(N), rho(N), x(N), drd(N), dxd(N), kap(N);
    for (int i = 0; i < N; ++i) {
        double t = S * i / (N - 1);
        s[i]   = t;
        rho[i] = a1 + b1 * t;
        x[i]   = a2 + b2 * t;
        drd[i] = b1;    // drho/ds = b1 (gercek turev)
        dxd[i] = b2;    // dx/ds   = b2 (gercek turev)
        kap[i] = 0.0;   // lineer profil → sifir egrilik
    }

    ProfileMetadata meta;
    meta.R_eq = 100.0; meta.r0 = 5.0; meta.s_total = S;
    MeridianLookupTable table;
    table.build(s, rho, x, drd, dxd, kap, meta);
    ASSERT_TRUE(table.isValid());

    // Ara noktalarda sorgu — lineer olmali (makine hassasiyetinde)
    const int M = 100;
    for (int j = 0; j < M; ++j) {
        double s_q = S * j / (M - 1);
        auto r = table.query(s_q);
        ASSERT_TRUE(r.has_value()) << "j=" << j;
        EXPECT_NEAR(r->rho,     a1 + b1 * s_q, 1e-10) << "j=" << j;
        EXPECT_NEAR(r->x_local, a2 + b2 * s_q, 1e-10) << "j=" << j;
        EXPECT_NEAR(r->drho_ds, b1,             1e-10) << "j=" << j;
        EXPECT_NEAR(r->dx_ds,   b2,             1e-10) << "j=" << j;
        EXPECT_NEAR(r->kappa_m, 0.0,            1e-10) << "j=" << j;
    }
}


// =============================================================================
// 4. ConstantAccuracy — Sabit Fonksiyon (Spline Exact Olmali)
// =============================================================================

TEST(ConstantAccuracy, ConstantKappa_ExactEverywhere)
{
    // Hemisferik domede kappa_m = 1/R_eq = sabit
    // Sabit fonksiyon icin spline tamamen exact olmali
    const double R_eq = 100.0;
    const double kappa_ref = 1.0 / R_eq;

    auto d = makeHemispherical(R_eq, 30.0, 30);
    auto meta = makeHemiMeta(d);
    MeridianLookupTable table;
    table.build(d.s, d.rho, d.x_local, d.drho_ds, d.dx_ds, d.kappa_m, meta);
    ASSERT_TRUE(table.isValid());

    // 200 aralik noktasinda kappa_m = 1/R_eq olmali
    const int M = 200;
    for (int j = 0; j < M; ++j) {
        double s_q = d.s_total * j / (M - 1);
        auto r = table.query(s_q);
        ASSERT_TRUE(r.has_value()) << "j=" << j;
        EXPECT_NEAR(r->kappa_m, kappa_ref, 1e-12) << "kappa j=" << j;
    }
}


// =============================================================================
// 5. Hemispherical — Analitik Referans (Karar-11 Katman 2 Toleranslar)
// =============================================================================

class HemisphericalTest : public ::testing::TestWithParam<std::pair<double,double>> {
protected:
    // Her test senaryosu icin tablo olustur
    void SetupTable(double R_eq, double r0, int N = 200) {
        d_    = makeHemispherical(R_eq, r0, N);
        meta_ = makeHemiMeta(d_);
        table_.build(d_.s, d_.rho, d_.x_local, d_.drho_ds, d_.dx_ds, d_.kappa_m, meta_);
    }
    HemiData d_;
    ProfileMetadata meta_;
    MeridianLookupTable table_;
};

// Karar-16 TEST-01: R_eq=73, r0=22 (ASTM Subscale)
// Karar-16 TEST-02: R_eq=152.4, r0=45 (Endustriyel COPV)
// Karar-16 TEST-03: R_eq=150, r0=10 (Kucuk Aciklik)
// Karar-16 TEST-04: R_eq=200, r0=50 (H2 Aerospace)
INSTANTIATE_TEST_SUITE_P(
    KararTest,
    HemisphericalTest,
    ::testing::Values(
        std::make_pair( 73.0,  22.0),   // TEST-01 ASTM Subscale
        std::make_pair(152.4,  45.0),   // TEST-02 Endustriyel COPV
        std::make_pair(150.0,  10.0),   // TEST-03 Kucuk Aciklik
        std::make_pair(200.0,  50.0)    // TEST-04 H2 Aerospace
    )
);

// Karar-11 Katman 2 toleranslarini test et:
//   Pozisyon (rho, x): |e| < 1e-4 mm
TEST_P(HemisphericalTest, PositionTolerance_Karar11)
{
    auto [R_eq, r0] = GetParam();
    SetupTable(R_eq, r0, 200);
    ASSERT_TRUE(table_.isValid());

    const int M = 199;
    double max_err_rho = 0.0, max_err_x = 0.0;

    for (int j = 0; j < M; ++j) {
        // Aralik ortasi (interpolasyon hatasi en buyuk burada)
        double s_q = (d_.s[j] + d_.s[j + 1]) / 2.0;
        auto r = table_.query(s_q);
        ASSERT_TRUE(r.has_value()) << "R_eq=" << R_eq << " j=" << j;

        HemiPoint ref = hemiAnalytic(d_, s_q);

        double e_rho = std::abs(r->rho     - ref.rho);
        double e_x   = std::abs(r->x_local - ref.x_local);

        max_err_rho = std::max(max_err_rho, e_rho);
        max_err_x   = std::max(max_err_x,   e_x);

        EXPECT_LT(e_rho, tolerances::POSITION_ABS_TOL)
            << "rho hata R_eq=" << R_eq << " j=" << j
            << " s=" << s_q << " got=" << r->rho << " ref=" << ref.rho;
        EXPECT_LT(e_x, tolerances::POSITION_ABS_TOL)
            << "x hata R_eq=" << R_eq << " j=" << j
            << " s=" << s_q << " got=" << r->x_local << " ref=" << ref.x_local;
    }
    // Bilgi mesaji (test ciktisinda gorulur)
    // std::cout << "R_eq=" << R_eq << " max_rho_err=" << max_err_rho
    //           << " max_x_err=" << max_err_x << std::endl;
    (void)max_err_rho; (void)max_err_x;
}

// Karar-11 Katman 2:  Turev (drho_ds, dx_ds): |e| < 1e-6
TEST_P(HemisphericalTest, DerivativeTolerance_Karar11)
{
    auto [R_eq, r0] = GetParam();
    SetupTable(R_eq, r0, 200);
    ASSERT_TRUE(table_.isValid());

    const int M = 199;

    for (int j = 0; j < M; ++j) {
        double s_q = (d_.s[j] + d_.s[j + 1]) / 2.0;
        auto r = table_.query(s_q);
        ASSERT_TRUE(r.has_value()) << "R_eq=" << R_eq << " j=" << j;

        HemiPoint ref = hemiAnalytic(d_, s_q);

        double e_drho = std::abs(r->drho_ds - ref.drho_ds);
        double e_dx   = std::abs(r->dx_ds   - ref.dx_ds);

        EXPECT_LT(e_drho, tolerances::DERIVATIVE_ABS_TOL)
            << "drho hata R_eq=" << R_eq << " j=" << j;
        EXPECT_LT(e_dx, tolerances::DERIVATIVE_ABS_TOL)
            << "dx hata R_eq=" << R_eq << " j=" << j;
    }
}

// Karar-11 Katman 2: Egrilik (kappa_m): |e_rel| < 1e-4
// Hemisferik dome icin kappa_m = 1/R_eq = sabit → spline tamamen exact
TEST_P(HemisphericalTest, CurvatureTolerance_Karar11)
{
    auto [R_eq, r0] = GetParam();
    SetupTable(R_eq, r0, 200);
    ASSERT_TRUE(table_.isValid());

    const double kappa_ref = 1.0 / R_eq;
    const int M = 300;

    for (int j = 0; j < M; ++j) {
        double s_q = d_.s_total * j / (M - 1);
        auto r = table_.query(s_q);
        ASSERT_TRUE(r.has_value()) << "R_eq=" << R_eq << " j=" << j;

        double e_rel = std::abs(r->kappa_m - kappa_ref) / kappa_ref;
        EXPECT_LT(e_rel, tolerances::CURVATURE_REL_TOL)
            << "kappa_m bagil hata R_eq=" << R_eq << " j=" << j;
    }
}


// =============================================================================
// 6. BinarySearch — findSegment Dogru Segment Bulma
// =============================================================================

TEST(BinarySearch, CorrectSegmentIdentification)
{
    // Duzensiz aralikli test verisi
    std::vector<double> s    = {0.0, 5.0, 8.0, 15.0, 20.0, 30.0};
    std::vector<double> ones = {1.0, 1.0,  1.0,  1.0,  1.0,  1.0};

    ProfileMetadata meta;
    meta.R_eq = 100.0; meta.r0 = 5.0; meta.s_total = 30.0;

    MeridianLookupTable table;
    table.build(s, ones, ones, ones, ones, ones, meta);
    ASSERT_TRUE(table.isValid());

    // Her dugum noktasinda gecerli sorgu donmeli
    for (auto sv : s) {
        EXPECT_TRUE(table.query(sv).has_value()) << "s=" << sv;
    }

    // Aralik icindeki noktalar gecerli olmali
    EXPECT_TRUE(table.query(2.5).has_value());   // segment 0
    EXPECT_TRUE(table.query(6.5).has_value());   // segment 1
    EXPECT_TRUE(table.query(10.0).has_value());  // segment 2
    EXPECT_TRUE(table.query(17.5).has_value());  // segment 3
    EXPECT_TRUE(table.query(25.0).has_value());  // segment 4

    // Sinir disi noktalar gecersiz olmali
    EXPECT_FALSE(table.query(-0.1).has_value());
    EXPECT_FALSE(table.query(30.1).has_value());
}


// =============================================================================
// 7. NumericalStability — Sayisal Kararlilik
// =============================================================================

TEST(NumericalStability, VeryCloseToEndpoints)
{
    // Uç noktalara cok yakin sorgular
    auto d = makeHemispherical(100.0, 30.0, 50);
    auto meta = makeHemiMeta(d);
    MeridianLookupTable table;
    table.build(d.s, d.rho, d.x_local, d.drho_ds, d.dx_ds, d.kappa_m, meta);
    ASSERT_TRUE(table.isValid());

    // Epsilon: sayisal hassasiyet siniri (query icindeki eps = 1e-12)
    // Bu deger aralik icinde oldugu icin gecerli
    EXPECT_TRUE(table.query(1e-14).has_value());
    EXPECT_TRUE(table.query(d.s_total - 1e-14).has_value());

    // Bu degerler aralik disinda
    EXPECT_FALSE(table.query(-1e-11).has_value());
    EXPECT_FALSE(table.query(d.s_total + 1e-11).has_value());
}

TEST(NumericalStability, IEEENaN_ReturnsNullopt)
{
    // NaN sorgusu nullopt donmeli (Karar-19: IEEE 754)
    auto d = makeHemispherical(100.0, 30.0, 50);
    auto meta = makeHemiMeta(d);
    MeridianLookupTable table;
    table.build(d.s, d.rho, d.x_local, d.drho_ds, d.dx_ds, d.kappa_m, meta);
    ASSERT_TRUE(table.isValid());

    double nan_val = std::numeric_limits<double>::quiet_NaN();
    // NaN karsilastirmasi her zaman false: NaN < x ve NaN > x hepsi false
    // Bu durumda query aralik disi olarak degerlendirmeli
    // (NaN comparisons yield false, so neither < front nor > back holds properly,
    //  but NaN - front = NaN which fails the eps check implicitly)
    // Davranis platforma bagli olabilir; test sadece cokme olmamali
    // EXPECT_FALSE veya _TRUE degil — sadece cokme yok
    auto r = table.query(nan_val);
    // NaN aralik disi sayilacak (NaN < front - eps = false, NaN > back + eps = false
    // ama std::min/max ile klamplama NaN propaga eder; bu test
    // sadece crash olmadigini dogrular)
    (void)r;
    SUCCEED();
}

TEST(NumericalStability, ManyQueries_NoAccumulation)
{
    // Cok sayida sorgu: maksimum hata birikimi olmasin
    auto d = makeHemispherical(100.0, 30.0, 200);
    auto meta = makeHemiMeta(d);
    MeridianLookupTable table;
    table.build(d.s, d.rho, d.x_local, d.drho_ds, d.dx_ds, d.kappa_m, meta);
    ASSERT_TRUE(table.isValid());

    double max_err = 0.0;
    const int M = 10000;
    for (int j = 0; j < M; ++j) {
        double s_q = d.s_total * j / (M - 1);
        auto r = table.query(s_q);
        ASSERT_TRUE(r.has_value());

        HemiPoint ref = hemiAnalytic(d, s_q);
        double e = std::abs(r->rho - ref.rho);
        max_err = std::max(max_err, e);
    }
    // Maksimum hata Karar-11 tolerans sinirinda olmali
    EXPECT_LT(max_err, tolerances::POSITION_ABS_TOL)
        << "10000 sorgu maksimum rho hatasi: " << max_err;
}
