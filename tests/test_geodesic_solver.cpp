// =============================================================================
// test_geodesic_solver.cpp — GeodesicSolver Birim Testleri
// =============================================================================
// Phase-2a S4: Geodesic ODE cozucu dogrulamasi
//
// MATLAB referans CSV dosyalari ile delta_phi karsilastirmasi:
//   |delta_phi_cpp - delta_phi_matlab| < 1e-4 rad
//
// 3 dome tipi x S1 senaryosu = 3 test
// S1 parametreleri: R_eq=152.4, r0=45, L_cyl=300, BW_eff=10, k_ell=0.7
// =============================================================================

#include <gtest/gtest.h>
#include "geodesic/geodesic_solver.h"
#include "geometry/filament_types.h"
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <cmath>

using namespace filament;
using namespace filament::geodesic;

namespace {

// S1-COPV-Standart parametreleri (verify_geodesic_path.m ile ayni)
constexpr double S1_R_EQ   = 152.4;
constexpr double S1_R0     = 45.0;
constexpr double S1_L_CYL  = 300.0;
constexpr double S1_BW_EFF = 10.0;
constexpr double S1_K_ELL  = 0.7;
constexpr double S1_R_E    = S1_R0 + S1_BW_EFF / 2.0;  // = 50.0

constexpr double DELTA_PHI_TOL = 1e-4;  // rad

// CSV referans dosyasindan son satirdaki phi degerini oku (= delta_phi_circuit)
double readDeltaPhiFromCSV(const std::string& filepath)
{
    std::ifstream file(filepath);
    if (!file.is_open()) {
        throw std::runtime_error("CSV dosyasi acilamadi: " + filepath);
    }

    std::string line;
    std::string last_line;

    // Header atla
    std::getline(file, line);

    // Son satiri bul
    while (std::getline(file, line)) {
        if (!line.empty()) {
            last_line = line;
        }
    }

    if (last_line.empty()) {
        throw std::runtime_error("CSV dosyasi bos: " + filepath);
    }

    // CSV sutunlari: s,rho,x,phi,alpha,kn
    // phi = 4. sutun (0-indexed: 3)
    std::istringstream ss(last_line);
    std::string token;
    for (int col = 0; col < 4; ++col) {
        if (!std::getline(ss, token, ',')) {
            throw std::runtime_error("CSV parse hatasi: " + filepath);
        }
    }
    return std::stod(token);
}

std::string refDataPath(const std::string& filename)
{
#ifdef REFERENCE_DATA_DIR
    return std::string(REFERENCE_DATA_DIR) + "/" + filename;
#else
    return "MATLAB/phase2a_winding/reference_data/" + filename;
#endif
}

} // anonymous namespace

// =============================================================================
// Test 1: Hemispherical S1 — delta_phi dogrulamasi
// =============================================================================
TEST(GeodesicSolverTest, HemisphericalS1_DeltaPhi)
{
    // MATLAB referans delta_phi oku
    const std::string csv_path = refDataPath("geodesic_ref_hemispherical_S1.csv");
    const double delta_phi_matlab = readDeltaPhiFromCSV(csv_path);

    // GeodesicSolver ile hesapla
    GeodesicParams params;
    params.dome_type = geometry::DomeType::Hemispherical;
    params.R_eq      = S1_R_EQ;
    params.r0        = S1_R0;
    params.L_cyl     = S1_L_CYL;
    params.R_E       = S1_R_E;
    params.k         = 1.0;  // hemispherical icin kullanilmaz

    auto result = GeodesicSolver::solve(params);
    ASSERT_TRUE(result.has_value()) << "GeodesicSolver hemispherical S1 basarisiz";

    const double delta_phi_cpp = result->delta_phi;
    const double err = std::abs(delta_phi_cpp - delta_phi_matlab);

    // Diagnostik bilgi
    std::printf("[HEMI S1] delta_phi_cpp=%.8f, delta_phi_matlab=%.8f, fark=%.2e rad\n",
                delta_phi_cpp, delta_phi_matlab, err);
    std::printf("          phi_dome=%.6f, phi_cyl=%.6f, phi_tail=%.2e, phi_num=%.6f\n",
                result->phi_dome, result->phi_cyl, result->phi_tail, result->phi_numerical);

    EXPECT_LT(err, DELTA_PHI_TOL)
        << "Hemispherical S1 delta_phi hatasi: " << err << " rad (limit: " << DELTA_PHI_TOL << ")";
}

// =============================================================================
// Test 2: Ellipsoidal S1 — delta_phi dogrulamasi
// =============================================================================
TEST(GeodesicSolverTest, EllipsoidalS1_DeltaPhi)
{
    const std::string csv_path = refDataPath("geodesic_ref_ellipsoidal_S1.csv");
    const double delta_phi_matlab = readDeltaPhiFromCSV(csv_path);

    GeodesicParams params;
    params.dome_type = geometry::DomeType::Ellipsoidal;
    params.R_eq      = S1_R_EQ;
    params.r0        = S1_R0;
    params.L_cyl     = S1_L_CYL;
    params.R_E       = S1_R_E;
    params.k         = S1_K_ELL;

    auto result = GeodesicSolver::solve(params);
    ASSERT_TRUE(result.has_value()) << "GeodesicSolver ellipsoidal S1 basarisiz";

    const double delta_phi_cpp = result->delta_phi;
    const double err = std::abs(delta_phi_cpp - delta_phi_matlab);

    std::printf("[ELLIP S1] delta_phi_cpp=%.8f, delta_phi_matlab=%.8f, fark=%.2e rad\n",
                delta_phi_cpp, delta_phi_matlab, err);
    std::printf("           phi_dome=%.6f, phi_cyl=%.6f, phi_tail=%.2e, phi_num=%.6f\n",
                result->phi_dome, result->phi_cyl, result->phi_tail, result->phi_numerical);

    EXPECT_LT(err, DELTA_PHI_TOL)
        << "Ellipsoidal S1 delta_phi hatasi: " << err << " rad (limit: " << DELTA_PHI_TOL << ")";
}

// =============================================================================
// Test 3: Isotensoid S1 — delta_phi dogrulamasi
// =============================================================================
TEST(GeodesicSolverTest, IsotensoidS1_DeltaPhi)
{
    const std::string csv_path = refDataPath("geodesic_ref_isotensoid_S1.csv");
    const double delta_phi_matlab = readDeltaPhiFromCSV(csv_path);

    GeodesicParams params;
    params.dome_type = geometry::DomeType::Isotensoid;
    params.R_eq      = S1_R_EQ;
    params.r0        = S1_R0;
    params.L_cyl     = S1_L_CYL;
    params.R_E       = S1_R_E;
    params.k         = 1.0;  // isotensoid icin kullanilmaz

    auto result = GeodesicSolver::solve(params);
    ASSERT_TRUE(result.has_value()) << "GeodesicSolver isotensoid S1 basarisiz";

    const double delta_phi_cpp = result->delta_phi;
    const double err = std::abs(delta_phi_cpp - delta_phi_matlab);

    std::printf("[ISO S1] delta_phi_cpp=%.8f, delta_phi_matlab=%.8f, fark=%.2e rad\n",
                delta_phi_cpp, delta_phi_matlab, err);
    std::printf("         phi_dome=%.6f, phi_cyl=%.6f, phi_tail=%.2e, phi_num=%.6f\n",
                result->phi_dome, result->phi_cyl, result->phi_tail, result->phi_numerical);

    EXPECT_LT(err, DELTA_PHI_TOL)
        << "Isotensoid S1 delta_phi hatasi: " << err << " rad (limit: " << DELTA_PHI_TOL << ")";
}

// =============================================================================
// Test 4: Gecersiz girisler — nullopt donmeli
// =============================================================================
TEST(GeodesicSolverTest, InvalidInputs)
{
    GeodesicParams params;
    params.dome_type = geometry::DomeType::Hemispherical;
    params.R_eq = 100.0;
    params.r0 = 30.0;
    params.L_cyl = 200.0;

    // R_E >= R_eq
    params.R_E = 100.0;
    EXPECT_FALSE(GeodesicSolver::solve(params).has_value());

    // R_E <= 0
    params.R_E = -5.0;
    EXPECT_FALSE(GeodesicSolver::solve(params).has_value());

    // R_eq <= r0
    params.R_E = 35.0;
    params.R_eq = 30.0;
    EXPECT_FALSE(GeodesicSolver::solve(params).has_value());
}

// =============================================================================
// Test 5: Temel fiziksel tutarlilik kontrolleri
// =============================================================================
TEST(GeodesicSolverTest, PhysicalConsistency)
{
    GeodesicParams params;
    params.dome_type = geometry::DomeType::Ellipsoidal;
    params.R_eq      = S1_R_EQ;
    params.r0        = S1_R0;
    params.L_cyl     = S1_L_CYL;
    params.R_E       = S1_R_E;
    params.k         = S1_K_ELL;

    auto result = GeodesicSolver::solve(params);
    ASSERT_TRUE(result.has_value());

    // delta_phi = 4*phi_dome + 2*phi_cyl
    const double expected = 4.0 * result->phi_dome + 2.0 * result->phi_cyl;
    EXPECT_NEAR(result->delta_phi, expected, 1e-12);

    // phi_dome = phi_numerical + phi_tail
    const double dome_sum = result->phi_numerical + result->phi_tail;
    EXPECT_NEAR(result->phi_dome, dome_sum, 1e-12);

    // phi_dome > 0, phi_cyl > 0
    EXPECT_GT(result->phi_dome, 0.0);
    EXPECT_GT(result->phi_cyl, 0.0);

    // alpha_eq = asin(R_E / R_eq)
    EXPECT_NEAR(result->alpha_eq, std::asin(S1_R_E / S1_R_EQ), 1e-12);

    // Yol vektorleri bos degil
    EXPECT_GT(result->s.size(), 10u);
    EXPECT_EQ(result->s.size(), result->rho.size());
    EXPECT_EQ(result->s.size(), result->phi.size());
    EXPECT_EQ(result->s.size(), result->alpha.size());

    // Ilk nokta: s=0, phi=0, rho=R_eq
    EXPECT_NEAR(result->s.front(), 0.0, 1e-12);
    EXPECT_NEAR(result->phi.front(), 0.0, 1e-12);
    EXPECT_NEAR(result->rho.front(), S1_R_EQ, 0.1);  // spline boundary

    // Son nokta: rho yaklasik R_E + epsilon
    EXPECT_NEAR(result->rho.back(), S1_R_E + 1e-3, 0.1);

    // alpha monoton artmali (dome'da ekvatordan turnaround'a)
    for (std::size_t i = 1; i < result->alpha.size(); ++i) {
        EXPECT_GE(result->alpha[i], result->alpha[i-1] - 1e-10);
    }
}

// =============================================================================
// Test 6: PCHIP Resample — nokta sayisi
// =============================================================================
TEST(GeodesicSolverResampleTest, PointCount)
{
    GeodesicParams params;
    params.dome_type = geometry::DomeType::Ellipsoidal;
    params.R_eq      = S1_R_EQ;
    params.r0        = S1_R0;
    params.L_cyl     = S1_L_CYL;
    params.R_E       = S1_R_E;
    params.k         = S1_K_ELL;

    auto path = GeodesicSolver::solve(params);
    ASSERT_TRUE(path.has_value());

    constexpr double DELTA_S = 0.5;  // mm
    auto resampled = resample(*path, DELTA_S);

    // Beklenen nokta sayisi: floor(range / delta_s) + 1 (+ olasi son nokta)
    const double range = path->s.back() - path->s.front();
    const std::size_t n_expected = static_cast<std::size_t>(std::floor(range / DELTA_S)) + 1;

    std::printf("[RESAMPLE] Orijinal: %zu nokta, Resampled: %zu nokta, Beklenen: ~%zu\n",
                path->s.size(), resampled.s.size(), n_expected);

    // Nokta sayisi makul aralikta: n_expected veya n_expected+1 (son nokta)
    EXPECT_GE(resampled.s.size(), n_expected);
    EXPECT_LE(resampled.s.size(), n_expected + 1);

    // Tum vektorler ayni boyut
    EXPECT_EQ(resampled.s.size(), resampled.rho.size());
    EXPECT_EQ(resampled.s.size(), resampled.x.size());
    EXPECT_EQ(resampled.s.size(), resampled.phi.size());
    EXPECT_EQ(resampled.s.size(), resampled.alpha.size());

    // Esit aralikli (ilk ve son haric)
    for (std::size_t i = 1; i < resampled.s.size() - 1; ++i) {
        const double ds = resampled.s[i] - resampled.s[i - 1];
        EXPECT_NEAR(ds, DELTA_S, 1e-10)
            << "Esit aralik bozuldu: i=" << i << ", ds=" << ds;
    }
}

// =============================================================================
// Test 7: PCHIP Resample — Clairaut korunumu
// =============================================================================
TEST(GeodesicSolverResampleTest, ClairautConservation)
{
    GeodesicParams params;
    params.dome_type = geometry::DomeType::Ellipsoidal;
    params.R_eq      = S1_R_EQ;
    params.r0        = S1_R0;
    params.L_cyl     = S1_L_CYL;
    params.R_E       = S1_R_E;
    params.k         = S1_K_ELL;

    auto path = GeodesicSolver::solve(params);
    ASSERT_TRUE(path.has_value());

    constexpr double DELTA_S = 0.5;
    auto resampled = resample(*path, DELTA_S);

    // Clairaut korunumu: |rho * sin(alpha) - R_E| < 1e-4 mm
    double max_clairaut_err = 0.0;
    for (std::size_t i = 0; i < resampled.s.size(); ++i) {
        const double clairaut = resampled.rho[i] * std::sin(resampled.alpha[i]);
        const double err = std::abs(clairaut - S1_R_E);
        if (err > max_clairaut_err) max_clairaut_err = err;
    }

    std::printf("[CLAIRAUT] Maks |rho*sin(alpha) - R_E| = %.2e mm (limit: 1e-4)\n",
                max_clairaut_err);

    EXPECT_LT(max_clairaut_err, 1e-4)
        << "Clairaut korunumu ihlali: " << max_clairaut_err << " mm";
}

// =============================================================================
// Test 8: PCHIP Resample — delta_phi korunumu
// =============================================================================
TEST(GeodesicSolverResampleTest, DeltaPhiPreserved)
{
    GeodesicParams params;
    params.dome_type = geometry::DomeType::Ellipsoidal;
    params.R_eq      = S1_R_EQ;
    params.r0        = S1_R0;
    params.L_cyl     = S1_L_CYL;
    params.R_E       = S1_R_E;
    params.k         = S1_K_ELL;

    auto path = GeodesicSolver::solve(params);
    ASSERT_TRUE(path.has_value());

    constexpr double DELTA_S = 0.5;
    auto resampled = resample(*path, DELTA_S);

    // Skaler delta_phi degismemeli (birebir kopyalanir)
    EXPECT_NEAR(resampled.delta_phi, path->delta_phi, 1e-15)
        << "delta_phi kopyalama hatasi";

    // phi_dome, phi_cyl de korunmali
    EXPECT_NEAR(resampled.phi_dome, path->phi_dome, 1e-15);
    EXPECT_NEAR(resampled.phi_cyl, path->phi_cyl, 1e-15);

    // Son phi degeri yaklasik ayni olmali (PCHIP interpolasyon farki < 1e-6)
    const double phi_end_orig = path->phi.back();
    const double phi_end_res  = resampled.phi.back();
    const double phi_diff = std::abs(phi_end_res - phi_end_orig);

    std::printf("[PHI END] Orijinal=%.8f, Resampled=%.8f, Fark=%.2e rad\n",
                phi_end_orig, phi_end_res, phi_diff);

    EXPECT_LT(phi_diff, 1e-6)
        << "Son phi degeri farki: " << phi_diff << " rad (limit: 1e-6)";
}

// =============================================================================
// GATE-2a Cross-Validation: 3 dome x 4 senaryo (mevcut CSV'ler)
// =============================================================================
struct ScenarioConfig {
    const char* name;
    const char* dome_name;
    geometry::DomeType dome_type;
    double R_eq, r0, L_cyl, BW_eff, k;
    const char* csv_file;  // nullptr ise CSV yok, sadece hesapla
};

class GeodesicCrossValidation : public ::testing::TestWithParam<ScenarioConfig> {};

TEST_P(GeodesicCrossValidation, DeltaPhiMatch)
{
    const auto& cfg = GetParam();
    const double R_E = cfg.r0 + cfg.BW_eff / 2.0;

    GeodesicParams params;
    params.dome_type = cfg.dome_type;
    params.R_eq = cfg.R_eq;
    params.r0   = cfg.r0;
    params.L_cyl = cfg.L_cyl;
    params.R_E   = R_E;
    params.k     = cfg.k;

    auto result = GeodesicSolver::solve(params);
    ASSERT_TRUE(result.has_value())
        << cfg.dome_name << " " << cfg.name << " basarisiz";

    if (cfg.csv_file) {
        const std::string csv_path = refDataPath(cfg.csv_file);
        const double delta_phi_matlab = readDeltaPhiFromCSV(csv_path);
        const double err = std::abs(result->delta_phi - delta_phi_matlab);

        std::printf("| %-14s | %-4s | %12.8f | %12.8f | %9.2e |\n",
                    cfg.dome_name, cfg.name, result->delta_phi, delta_phi_matlab, err);

        EXPECT_LT(err, DELTA_PHI_TOL)
            << cfg.dome_name << " " << cfg.name
            << " delta_phi hatasi: " << err << " rad";
    } else {
        std::printf("| %-14s | %-4s | %12.8f | %12s | %9s |\n",
                    cfg.dome_name, cfg.name, result->delta_phi, "N/A", "N/A");
    }
}

static const ScenarioConfig kAllScenarios[] = {
    // Hemispherical
    {"S1", "hemispherical",  geometry::DomeType::Hemispherical, 152.4, 45, 300, 10, 1.0,
     "geodesic_ref_hemispherical_S1.csv"},
    {"S2", "hemispherical",  geometry::DomeType::Hemispherical, 200,   60, 400, 12, 1.0,
     "geodesic_ref_hemispherical_S2.csv"},
    {"S3", "hemispherical",  geometry::DomeType::Hemispherical, 100,   20, 150,  6, 1.0,
     "geodesic_ref_hemispherical_S3.csv"},
    {"S4", "hemispherical",  geometry::DomeType::Hemispherical, 120,   50, 250,  8, 1.0,
     "geodesic_ref_hemispherical_S4.csv"},
    // Ellipsoidal
    {"S1", "ellipsoidal",    geometry::DomeType::Ellipsoidal,   152.4, 45, 300, 10, 0.7,
     "geodesic_ref_ellipsoidal_S1.csv"},
    {"S2", "ellipsoidal",    geometry::DomeType::Ellipsoidal,   200,   60, 400, 12, 0.5,
     "geodesic_ref_ellipsoidal_S2.csv"},
    {"S3", "ellipsoidal",    geometry::DomeType::Ellipsoidal,   100,   20, 150,  6, 1.2,
     "geodesic_ref_ellipsoidal_S3.csv"},
    {"S4", "ellipsoidal",    geometry::DomeType::Ellipsoidal,   120,   50, 250,  8, 0.8,
     "geodesic_ref_ellipsoidal_S4.csv"},
    // Isotensoid
    {"S1", "isotensoid",     geometry::DomeType::Isotensoid,    152.4, 45, 300, 10, 1.0,
     "geodesic_ref_isotensoid_S1.csv"},
    {"S2", "isotensoid",     geometry::DomeType::Isotensoid,    200,   60, 400, 12, 1.0,
     "geodesic_ref_isotensoid_S2.csv"},
    {"S3", "isotensoid",     geometry::DomeType::Isotensoid,    100,   20, 150,  6, 1.0,
     "geodesic_ref_isotensoid_S3.csv"},
    {"S4", "isotensoid",     geometry::DomeType::Isotensoid,    120,   50, 250,  8, 1.0,
     "geodesic_ref_isotensoid_S4.csv"},
};

INSTANTIATE_TEST_SUITE_P(
    GATE2a,
    GeodesicCrossValidation,
    ::testing::ValuesIn(kAllScenarios),
    [](const ::testing::TestParamInfo<ScenarioConfig>& info) {
        return std::string(info.param.dome_name) + "_" + info.param.name;
    }
);
