// =============================================================================
// test_config_parser.cpp — JSON Konfigurasyon Parser Birim Testleri
// =============================================================================
// Phase-1b S8: config_parser.h / config_parser.cpp testleri
// Karar-20: JSON format dogrulama
// Karar-8:  Derece → radyan donusumu dogrulama
// Karar-5:  Input validation kisitlari
// =============================================================================

#include <gtest/gtest.h>
#include "geometry/config_parser.h"
#include "geometry/filament_types.h"
#include <cmath>
#include <fstream>
#include <string>

using namespace filament::geometry;

// ===========================================================================
// Yardimci: Karar-20 ornek JSON sablon
// ===========================================================================
static const char* FULL_VALID_JSON = R"({
    "_units": "Uzunluk: mm, Aci: derece, Hiz: mm/s, Kuvvet: N",

    "mandrel": {
        "R_eq": 152.4,
        "r0": 45.0,
        "L_cyl": 300.0,
        "dome_type": "isotensoid",
        "k": null
    },

    "tow": {
        "BW": 6.35,
        "BT": 0.25,
        "N_tow": 1,
        "Fiber_tension": 50.0,
        "Winding_type": "wet"
    },

    "winding_sequence": [
        {"winding_type": "helical", "alpha_deg": 30.0, "N_layers": 2},
        {"winding_type": "hoop",    "alpha_deg": 89.0, "N_layers": 1},
        {"winding_type": "helical", "alpha_deg": 30.0, "N_layers": 2},
        {"winding_type": "hoop",    "alpha_deg": 89.0, "N_layers": 1}
    ]
})";


// ===========================================================================
// A. Gecerli Parse Testleri
// ===========================================================================

TEST(ConfigParserValid, FullValidJSON_Karar20Template)
{
    auto cfg = loadConfigFromString(FULL_VALID_JSON);

    // Mandrel
    EXPECT_DOUBLE_EQ(cfg.mandrel.R_eq, 152.4);
    EXPECT_DOUBLE_EQ(cfg.mandrel.r0, 45.0);
    EXPECT_DOUBLE_EQ(cfg.mandrel.L_cyl, 300.0);
    EXPECT_EQ(cfg.mandrel.dome_type, DomeType::Isotensoid);
    EXPECT_DOUBLE_EQ(cfg.mandrel.k, 1.0);  // null → varsayilan 1.0

    // Tow
    EXPECT_DOUBLE_EQ(cfg.tow.BW, 6.35);
    EXPECT_DOUBLE_EQ(cfg.tow.BT, 0.25);
    EXPECT_EQ(cfg.tow.N_tow, 1);
    EXPECT_DOUBLE_EQ(cfg.tow.Fiber_tension, 50.0);
    EXPECT_EQ(cfg.tow.Winding_type, "wet");

    // Winding sequence
    ASSERT_EQ(cfg.winding_sequence.size(), 4u);
    EXPECT_EQ(cfg.winding_sequence[0].winding_type, "helical");
    EXPECT_EQ(cfg.winding_sequence[1].winding_type, "hoop");
    EXPECT_EQ(cfg.winding_sequence[0].N_layers, 2);
    EXPECT_EQ(cfg.winding_sequence[1].N_layers, 1);
}

TEST(ConfigParserValid, HemisphericalDome)
{
    const char* json = R"({
        "mandrel": {
            "R_eq": 73.0,
            "r0": 22.0,
            "L_cyl": 200.0,
            "dome_type": "hemispherical"
        }
    })";

    auto cfg = loadConfigFromString(json);
    EXPECT_EQ(cfg.mandrel.dome_type, DomeType::Hemispherical);
    EXPECT_DOUBLE_EQ(cfg.mandrel.k, 1.0);  // varsayilan
    EXPECT_DOUBLE_EQ(cfg.mandrel.R_eq, 73.0);
    EXPECT_DOUBLE_EQ(cfg.mandrel.r0, 22.0);
    EXPECT_DOUBLE_EQ(cfg.mandrel.L_cyl, 200.0);
}

TEST(ConfigParserValid, EllipsoidalDomeWithK)
{
    const char* json = R"({
        "mandrel": {
            "R_eq": 152.4,
            "r0": 45.0,
            "L_cyl": 300.0,
            "dome_type": "ellipsoidal",
            "k": 0.7
        }
    })";

    auto cfg = loadConfigFromString(json);
    EXPECT_EQ(cfg.mandrel.dome_type, DomeType::Ellipsoidal);
    EXPECT_DOUBLE_EQ(cfg.mandrel.k, 0.7);
}

TEST(ConfigParserValid, IsotensoidDomeNullK)
{
    const char* json = R"({
        "mandrel": {
            "R_eq": 150.0,
            "r0": 10.0,
            "L_cyl": 0.0,
            "dome_type": "isotensoid",
            "k": null
        }
    })";

    auto cfg = loadConfigFromString(json);
    EXPECT_EQ(cfg.mandrel.dome_type, DomeType::Isotensoid);
    EXPECT_DOUBLE_EQ(cfg.mandrel.k, 1.0);     // null → 1.0
    EXPECT_DOUBLE_EQ(cfg.mandrel.L_cyl, 0.0); // L_cyl=0 gecerli
}

TEST(ConfigParserValid, MandrelOnlyNoTowNoSequence)
{
    const char* json = R"({
        "mandrel": {
            "R_eq": 200.0,
            "r0": 50.0,
            "L_cyl": 400.0,
            "dome_type": "hemispherical"
        }
    })";

    auto cfg = loadConfigFromString(json);
    EXPECT_DOUBLE_EQ(cfg.mandrel.R_eq, 200.0);
    EXPECT_TRUE(cfg.winding_sequence.empty());
    // Tow defaults
    EXPECT_DOUBLE_EQ(cfg.tow.BW, 0.0);
    EXPECT_DOUBLE_EQ(cfg.tow.BT, 0.0);
}

TEST(ConfigParserValid, DryWindingType)
{
    const char* json = R"({
        "mandrel": {
            "R_eq": 100.0, "r0": 30.0, "L_cyl": 150.0,
            "dome_type": "hemispherical"
        },
        "tow": {
            "BW": 3.175, "BT": 0.15, "N_tow": 4,
            "Fiber_tension": 80.0, "Winding_type": "dry"
        }
    })";

    auto cfg = loadConfigFromString(json);
    EXPECT_EQ(cfg.tow.Winding_type, "dry");
    EXPECT_EQ(cfg.tow.N_tow, 4);
}

TEST(ConfigParserValid, EllipsoidalMinK)
{
    // k = ELLIPSOIDAL_K_MIN sinirinda → gecerli olmali
    const char* json = R"({
        "mandrel": {
            "R_eq": 100.0, "r0": 30.0, "L_cyl": 100.0,
            "dome_type": "ellipsoidal",
            "k": 0.15
        }
    })";

    EXPECT_NO_THROW(loadConfigFromString(json));
}

TEST(ConfigParserValid, PolarWindingType)
{
    const char* json = R"({
        "mandrel": {
            "R_eq": 100.0, "r0": 30.0, "L_cyl": 100.0,
            "dome_type": "hemispherical"
        },
        "winding_sequence": [
            {"winding_type": "polar", "alpha_deg": 10.0, "N_layers": 1}
        ]
    })";

    auto cfg = loadConfigFromString(json);
    ASSERT_EQ(cfg.winding_sequence.size(), 1u);
    EXPECT_EQ(cfg.winding_sequence[0].winding_type, "polar");
}


// ===========================================================================
// B. Karar-8: Derece → Radyan Donusumu
// ===========================================================================

TEST(ConfigParserDegreeConversion, AlphaDegToRad_30)
{
    const char* json = R"({
        "mandrel": {
            "R_eq": 100.0, "r0": 30.0, "L_cyl": 100.0,
            "dome_type": "hemispherical"
        },
        "winding_sequence": [
            {"winding_type": "helical", "alpha_deg": 30.0, "N_layers": 1}
        ]
    })";

    auto cfg = loadConfigFromString(json);
    double expected_rad = 30.0 * constants::PI / 180.0;
    EXPECT_NEAR(cfg.winding_sequence[0].alpha_rad, expected_rad, 1e-15);
}

TEST(ConfigParserDegreeConversion, AlphaDegToRad_89)
{
    const char* json = R"({
        "mandrel": {
            "R_eq": 100.0, "r0": 30.0, "L_cyl": 100.0,
            "dome_type": "hemispherical"
        },
        "winding_sequence": [
            {"winding_type": "hoop", "alpha_deg": 89.0, "N_layers": 1}
        ]
    })";

    auto cfg = loadConfigFromString(json);
    double expected_rad = 89.0 * constants::PI / 180.0;
    EXPECT_NEAR(cfg.winding_sequence[0].alpha_rad, expected_rad, 1e-15);
}

TEST(ConfigParserDegreeConversion, AlphaDegToRad_45)
{
    const char* json = R"({
        "mandrel": {
            "R_eq": 100.0, "r0": 30.0, "L_cyl": 100.0,
            "dome_type": "hemispherical"
        },
        "winding_sequence": [
            {"winding_type": "helical", "alpha_deg": 45.0, "N_layers": 2}
        ]
    })";

    auto cfg = loadConfigFromString(json);
    double expected_rad = constants::PI / 4.0;
    EXPECT_NEAR(cfg.winding_sequence[0].alpha_rad, expected_rad, 1e-15);
}

TEST(ConfigParserDegreeConversion, MultipleLayersConversion)
{
    auto cfg = loadConfigFromString(FULL_VALID_JSON);

    // 30 derece → PI/6 rad
    double rad_30 = 30.0 * constants::PI / 180.0;
    double rad_89 = 89.0 * constants::PI / 180.0;

    ASSERT_EQ(cfg.winding_sequence.size(), 4u);
    EXPECT_NEAR(cfg.winding_sequence[0].alpha_rad, rad_30, 1e-15);
    EXPECT_NEAR(cfg.winding_sequence[1].alpha_rad, rad_89, 1e-15);
    EXPECT_NEAR(cfg.winding_sequence[2].alpha_rad, rad_30, 1e-15);
    EXPECT_NEAR(cfg.winding_sequence[3].alpha_rad, rad_89, 1e-15);
}


// ===========================================================================
// C. Mandrel Validation — Gecersiz Parametreler
// ===========================================================================

TEST(ConfigParserMandrelValidation, MissingMandrelSection)
{
    EXPECT_THROW(loadConfigFromString(R"({})"), std::invalid_argument);
}

TEST(ConfigParserMandrelValidation, MissingReq)
{
    const char* json = R"({
        "mandrel": {
            "r0": 22.0, "L_cyl": 200.0, "dome_type": "hemispherical"
        }
    })";
    EXPECT_THROW(loadConfigFromString(json), std::invalid_argument);
}

TEST(ConfigParserMandrelValidation, MissingR0)
{
    const char* json = R"({
        "mandrel": {
            "R_eq": 73.0, "L_cyl": 200.0, "dome_type": "hemispherical"
        }
    })";
    EXPECT_THROW(loadConfigFromString(json), std::invalid_argument);
}

TEST(ConfigParserMandrelValidation, MissingLcyl)
{
    const char* json = R"({
        "mandrel": {
            "R_eq": 73.0, "r0": 22.0, "dome_type": "hemispherical"
        }
    })";
    EXPECT_THROW(loadConfigFromString(json), std::invalid_argument);
}

TEST(ConfigParserMandrelValidation, MissingDomeType)
{
    const char* json = R"({
        "mandrel": {
            "R_eq": 73.0, "r0": 22.0, "L_cyl": 200.0
        }
    })";
    EXPECT_THROW(loadConfigFromString(json), std::invalid_argument);
}

TEST(ConfigParserMandrelValidation, InvalidDomeType)
{
    const char* json = R"({
        "mandrel": {
            "R_eq": 73.0, "r0": 22.0, "L_cyl": 200.0,
            "dome_type": "torispherical"
        }
    })";
    EXPECT_THROW(loadConfigFromString(json), std::invalid_argument);
}

TEST(ConfigParserMandrelValidation, ReqZero)
{
    const char* json = R"({
        "mandrel": {
            "R_eq": 0.0, "r0": 22.0, "L_cyl": 200.0,
            "dome_type": "hemispherical"
        }
    })";
    EXPECT_THROW(loadConfigFromString(json), std::invalid_argument);
}

TEST(ConfigParserMandrelValidation, ReqNegative)
{
    const char* json = R"({
        "mandrel": {
            "R_eq": -73.0, "r0": 22.0, "L_cyl": 200.0,
            "dome_type": "hemispherical"
        }
    })";
    EXPECT_THROW(loadConfigFromString(json), std::invalid_argument);
}

TEST(ConfigParserMandrelValidation, R0Zero)
{
    const char* json = R"({
        "mandrel": {
            "R_eq": 73.0, "r0": 0.0, "L_cyl": 200.0,
            "dome_type": "hemispherical"
        }
    })";
    EXPECT_THROW(loadConfigFromString(json), std::invalid_argument);
}

TEST(ConfigParserMandrelValidation, R0Negative)
{
    const char* json = R"({
        "mandrel": {
            "R_eq": 73.0, "r0": -5.0, "L_cyl": 200.0,
            "dome_type": "hemispherical"
        }
    })";
    EXPECT_THROW(loadConfigFromString(json), std::invalid_argument);
}

TEST(ConfigParserMandrelValidation, R0EqualReq)
{
    const char* json = R"({
        "mandrel": {
            "R_eq": 73.0, "r0": 73.0, "L_cyl": 200.0,
            "dome_type": "hemispherical"
        }
    })";
    EXPECT_THROW(loadConfigFromString(json), std::invalid_argument);
}

TEST(ConfigParserMandrelValidation, R0GreaterThanReq)
{
    const char* json = R"({
        "mandrel": {
            "R_eq": 73.0, "r0": 100.0, "L_cyl": 200.0,
            "dome_type": "hemispherical"
        }
    })";
    EXPECT_THROW(loadConfigFromString(json), std::invalid_argument);
}

TEST(ConfigParserMandrelValidation, LcylNegative)
{
    const char* json = R"({
        "mandrel": {
            "R_eq": 73.0, "r0": 22.0, "L_cyl": -10.0,
            "dome_type": "hemispherical"
        }
    })";
    EXPECT_THROW(loadConfigFromString(json), std::invalid_argument);
}

TEST(ConfigParserMandrelValidation, EllipsoidalKBelowMin)
{
    // k = 0.1 < ELLIPSOIDAL_K_MIN = 0.15 → S-GEO-03
    const char* json = R"({
        "mandrel": {
            "R_eq": 73.0, "r0": 22.0, "L_cyl": 200.0,
            "dome_type": "ellipsoidal",
            "k": 0.1
        }
    })";
    EXPECT_THROW(loadConfigFromString(json), std::invalid_argument);
}


// ===========================================================================
// D. Tow Validation
// ===========================================================================

TEST(ConfigParserTowValidation, InvalidWindingType)
{
    const char* json = R"({
        "mandrel": {
            "R_eq": 73.0, "r0": 22.0, "L_cyl": 200.0,
            "dome_type": "hemispherical"
        },
        "tow": {
            "BW": 6.35, "BT": 0.25, "N_tow": 1,
            "Fiber_tension": 50.0, "Winding_type": "prepreg"
        }
    })";
    EXPECT_THROW(loadConfigFromString(json), std::invalid_argument);
}


// ===========================================================================
// E. Winding Sequence Validation
// ===========================================================================

TEST(ConfigParserSequenceValidation, MissingWindingType)
{
    const char* json = R"({
        "mandrel": {
            "R_eq": 73.0, "r0": 22.0, "L_cyl": 200.0,
            "dome_type": "hemispherical"
        },
        "winding_sequence": [
            {"alpha_deg": 30.0, "N_layers": 2}
        ]
    })";
    EXPECT_THROW(loadConfigFromString(json), std::invalid_argument);
}

TEST(ConfigParserSequenceValidation, MissingAlphaDeg)
{
    const char* json = R"({
        "mandrel": {
            "R_eq": 73.0, "r0": 22.0, "L_cyl": 200.0,
            "dome_type": "hemispherical"
        },
        "winding_sequence": [
            {"winding_type": "helical", "N_layers": 2}
        ]
    })";
    EXPECT_THROW(loadConfigFromString(json), std::invalid_argument);
}

TEST(ConfigParserSequenceValidation, MissingNLayers)
{
    const char* json = R"({
        "mandrel": {
            "R_eq": 73.0, "r0": 22.0, "L_cyl": 200.0,
            "dome_type": "hemispherical"
        },
        "winding_sequence": [
            {"winding_type": "helical", "alpha_deg": 30.0}
        ]
    })";
    EXPECT_THROW(loadConfigFromString(json), std::invalid_argument);
}

TEST(ConfigParserSequenceValidation, NLayersZero)
{
    const char* json = R"({
        "mandrel": {
            "R_eq": 73.0, "r0": 22.0, "L_cyl": 200.0,
            "dome_type": "hemispherical"
        },
        "winding_sequence": [
            {"winding_type": "helical", "alpha_deg": 30.0, "N_layers": 0}
        ]
    })";
    EXPECT_THROW(loadConfigFromString(json), std::invalid_argument);
}

TEST(ConfigParserSequenceValidation, AlphaDegZero)
{
    const char* json = R"({
        "mandrel": {
            "R_eq": 73.0, "r0": 22.0, "L_cyl": 200.0,
            "dome_type": "hemispherical"
        },
        "winding_sequence": [
            {"winding_type": "helical", "alpha_deg": 0.0, "N_layers": 1}
        ]
    })";
    EXPECT_THROW(loadConfigFromString(json), std::invalid_argument);
}

TEST(ConfigParserSequenceValidation, AlphaDeg90)
{
    const char* json = R"({
        "mandrel": {
            "R_eq": 73.0, "r0": 22.0, "L_cyl": 200.0,
            "dome_type": "hemispherical"
        },
        "winding_sequence": [
            {"winding_type": "hoop", "alpha_deg": 90.0, "N_layers": 1}
        ]
    })";
    EXPECT_THROW(loadConfigFromString(json), std::invalid_argument);
}

TEST(ConfigParserSequenceValidation, InvalidSequenceWindingType)
{
    const char* json = R"({
        "mandrel": {
            "R_eq": 73.0, "r0": 22.0, "L_cyl": 200.0,
            "dome_type": "hemispherical"
        },
        "winding_sequence": [
            {"winding_type": "geodesic", "alpha_deg": 30.0, "N_layers": 1}
        ]
    })";
    EXPECT_THROW(loadConfigFromString(json), std::invalid_argument);
}


// ===========================================================================
// F. JSON Hata Durumlari
// ===========================================================================

TEST(ConfigParserJsonErrors, InvalidJsonSyntax)
{
    EXPECT_THROW(loadConfigFromString("{invalid json}"), std::runtime_error);
}

TEST(ConfigParserJsonErrors, EmptyString)
{
    EXPECT_THROW(loadConfigFromString(""), std::runtime_error);
}

TEST(ConfigParserJsonErrors, NonExistentFile)
{
    EXPECT_THROW(loadConfig("nonexistent_file_12345.json"), std::runtime_error);
}


// ===========================================================================
// G. Dosya Okuma Testi
// ===========================================================================

TEST(ConfigParserFileIO, WriteAndReadJsonFile)
{
    // Gecici JSON dosyasi olustur
    const std::string tmpfile = "test_config_temp.json";
    {
        std::ofstream f(tmpfile);
        f << R"({
            "mandrel": {
                "R_eq": 100.0,
                "r0": 25.0,
                "L_cyl": 180.0,
                "dome_type": "ellipsoidal",
                "k": 0.6
            },
            "tow": {
                "BW": 3.175,
                "BT": 0.2,
                "N_tow": 2,
                "Fiber_tension": 40.0,
                "Winding_type": "wet"
            },
            "winding_sequence": [
                {"winding_type": "helical", "alpha_deg": 25.0, "N_layers": 3}
            ]
        })";
    }

    auto cfg = loadConfig(tmpfile);
    EXPECT_DOUBLE_EQ(cfg.mandrel.R_eq, 100.0);
    EXPECT_DOUBLE_EQ(cfg.mandrel.r0, 25.0);
    EXPECT_DOUBLE_EQ(cfg.mandrel.L_cyl, 180.0);
    EXPECT_EQ(cfg.mandrel.dome_type, DomeType::Ellipsoidal);
    EXPECT_DOUBLE_EQ(cfg.mandrel.k, 0.6);

    EXPECT_EQ(cfg.tow.N_tow, 2);
    EXPECT_EQ(cfg.tow.Winding_type, "wet");

    ASSERT_EQ(cfg.winding_sequence.size(), 1u);
    double rad_25 = 25.0 * constants::PI / 180.0;
    EXPECT_NEAR(cfg.winding_sequence[0].alpha_rad, rad_25, 1e-15);

    // Temizle
    std::remove(tmpfile.c_str());
}
