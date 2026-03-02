// =============================================================================
// config_parser.h — JSON Konfigurasyon Dosya Parser
// =============================================================================
// Karar-20: JSON formati (nlohmann/json)
// Karar-8:  Derece → radyan donusumu burada (I/O katmani kurali)
// Detayli implementasyon: Phase-1b S8
// =============================================================================

#ifndef FILAMENT_GEOMETRY_CONFIG_PARSER_H
#define FILAMENT_GEOMETRY_CONFIG_PARSER_H

#include "geometry/filament_types.h"
#include <string>

namespace filament {
namespace geometry {

struct MandrelConfig {
    double R_eq;
    double r0;
    double L_cyl;
    DomeType dome_type;
    double k;
};

struct TowConfig {
    double BW;
    double BT;
    int N_tow;
    double Fiber_tension;
    std::string Winding_type;
};

struct WindingConfig {
    MandrelConfig mandrel;
    TowConfig tow;
};

// JSON dosyasindan konfigurasyon oku
// Ust katman (Karar-17): gecersiz dosya/parametre → exception
WindingConfig loadConfig(const std::string& filepath);

} // namespace geometry
} // namespace filament

#endif // FILAMENT_GEOMETRY_CONFIG_PARSER_H
