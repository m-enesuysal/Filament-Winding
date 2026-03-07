// =============================================================================
// config_parser.h — JSON Konfigurasyon Dosya Parser
// =============================================================================
// Karar-20: JSON formati (nlohmann/json)
// Karar-8:  Derece → radyan donusumu burada (I/O katmani kurali)
// Karar-17: Ust katman — gecersiz dosya/parametre → exception
// =============================================================================

#ifndef FILAMENT_GEOMETRY_CONFIG_PARSER_H
#define FILAMENT_GEOMETRY_CONFIG_PARSER_H

#include "geometry/filament_types.h"
#include <string>
#include <vector>

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

struct WindingLayer {
    std::string winding_type;  // "helical", "hoop", "polar"
    double alpha_rad;          // Winding acisi [rad] (Karar-8: dahili radyan)
    int N_layers;
};

struct WindingConfig {
    MandrelConfig mandrel;
    TowConfig tow;
    std::vector<WindingLayer> winding_sequence;
};

// JSON dosyasindan konfigurasyon oku
// Ust katman (Karar-17): gecersiz dosya/parametre → exception
WindingConfig loadConfig(const std::string& filepath);

// JSON string'den konfigurasyon oku (test amacli)
WindingConfig loadConfigFromString(const std::string& json_str);

} // namespace geometry
} // namespace filament

#endif // FILAMENT_GEOMETRY_CONFIG_PARSER_H
