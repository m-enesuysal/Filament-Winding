// =============================================================================
// config_parser.cpp — Placeholder (Phase-1b S8'de implement edilecek)
// =============================================================================

#include "geometry/config_parser.h"
#include <nlohmann/json.hpp>
#include <stdexcept>

namespace filament {
namespace geometry {

WindingConfig loadConfig(const std::string& /*filepath*/)
{
    throw std::runtime_error(
        "loadConfig() henuz implement edilmedi. "
        "Phase-1b S8 session'inda yapilacak.");
}

} // namespace geometry
} // namespace filament
