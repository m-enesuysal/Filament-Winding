# =============================================================================
# filament_core — Statik Kütüphane (Karar-19)
# =============================================================================
# Phase-1b: geometry modülü
# Sonraki fazlarda: winding, kinematics, trajectory, postprocessor eklenir.
# =============================================================================

set(GEOMETRY_SOURCES
    meridian_lookup_table.cpp
    hemispherical_profile.cpp
    ellipsoidal_profile.cpp
    isotensoid_profile.cpp
    mandrel_geometry.cpp
    config_parser.cpp
    ${CMAKE_SOURCE_DIR}/src/geodesic/geodesic_solver.cpp
)

add_library(filament_core STATIC ${GEOMETRY_SOURCES})

target_include_directories(filament_core
    PUBLIC
        ${FILAMENT_INCLUDE_DIR}
)

target_link_libraries(filament_core
    PUBLIC
        Boost::boost
        Eigen3::Eigen
        nlohmann_json::nlohmann_json
)
