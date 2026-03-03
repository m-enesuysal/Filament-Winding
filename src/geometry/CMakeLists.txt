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
)

add_library(filament_core STATIC ${GEOMETRY_SOURCES})

target_include_directories(filament_core
    PUBLIC
        ${FILAMENT_INCLUDE_DIR}
)

target_link_libraries(filament_core
    PUBLIC
        Eigen3::Eigen
        nlohmann_json::nlohmann_json
)

# Boost hedef adı kurulum yöntemine göre farklı:
#   - System/vcpkg: Boost::boost  (FindBoost veya BoostConfig)
#   - FetchContent tarball (Boost 1.80+): Boost::headers
if(TARGET Boost::boost)
    target_link_libraries(filament_core PUBLIC Boost::boost)
elseif(TARGET Boost::headers)
    target_link_libraries(filament_core PUBLIC Boost::headers)
else()
    message(WARNING
        "Boost hedefi bulunamadı (Boost::boost veya Boost::headers). "
        "IsotensoidProfile (S5) Boost.Odeint gerektirdiğinden bu aşamada "
        "placeholder olarak derleniyor. Boost'u kurun veya CMAKE_PREFIX_PATH ayarlayın.")
endif()
