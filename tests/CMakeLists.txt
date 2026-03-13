# =============================================================================
# filament_tests — Google Test Birim Testleri (Karar-19)
# =============================================================================

enable_testing()

set(TEST_SOURCES
    test_placeholder.cpp
    test_meridian_lookup_table.cpp
    test_hemispherical_profile.cpp
    test_ellipsoidal_profile.cpp
    test_isotensoid_profile.cpp
    test_mandrel_geometry.cpp
    test_gate_1b_01.cpp
    test_config_parser.cpp
    test_geodesic_solver.cpp
)

# Referans veri dizini (Phase-2a CSV dosyalari)
set(REFERENCE_DATA_DIR "${CMAKE_SOURCE_DIR}/MATLAB/phase2a_winding/reference_data")

add_executable(filament_tests ${TEST_SOURCES})

target_link_libraries(filament_tests
    PRIVATE
        filament_core
        GTest::gtest
        GTest::gtest_main
)

target_compile_definitions(filament_tests
    PRIVATE
        REFERENCE_DATA_DIR="${REFERENCE_DATA_DIR}"
)

include(GoogleTest)
gtest_discover_tests(filament_tests)
