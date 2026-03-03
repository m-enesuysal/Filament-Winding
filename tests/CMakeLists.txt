# =============================================================================
# filament_tests — Google Test Birim Testleri (Karar-19)
# =============================================================================

enable_testing()

set(TEST_SOURCES
    test_placeholder.cpp
    test_meridian_lookup_table.cpp
)

add_executable(filament_tests ${TEST_SOURCES})

target_link_libraries(filament_tests
    PRIVATE
        filament_core
        GTest::gtest
        GTest::gtest_main
)

include(GoogleTest)
gtest_discover_tests(filament_tests)
