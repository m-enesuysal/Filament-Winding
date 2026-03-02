# =============================================================================
# Filament Winding Yazılımı — Kök CMake Yapılandırması
# =============================================================================
#
# Proje: 4-Eksen Lathe-Type COPV Filament Sarım Makinesi Yazılımı
# Faz:   Phase-1b (Geometry Kernel)
# Referans: docs/phase0_decisions.md — Karar-19
#
# Targetlar:
#   filament_core   — Statik kütüphane (tüm numerik hesap çekirdeği)
#   filament_tests  — Google Test birim testleri
#   filament_python — pybind11 binding modülü (Phase-1b'de placeholder)
#
# KRİTİK KURAL (Karar-19):
#   -ffast-math HİÇBİR KOŞULDA KULLANILMAYACAK.
#   IEEE 754 uyumluluğu zorunludur. Gerekçe: singülarite tespiti (Karar-15)
#   ve tolerans yönetimi (Karar-11) bu flag ile güvenilir çalışamaz.
#   Bu kuralı değiştirmeye çalışmayın.
#
# =============================================================================

cmake_minimum_required(VERSION 3.16)

project(FilamentWinding
    VERSION   0.1.0
    LANGUAGES CXX
    DESCRIPTION "4-Eksen Lathe-Type COPV Filament Sarim Yazilimi"
)

# ---------------------------------------------------------------------------
# C++ Standardı (Karar-19: C++17 zorunlu)
# ---------------------------------------------------------------------------
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# ---------------------------------------------------------------------------
# Compiler Uyarıları (Karar-19: -Wall -Wextra -Wpedantic)
# ---------------------------------------------------------------------------
add_compile_options(
    $<$<CXX_COMPILER_ID:GNU,Clang,AppleClang>:-Wall>
    $<$<CXX_COMPILER_ID:GNU,Clang,AppleClang>:-Wextra>
    $<$<CXX_COMPILER_ID:GNU,Clang,AppleClang>:-Wpedantic>
)

# ---------------------------------------------------------------------------
# KRİTİK: -ffast-math YASAĞI (Karar-19)
# IEEE 754 uyumluluğu zorunlu. NaN/Inf kontrolü, singülarite tespiti
# (Karar-15) ve tolerans yönetimi (Karar-11) bu flag ile çalışamaz.
# ---------------------------------------------------------------------------
add_compile_options(
    $<$<CXX_COMPILER_ID:GNU,Clang,AppleClang>:-fno-fast-math>
)

# ---------------------------------------------------------------------------
# Build Tipi Varsayılanı
# ---------------------------------------------------------------------------
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
    set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS
        "Debug" "Release" "RelWithDebInfo" "MinSizeRel")
endif()

# Karar-19: Release -O2 (CMake varsayılan -O3 yerine)
set(CMAKE_CXX_FLAGS_RELEASE "-O2 -DNDEBUG")

# ---------------------------------------------------------------------------
# Bağımlılıklar (Karar-19)
# ---------------------------------------------------------------------------

# Boost (Odeint — header-only, Karar-6: ODE çözücü)
find_package(Boost REQUIRED)

# Eigen3 (Lineer cebir, Karar-19)
find_package(Eigen3 3.3 REQUIRED NO_MODULE)

# Google Test (Karar-19)
include(FetchContent)
FetchContent_Declare(
    googletest
    GIT_REPOSITORY https://github.com/google/googletest.git
    GIT_TAG        v1.14.0
)
set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
FetchContent_MakeAvailable(googletest)

# nlohmann/json (header-only, Karar-20)
FetchContent_Declare(
    nlohmann_json
    GIT_REPOSITORY https://github.com/nlohmann/json.git
    GIT_TAG        v3.11.3
)
FetchContent_MakeAvailable(nlohmann_json)

# ---------------------------------------------------------------------------
# Ortak include dizini
# ---------------------------------------------------------------------------
set(FILAMENT_INCLUDE_DIR ${CMAKE_SOURCE_DIR}/include)

# ---------------------------------------------------------------------------
# Alt dizinler
# ---------------------------------------------------------------------------
add_subdirectory(src/geometry)
add_subdirectory(tests)

# ---------------------------------------------------------------------------
# filament_python (Karar-18: pybind11 — Phase-1b'de placeholder)
# ---------------------------------------------------------------------------
option(BUILD_PYTHON_BINDINGS "pybind11 Python binding derle" OFF)
if(BUILD_PYTHON_BINDINGS)
    find_package(pybind11 REQUIRED)
    pybind11_add_module(filament_python python/visualization/bindings.cpp)
    target_link_libraries(filament_python PRIVATE filament_core)
endif()

# ---------------------------------------------------------------------------
# Build bilgisi
# ---------------------------------------------------------------------------
message(STATUS "")
message(STATUS "=== Filament Winding Build Yapilandirmasi ===")
message(STATUS "  C++ Standardi  : C++${CMAKE_CXX_STANDARD}")
message(STATUS "  Build Tipi     : ${CMAKE_BUILD_TYPE}")
message(STATUS "  Compiler       : ${CMAKE_CXX_COMPILER_ID} ${CMAKE_CXX_COMPILER_VERSION}")
message(STATUS "  Boost          : ${Boost_VERSION}")
message(STATUS "  Eigen3         : ${Eigen3_VERSION}")
message(STATUS "  ffast-math     : YASAK (Karar-19, IEEE 754 zorunlu)")
message(STATUS "  Python binding : ${BUILD_PYTHON_BINDINGS}")
message(STATUS "=============================================")
message(STATUS "")
