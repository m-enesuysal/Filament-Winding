// =============================================================================
// config_parser.cpp — JSON Konfigurasyon Parser Implementasyonu
// =============================================================================
// Karar-20: JSON formati (nlohmann/json)
// Karar-8:  Derece → radyan donusumu I/O katmaninda
// Karar-5:  Input validation kisitlari
// Karar-17: Ust katman — exception-based hata yonetimi
// =============================================================================

#include "geometry/config_parser.h"
#include <nlohmann/json.hpp>
#include <fstream>
#include <stdexcept>

namespace filament {
namespace geometry {

// ---------------------------------------------------------------------------
// Yardimci: dome_type string → DomeType enum
// ---------------------------------------------------------------------------
static DomeType parseDomeType(const std::string& s)
{
    if (s == "hemispherical") return DomeType::Hemispherical;
    if (s == "ellipsoidal")   return DomeType::Ellipsoidal;
    if (s == "isotensoid")    return DomeType::Isotensoid;

    throw std::invalid_argument(
        "Gecersiz dome_type: \"" + s + "\". "
        "Gecerli degerler: \"hemispherical\", \"ellipsoidal\", \"isotensoid\".");
}

// ---------------------------------------------------------------------------
// Yardimci: winding_type string dogrulama
// ---------------------------------------------------------------------------
static void validateWindingType(const std::string& s)
{
    if (s == "helical" || s == "hoop" || s == "polar") return;

    throw std::invalid_argument(
        "Gecersiz winding_type: \"" + s + "\". "
        "Gecerli degerler: \"helical\", \"hoop\", \"polar\".");
}

// ---------------------------------------------------------------------------
// Yardimci: Winding_type (tow) string dogrulama
// ---------------------------------------------------------------------------
static void validateTowWindingType(const std::string& s)
{
    if (s == "wet" || s == "dry") return;

    throw std::invalid_argument(
        "Gecersiz Winding_type: \"" + s + "\". "
        "Gecerli degerler: \"wet\", \"dry\".");
}

// ---------------------------------------------------------------------------
// Ortak parse fonksiyonu (JSON nesnesinden)
// ---------------------------------------------------------------------------
static WindingConfig parseJson(const nlohmann::json& j)
{
    WindingConfig config{};

    // -----------------------------------------------------------------------
    // mandrel bolumu (zorunlu)
    // -----------------------------------------------------------------------
    if (!j.contains("mandrel")) {
        throw std::invalid_argument("JSON'da 'mandrel' bolumu eksik.");
    }
    const auto& jm = j["mandrel"];

    // Zorunlu alanlar
    if (!jm.contains("R_eq"))
        throw std::invalid_argument("mandrel.R_eq alani eksik.");
    if (!jm.contains("r0"))
        throw std::invalid_argument("mandrel.r0 alani eksik.");
    if (!jm.contains("L_cyl"))
        throw std::invalid_argument("mandrel.L_cyl alani eksik.");
    if (!jm.contains("dome_type"))
        throw std::invalid_argument("mandrel.dome_type alani eksik.");

    config.mandrel.R_eq = jm["R_eq"].get<double>();
    config.mandrel.r0   = jm["r0"].get<double>();
    config.mandrel.L_cyl = jm["L_cyl"].get<double>();
    config.mandrel.dome_type = parseDomeType(jm["dome_type"].get<std::string>());

    // k alani: null veya eksik → varsayilan 1.0
    if (jm.contains("k") && !jm["k"].is_null()) {
        config.mandrel.k = jm["k"].get<double>();
    } else {
        config.mandrel.k = 1.0;
    }

    // Karar-5 input validation
    if (config.mandrel.R_eq <= 0.0)
        throw std::invalid_argument("mandrel.R_eq pozitif olmali (R_eq > 0).");
    if (config.mandrel.r0 <= 0.0)
        throw std::invalid_argument("mandrel.r0 pozitif olmali (r0 > 0).");
    if (config.mandrel.r0 >= config.mandrel.R_eq)
        throw std::invalid_argument("mandrel.r0, R_eq'den kucuk olmali (r0 < R_eq).");
    if (config.mandrel.L_cyl < 0.0)
        throw std::invalid_argument("mandrel.L_cyl negatif olamaz (L_cyl >= 0).");

    // Elipsoidal icin k kontrolu (Karar-15 S-GEO-03)
    if (config.mandrel.dome_type == DomeType::Ellipsoidal) {
        if (config.mandrel.k < limits::ELLIPSOIDAL_K_MIN)
            throw std::invalid_argument(
                "mandrel.k, minimum " + std::to_string(limits::ELLIPSOIDAL_K_MIN) +
                " olmali (S-GEO-03).");
    }

    // -----------------------------------------------------------------------
    // tow bolumu (opsiyonel — Phase-2 hazirlik)
    // -----------------------------------------------------------------------
    if (j.contains("tow")) {
        const auto& jt = j["tow"];

        if (jt.contains("BW"))
            config.tow.BW = jt["BW"].get<double>();
        if (jt.contains("BT"))
            config.tow.BT = jt["BT"].get<double>();
        if (jt.contains("N_tow"))
            config.tow.N_tow = jt["N_tow"].get<int>();
        if (jt.contains("Fiber_tension"))
            config.tow.Fiber_tension = jt["Fiber_tension"].get<double>();
        if (jt.contains("Winding_type")) {
            config.tow.Winding_type = jt["Winding_type"].get<std::string>();
            validateTowWindingType(config.tow.Winding_type);
        }

        // Tow validation
        if (config.tow.BW < 0.0)
            throw std::invalid_argument("tow.BW negatif olamaz.");
        if (config.tow.BT < 0.0)
            throw std::invalid_argument("tow.BT negatif olamaz.");
        if (config.tow.N_tow < 0)
            throw std::invalid_argument("tow.N_tow negatif olamaz.");
        if (config.tow.Fiber_tension < 0.0)
            throw std::invalid_argument("tow.Fiber_tension negatif olamaz.");
    }

    // -----------------------------------------------------------------------
    // winding_sequence bolumu (opsiyonel — Phase-2 hazirlik)
    // Karar-8: alpha_deg → radyan donusumu burada yapilir
    // -----------------------------------------------------------------------
    if (j.contains("winding_sequence") && j["winding_sequence"].is_array()) {
        for (const auto& jw : j["winding_sequence"]) {
            WindingLayer layer{};

            if (!jw.contains("winding_type"))
                throw std::invalid_argument(
                    "winding_sequence elemani: 'winding_type' alani eksik.");
            layer.winding_type = jw["winding_type"].get<std::string>();
            validateWindingType(layer.winding_type);

            if (!jw.contains("alpha_deg"))
                throw std::invalid_argument(
                    "winding_sequence elemani: 'alpha_deg' alani eksik.");
            double alpha_deg = jw["alpha_deg"].get<double>();

            // Karar-8: Derece → radyan donusumu (I/O katmani)
            layer.alpha_rad = alpha_deg * constants::PI / 180.0;

            if (!jw.contains("N_layers"))
                throw std::invalid_argument(
                    "winding_sequence elemani: 'N_layers' alani eksik.");
            layer.N_layers = jw["N_layers"].get<int>();

            if (layer.N_layers <= 0)
                throw std::invalid_argument(
                    "winding_sequence: N_layers pozitif olmali.");

            // Winding acisi araligi dogrulamasi (Karar-7)
            if (alpha_deg <= 0.0 || alpha_deg >= 90.0)
                throw std::invalid_argument(
                    "winding_sequence: alpha_deg (0, 90) araliginda olmali, "
                    "verilen: " + std::to_string(alpha_deg) + " derece.");

            config.winding_sequence.push_back(layer);
        }
    }

    return config;
}

// ---------------------------------------------------------------------------
// loadConfig — dosyadan okuma
// ---------------------------------------------------------------------------
WindingConfig loadConfig(const std::string& filepath)
{
    std::ifstream file(filepath);
    if (!file.is_open()) {
        throw std::runtime_error(
            "Konfigurasyon dosyasi acilamadi: " + filepath);
    }

    nlohmann::json j;
    try {
        file >> j;
    } catch (const nlohmann::json::parse_error& e) {
        throw std::runtime_error(
            "JSON parse hatasi: " + std::string(e.what()));
    }

    return parseJson(j);
}

// ---------------------------------------------------------------------------
// loadConfigFromString — string'den okuma (test amacli)
// ---------------------------------------------------------------------------
WindingConfig loadConfigFromString(const std::string& json_str)
{
    nlohmann::json j;
    try {
        j = nlohmann::json::parse(json_str);
    } catch (const nlohmann::json::parse_error& e) {
        throw std::runtime_error(
            "JSON parse hatasi: " + std::string(e.what()));
    }

    return parseJson(j);
}

} // namespace geometry
} // namespace filament
