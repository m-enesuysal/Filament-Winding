Phase-1b Bağımlılık Grafiği
S1: Build Altyapısı + Dizin İskeleti
│
├──► S2: MeridianLookupTable (veri yapısı + interpolasyon)
│    │
│    ├──► S3: IMeridianProfile arayüzü + HemisphericalProfile
│    │    │
│    │    ├──► S4: EllipsoidalProfile (+ S-GEO-04 k=1 testi)
│    │    │
│    │    └──► S5: IsotensoidProfile (Boost.Odeint ODE çözücü)
│    │         │
│    │         └──► S6: MandrelGeometry (silindir + dome birleşimi)
│    │              │
│    │              └──► S7: GATE-1b-01 Cross-Validation
│    │
│    └──► S8: JSON Konfigürasyon Parser
Kritik yol: S1 → S2 → S3 → S5 → S6 → S7. S4 ve S8 paralel çalışabilir ama session bazında sıralı gideceğiz.
________________________________________
Session Planı
S1 — Build Altyapısı ve Dizin İskeleti
Hedef: Boş ama derlenebilir proje yapısı.
Kapsam:
•	Repo dizin yapısının oluşturulması (system prompt'taki yapıya uygun: src/geometry/, include/geometry/, tests/, docs/)
•	Kök CMakeLists.txt + alt dizin CMakeLists.txt dosyaları
•	Karar-19 bağımlılıklarının tanımlanması: Boost.Odeint (header-only), Eigen3, Google Test, nlohmann/json (header-only)
•	-ffast-math yasağının CMake yorumlarına işlenmesi (Karar-19 kritik kuralı)
•	Compiler flags: C++17, -Wall -Wextra -Wpedantic
•	Üç target'ın iskeleti: filament_core (statik lib), filament_tests (executable), filament_python (placeholder)
•	Tek bir "hello world" düzeyinde derleme testi — cmake --build başarılı olmalı
Tamamlanma kriteri: cmake -B build && cmake --build build hatasız. Google Test "boş test" çalışıyor.
Çıktılar: CMakeLists.txt (kök + alt dizinler), boş placeholder header/source dosyaları, README güncelleme
________________________________________
S2 — MeridianLookupTable
Hedef: Tüm dome profillerinin ortak veri taşıyıcısı ve interpolasyon altyapısı.
Kapsam:
•	include/geometry/meridian_lookup_table.h + src/geometry/meridian_lookup_table.cpp
•	Tablo struct tanımı: {s_i, ρ_i, x_i, dρ/ds_i, dx/ds_i, κ_m_i} (Karar-9)
•	Skaler meta-veri alanları: R_eq, r₀, s_total, h_dome, A_dome, kappa_eq, kappa_pol, alpha_w, aspect_r (struct arayüz yaması uyumlu)
•	Interpolasyon motoru: kübik spline (GATE-1a'da seçilen yöntem) — s parametresinden tüm büyüklükleri sorgulama
•	Binary search ile O(log n) segment tespiti
•	Hata yönetimi: alt katman kuralları (Karar-17) — return-code / std::optional, exception yok
•	Birim testleri: bilinen analitik veriden (hemispherical R=100, r₀=30 gibi) tablo oluştur, interpolasyon sorgula, Karar-11 Katman 2 toleranslarını kontrol et
Tamamlanma kriteri: Tüm birim testleri PASS. Bilinen analitik referansa karşı interpolasyon hataları Karar-11 toleransları içinde.
Bağımlılık: S1
________________________________________
S3 — IMeridianProfile Arayüzü + HemisphericalProfile
Hedef: Soyut arayüz tanımı ve en basit implementasyon.
Kapsam:
•	include/geometry/i_meridian_profile.h — soyut sınıf (Karar-10 zorunlu metodları)
•	include/geometry/hemispherical_profile.h + src/geometry/hemispherical_profile.cpp
•	generateProfile(R_eq, r₀, N_points) → MeridianLookupTable implementasyonu
•	Kapalı-form hesap: ρ(θ), x(θ), dρ/ds, dx/ds, κ_m = 1/R_eq (Phase-1a MATLAB matematiğinin C++ karşılığı)
•	Input validation: r₀ < R_eq, r₀ > 0 (Karar-5 kısıtları) — üst katman exception (Karar-17)
•	Birim testleri: MATLAB hemispherical_dome_profile.m ile aynı parametreler → aynı çıktılar (TEST-01 ila TEST-04)
•	Testlerde tolerans referansı: Karar-11 Katman 2
Tamamlanma kriteri: TEST-01..04 parametreleriyle C++ çıktısı MATLAB referans değerlerine Karar-11 toleransları içinde eşit. Birim testler PASS.
Bağımlılık: S2
Not: Bu session'da MATLAB referans değerleri .csv veya hardcoded expected values olarak test dosyasına gömülecek. Tam GATE-1b-01 cross-validation S7'de yapılacak.
________________________________________
S4 — EllipsoidalProfile (+ S-GEO-04)
Hedef: Elipsoidal dome implementasyonu ve k=1 örtüşme doğrulaması.
Kapsam:
•	include/geometry/ellipsoidal_profile.h + src/geometry/ellipsoidal_profile.cpp
•	Parametrik hesap: ρ(θ) = R_eq·cos(θ), x(θ) = k·R_eq·sin(θ), f(θ) = √(sin²θ + k²·cos²θ)
•	Yay uzunluğu: sayısal integrasyon (eliptik integral — Simpson veya adaptif quadrature)
•	Eğrilik: analitik κ_m(θ) = k / (R_eq · f³(θ))
•	Input validation: k > 0, k_min = 0.15 alt sınır (S-GEO-03), r₀ < R_eq
•	S-GEO-04 testi: EllipsoidalProfile(k=1.0) çıktısı HemisphericalProfile çıktısıyla karşılaştırma — Karar-11 toleransları içinde birebir örtüşme
•	Birim testleri: TEST-01..04 (farklı k değerleriyle), S-GEO-03 uç durum (k=0.15 sınırda, k=0.1 hata fırlatma)
Tamamlanma kriteri: Tüm birim testler PASS. S-GEO-04 örtüşme doğrulanmış.
Bağımlılık: S3 (HemisphericalProfile referans olarak gerekli)
________________________________________
S5 — IsotensoidProfile (ODE Çözücü)
Hedef: En karmaşık dome profili — Boost.Odeint ile numerik ODE çözümü.
Kapsam:
•	include/geometry/isotensoid_profile.h + src/geometry/isotensoid_profile.cpp
•	Netting theory ODE sistemi: Phase-1a'da MATLAB'da doğrulanmış diferansiyel denklemlerin C++ karşılığı
•	Boost.Odeint Dormand-Prince (RK45) adaptif çözücü entegrasyonu
•	Toleranslar: RelTol = 1e-8, AbsTol = 1e-10 (Karar-11 Katman 1)
•	S-GEO-01 yönetimi: polar açıklık yakınında adaptif adım davranışının kontrolü
•	Profil sonrası uniform s ızgarasına yeniden örnekleme (MeridianLookupTable'a yazım)
•	Eğrilik hesabı: ODE çözümünden analitik türetme (Phase-1a'da kanıtlanmış formüller)
•	Birim testleri: TEST-01..04, MATLAB referans değerleriyle karşılaştırma
Tamamlanma kriteri: Tüm birim testler PASS. Boost.Odeint çıktısı MATLAB ode45 referansıyla uyumlu.
Bağımlılık: S2 (MeridianLookupTable)
Risk: Boost.Odeint API'si MATLAB ode45'ten farklı — adaptif adım kontrol parametrelerinin eşlenmesi dikkat gerektirir. Session içinde bu eşleme belgelenmeli.
________________________________________
S6 — MandrelGeometry (Silindir + Dome Birleşimi)
Hedef: Tam mandrel yüzey temsili — Karar-10 üst katman sınıfı.
Kapsam:
•	include/geometry/mandrel_geometry.h + src/geometry/mandrel_geometry.cpp
•	Constructor: dome_type enum + parametreler → uygun IMeridianProfile oluşturma (factory pattern)
•	Silindirik bölge: analitik tanım (ρ = R_eq sabit, κ_m = 0, β = 0) — Phase-0 kararı uyarınca ayrı lookup table gereksiz
•	Global yay uzunluğu koordinat sistemi: s_global ∈ [0, s_dome + L_cyl + s_dome]
•	Bölge tespiti: isOnDome1 / isOnCylinder / isOnDome2 — Karar-11 Katman 3 tolerans bandı (±1e-6 mm, dome tarafına atama)
•	Sorgulama: point(s_global, φ), windingAngleToClairaut(α, ρ), curvature(s_global)
•	C¹ süreklilik: silindir-dome geçişinde teğet sürekliliğinin runtime doğrulaması (assertion)
•	Immutable tasarım (Karar-21): constructor'da tüm hesap yapılır, sonra sadece sorgulama
•	Birim testleri: geçiş noktası sürekliliği, global-lokal koordinat dönüşümü tutarlılığı, tüm dome tipleriyle entegrasyon testi
Tamamlanma kriteri: Tüm dome tipleri MandrelGeometry üzerinden sorgulanabiliyor. Geçiş noktası C¹ sürekliliği birim testlerle doğrulanmış. Birim testler PASS.
Bağımlılık: S3, S4, S5 (tüm profil implementasyonları)
________________________________________
S7 — GATE-1b-01 Cross-Validation
Hedef: Phase-1b tamamlanma kapısı — C++ vs MATLAB karşılaştırması.
Kapsam:
•	MATLAB referans verilerinin dışa aktarılması: her profil tipi × TEST-01..04 = 12 veri seti → .csv dosyaları (s, ρ, x, dρ/ds, dx/ds, κ_m)
•	C++ tarafında aynı parametrelerle profil üretimi → nokta bazında karşılaştırma
•	GATE-1b-01 koşullarının tamamı (phase0_decisions.md Bölüm 3): 
o	Meridyen profil noktaları (ρ(s), x(s)) karşılaştırma
o	Birinci türevler ve eğrilik karşılaştırma
o	Polar açıklık yakını sapma özel raporu
o	Maksimum sapma ve RMS sapma hesabı
•	Sonuç raporu: docs/gate_1b_01_report.md — Türkçe
•	MandrelGeometry entegrasyon testi: tam mandrel üzerinden global sorgulama doğrulaması
Tamamlanma kriteri: GATE-1b-01 tüm koşulları PASS. Rapor üretilmiş. Phase-1b kapatılabilir durumda.
Bağımlılık: S6, MATLAB .csv referans verileri
________________________________________
S8 — JSON Konfigürasyon Parser (paralel çalışabilir)
Hedef: Karar-20 konfigürasyon dosya formatının C++ tarafında okunması.
Kapsam:
•	include/geometry/config_parser.h + src/geometry/config_parser.cpp
•	nlohmann/json ile JSON okuma
•	mandrel bölümü: R_eq, r₀, L_cyl, dome_type, k → MandrelGeometry oluşturma
•	tow bölümü: BW, BT, N_tow, Fiber_tension, Winding_type → struct (Phase-2 için hazırlık)
•	winding_sequence bölümü: parse + validation (Phase-2 için hazırlık)
•	Derece → radyan dönüşümü burada (Karar-8 I/O katmanı kuralı)
•	Input validation: Karar-5 kısıtları, dome_type enum kontrolü
•	Birim testleri: geçerli JSON → doğru parse, geçersiz JSON → anlamlı hata mesajı
Tamamlanma kriteri: Örnek konfigürasyon dosyası (Karar-20'deki şablon) başarıyla parse edilip MandrelGeometry oluşturulabiliyor. Birim testler PASS.
Bağımlılık: S1 (build), S6 (MandrelGeometry — tam entegrasyon testi için)
Not: S8, S2-S5 ile paralel başlatılabilir; MandrelGeometry'ye bağımlı olan kısım sadece entegrasyon testi. Parser'ın kendisi bağımsız yazılabilir.
________________________________________
Özet Tablo
Session	Hedef	Bağımlılık	Tahmini Dosya Sayısı	Kritik Karar Referansları
S1	Build + dizin iskeleti	—	~8 (CMake + placeholder)	K-19
S2	MeridianLookupTable	S1	2 (.h/.cpp) + test	K-9, K-11
S3	IMeridianProfile + Hemispherical	S2	3 (.h×2, .cpp) + test	K-10, K-17
S4	EllipsoidalProfile + S-GEO-04	S3	2 (.h/.cpp) + test	K-10, S-GEO-03/04
S5	IsotensoidProfile (ODE)	S2	2 (.h/.cpp) + test	K-6, K-11, S-GEO-01
S6	MandrelGeometry	S3, S4, S5	2 (.h/.cpp) + test	K-10, K-21, K-11
S7	GATE-1b-01 cross-validation	S6	rapor + test	GATE-1b-01
S8	JSON config parser	S1 (S6 opsiyonel)	2 (.h/.cpp) + test	K-20, K-8
Kritik yol: S1 → S2 → S3 → S5 → S6 → S7 (6 ardışık session)
Paralel çalışabilir çift: S4 ∥ S5 (ikisi de sadece S2/S3'e bağlı, birbirine değil). S8 erken başlatılabilir.
________________________________________
Session Sıralaması Önerisi
S1 → S2 → S3 → S4 → S5 → S6 → S7 → S8
Gerekçe: S4 ve S5 teorik olarak paralel ama pratikte sıralı gideceğiz. S4'ü S5'ten önce yapmak mantıklı çünkü elipsoidal daha basit — S3'teki pattern'ı (kapalı form profil → tablo) pekiştirir, S5'teki ODE karmaşıklığına geçişte adım oluşturur. S8 en sona bırakılabilir çünkü Phase-2'ye kadar zorunlu değil ama Phase-1b kapsamında tamamlanması repo bütünlüğü açısından tercih edilir.
Her session sonunda .tar paket çıktısı verilecek, session başlangıcında önceki session'ın dosyaları repo'da mevcut olmalı.

