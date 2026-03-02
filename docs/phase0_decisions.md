# Phase-0: Sistem Mimarisi Tanımı ve Formal Planlama

**Proje:** Filament Winding Yazılımı — 4-Eksen Lathe-Type COPV Üretim Makinesi  
**Faz:** Phase-0 (Mimari Kararlar)  
**Durum:** TAMAMLANDI  
**Tarih:** 2026-02-25  
**Son Güncelleme:** 2026-03-02 (GATE-1a-01 sonuçları eklendi)  
**Repo:** https://github.com/m-enesuysal/Filament-Winding

---

## 1. Genel Bakış

Bu doküman, filament sarım yazılımının Phase-1a'ya geçmeden önce alınan tüm mimari
kararları içerir. Toplam 21 karar alınmış, 5 doğrulama kapısı (GATE) tanımlanmış,
2 açık madde kayıt altına alınmıştır.

Tüm kararlar üç mühendislik perspektifinden eşzamanlı olarak değerlendirilmiştir:
kompozit üretim mühendisi, hareket kontrol mühendisi ve kritik sistem denetçisi.

---

## 2. Mimari Kararlar

### KARAR-1: 4-Eksen Tanımı ve Konvansiyonu

**Durum:** Onaylandı

Makine eksen konfigürasyonu standart lathe-winder tanımıyla birebir örtüşmektedir:

| Eksen | Hareket Tipi | Açıklama |
|-------|-------------|----------|
| **X** | Lineer | Carriage öteleme — mandrel ekseni boyunca |
| **Y** | Lineer | Cross-carriage öteleme — mandrel eksenine dik (radyal) |
| **C** | Dönel | Spindle / mandrel rotasyonu — mandrel ekseni etrafında |
| **A** | Dönel | Feed eye eğim rotasyonu |

Yazılımdaki eksen isimlendirme konvansiyonu {X, Y, C, A} olarak sabitlenmiştir.
Tüm kinematik denklemler, G-code çıktısı ve dokümantasyon bu isimlendirmeyi kullanacaktır.

**Not (Phase-3 için):** Y ekseninin işaret konvansiyonu (mandrel'e doğru = pozitif mi,
uzaklaşma = pozitif mi) Phase-3 öncesinde kesinleştirilecektir. (Bkz. AÇIK-P3-01)

---

### KARAR-2: Mandrel Geometri — Simetri Varsayımı

**Durum:** Onaylandı  
**Seçim:** Simetrik mandrel

Ön dome ve arka dome geometrik olarak aynıdır. Tek bir dome profil fonksiyonu
tüm mandreli tanımlar.

**Sonuçları:**
- Clairaut sabiti her iki dome için aynı: r₀·sin(α₀) = sabit
- İleri ve geri stroklar mandrel orta düzlemi etrafında tam simetrik
- Pattern hesabında sadece yarım strok hesaplanıp yansıtılabilir
- Kinematik çözücüde dome geçiş mantığı tek tip olur

---

### KARAR-3: Dome Profil Parametrizasyonu — Matematiksel Temsil

**Durum:** Onaylandı  
**Seçim:** Birleşik parametrik meridyen temsili

Tüm dome tipleri (isotensoid, elipsoidal, hemispherical) ortak bir arayüzle temsil
edilir: meridyen eğrisi s(t) → (ρ(t), x(t)) parametrik formunda.

- Isotensoid: numerik ODE çözümünden üretilir
- Elipsoidal: analitik parametrizasyon (a = R_eq, b = k·R_eq)
- Hemispherical: analitik parametrizasyon (yarıçap = R_eq)

Phase-2 geodesic solver dome tipinden bağımsız çalışır — tek bir ortak API.

**Koşul:** Isotensoid için numerik türev kullanılacaksa, Phase-1a'da MATLAB'da önce
analitik türev ile karşılaştırmalı doğrulama yapılacaktır. İnterpolasyon yöntemi
(kübik spline veya Chebyshev) MATLAB gate'de kanıtlanmadan Phase-1b'ye geçilmeyecektir.
(Bkz. GATE-1a)

---

### KARAR-4: Koordinat Sistemi Konvansiyonu

**Durum:** Onaylandı

**Kartezyen referans çerçevesi {E₀}:**
- **Orijin:** Mandrel geometrik merkezi (silindir orta noktası, dönme ekseni üzerinde)
- **+X:** Headstock → Tailstock yönü (standart CNC konvansiyonu)
- **+Y:** Yatay düzlemde mandrel eksenine dik (sağ el kuralı ile)
- **+Z:** Dikey yukarı

**Silindirik koordinat sistemi:**
- **ρ (rho):** Mandrel dönme ekseninden radyal mesafe (ρ ≥ 0)
- **φ (phi):** Mandrel dönme açısı (C-ekseni rotasyonu), φ = 0 referansı: feed eye tarafı
- **x:** Mandrel ekseni boyunca konum (Kartezyen X ile aynı)

**Dönüşüm formülleri:**
- X = x
- Y = ρ·cos(φ)
- Z = ρ·sin(φ)

φ = 0 referansı +Y yönüne (feed eye tarafı) karşılık gelir.

**Meridyen eğrisi koordinatları:**
- Meridyen eğrisi (ρ, x) düzleminde tanımlanır — φ'den bağımsız (dönme simetrisi)
- s: Meridyen yay uzunluğu parametresi (ekvatordan polar açıklığa doğru artan)

Bu konvansiyon Phase-1a'dan itibaren tüm hesaplarda sabittir.

---

### KARAR-5: Mandrel Geometri Giriş Parametreleri

**Durum:** Onaylandı  
**Revizyon:** Karar-12 ve Karar-21 ile güncellenmiştir.

**Tüm dome tipleri için ortak parametreler:**

| Parametre | Birim | Açıklama |
|-----------|-------|----------|
| R_eq | mm | Ekvator (silindir) yarıçapı |
| r₀ | mm | Polar açıklık yarıçapı (boss/fitting yarıçapı) |
| L_cyl | mm | Silindirik gövde uzunluğu (dome başlangıç noktaları arası) |

**Dome tipine özel ek parametreler:**

| Dome Tipi | Ek Parametre | Açıklama |
|-----------|-------------|----------|
| Isotensoid | — | Profil R_eq ve r₀'dan tek (unique) olarak belirlenir |
| Elipsoidal | k [-] | Dome aspect ratio: k = h_dome / R_eq |
| Hemispherical | — | Yarıçap = R_eq (elipsoidal k = 1 özel hali) |

**Efektif polar açıklık yarıçapı (Phase-2 için):**

RE = r₀ + BW_eff / 2 = r₀ + (N_tow × BW) / 2

Burada RE, fiber merkezi dönüş yarıçapıdır. Geodesic solver Phase-2'de r₀ değil
RE'yi kullanacaktır. BW ve N_tow parametreleri Phase-2 öncesinde tanımlanması gereken
zorunlu girdilerdir. (Bkz. Karar-12)

**Türetilen büyüklükler (hesaplanan):**
- Dome yüksekliği (h_dome)
- Toplam mandrel uzunluğu (L_total = L_cyl + 2·h_dome)
- Meridyen yay uzunluğu (s_total)
- Yüzey alanı

**Doğrulama kısıtları (input validation):**
- r₀ < R_eq
- r₀ > 0
- L_cyl ≥ 0
- Elipsoidal için: k > 0
- RE < R_eq (yani r₀ + (N_tow × BW)/2 < R_eq)

---

### KARAR-6: Numerik Altyapı — ODE Çözücü Stratejisi

**Durum:** Onaylandı  
**Seçim:** Hibrit yaklaşım

- MATLAB doğrulamasında: ode45 (Dormand-Prince) kullanılacak
- C++ implementasyonunda: Dormand-Prince (RK45) adaptif çözücü — Boost.Odeint kütüphanesi
- Polar açıklık yakınında minimum/maksimum adım boyutu sınırları Phase-1a MATLAB
  gate'de belirlenecek

**Koşul:** Phase-1a MATLAB gate'inde Boost.Odeint çıktısı ile ode45 çıktısı aynı test
senaryosunda karşılaştırılacak, maksimum sapma belgelenecektir. Bu karşılaştırma
kanıtlanmadan C++ ODE altyapısı onaylanmış sayılmaz. (Bkz. GATE-1b-01)

---

### KARAR-7: Winding Açısı Tanımı ve Konvansiyonu

**Durum:** Onaylandı  
**Seçim:** Konvansiyon A — α, mandrel ekseni ile fiber arasındaki açı

- α = 0° → fiber mandrel eksenine paralel (polar winding yönü)
- α = 90° → fiber mandrel eksenine dik (hoop winding)
- Clairaut relation: ρ·sin(α) = sabit

Tüm faz dokümantasyonunda, C++ kodunda ve MATLAB scriptlerinde α sembolü ve
"winding angle" terminolojisi kullanılacaktır.

**Winding tipi açı aralıkları:**

| Winding Tipi | Açı Aralığı |
|-------------|-------------|
| Polar | α ≈ 5°–15° |
| Helical | α ≈ 15°–75° |
| Hoop | α ≈ 85°–90° |

**Minimum winding açısı (Phase-2'de belirlenecek):**

İki kısıttan büyük olanı efektif alt sınır olarak alınacaktır:
- Geometrik alt sınır: sin(α_min_geo) = RE / R_eff
- Makinesel alt sınır: tan(α_min_mach) = πD / P (carriage/spindle hız oranından)
- Efektif alt sınır: α_min = max(α_min_geo, α_min_mach)

---

### KARAR-8: Birim Sistemi Konvansiyonu

**Durum:** Onaylandı  
**Seçim:** SI-mm sistemi

| Büyüklük | Dahili Hesap | Kullanıcı I/O | Gerekçe |
|----------|-------------|---------------|---------|
| Uzunluk | mm | mm | CNC/G-code standardı |
| Açı | radyan | derece | Trigonometrik doğruluk |
| Lineer hız | mm/s | mm/s | Trajectory ve G-code uyumu |
| Açısal hız | rad/s | deg/s | Dahili tutarlılık |
| İvme | mm/s², rad/s² | mm/s², deg/s² | Trajectory planning |
| Kuvvet/Gerilme | N, MPa | N, MPa | SI standardı |
| Zaman | s | s | SI temel birimi |

**Kritik kurallar:**
- Dahili hesaplamalarda **sadece radyan** kullanılacak
- Derece ↔ radyan dönüşümü **yalnızca** giriş/çıkış katmanında yapılacak
- Hiçbir modül içinde derece kullanılmayacak
- G-code çıkışındaki birim dönüşümü **yalnızca** Phase-5 Post-Processor'da gerçekleşecek

**Not:** Kontrolcünün G-code açı beklentisi (derece/radyan) Phase-5 başlamadan önce
teyit edilecektir. (Bkz. AÇIK-P5-01)

---

### KARAR-9: Mandrel Yüzey Temsili — Dahili Veri Yapısı

**Durum:** Onaylandı  
**Seçim:** Önceden hesaplanmış lookup table + interpolasyon

Tüm dome tipleri için meridyen profili yeterli çözünürlükte hesaplanıp tabloda saklanır.

**Tablo yapısı:** {s_i, ρ_i, x_i, dρ/ds_i, dx/ds_i, κ_i}

**Tablo çözünürlüğü:** Adaptif örnekleme ile belirlenecek:
- Polar açıklık yakınında (yüksek eğrilik): sık örnekleme
- Silindir-dome geçişinde: seyrek örnekleme
- Adaptif örnekleme kriteri Phase-1a MATLAB gate'inde kanıtlanacak

**Sorgulama:** GATE-1a'da seçilecek interpolasyon yöntemiyle (spline veya Chebyshev)

**Koşul:** İnterpolasyon hatasının kabul edilebilir düzeyde olduğu GATE-1a'da
kanıtlanmadan Phase-1b'de bu veri yapısı uygulamaya geçilemez.

---

### KARAR-10: C++ Geometry Kernel — Sınıf Hiyerarşisi

**Durum:** Onaylandı  
**Revizyon:** Karar-21 ile MandrelGeometry immutable + factory pattern olarak güncellendi.

**Mimari:**

```
IMeridianProfile  (soyut arayüz / abstract interface)
│
├── IsotensoidProfile    : IMeridianProfile
├── EllipsoidalProfile   : IMeridianProfile
└── HemisphericalProfile : IMeridianProfile

MeridianLookupTable  (interpolasyonlu tablo — tüm profil tipleri burada saklanır)

MandrelGeometry  (üst katman — silindir + iki dome birleşimi, immutable)
```

**IMeridianProfile zorunlu metodları:**
- generateProfile(params) → MeridianLookupTable
- radius(s) → ρ
- axialPosition(s) → x
- slopeAngle(s) → β
- curvature(s) → κ
- surfaceNormal(s, φ) → (n_ρ, n_φ, n_x)

**MandrelGeometry zorunlu metodları:**
- point(s_global, φ) → (ρ, φ, x)
- windingAngleToClairaut(α, ρ) → c
- isOnCylinder(s_global) → bool
- isOnDome(s_global) → bool

**Immutable + Factory Pattern (Karar-21 güncellemesi):**
Her sekans için yeni MandrelGeometry nesnesi oluşturulacak. Mutation yaklaşımı
reddedildi — çoklu sekans debug sürecinde state takibi için immutable tasarım zorunlu.

**C¹ Süreklilik Koşulu:** Silindir-dome geçişinde teğet süreklilik her dome tipi için
ekvator noktasında MATLAB'da sayısal olarak kanıtlanacak. Bu kanıt sağlanmadan
Phase-1b'de MandrelGeometry geçiş mantığı uygulamaya geçilemez. (Bkz. GATE-1a)

---

### KARAR-11: Numerik Tolerans ve Yakınsama Kriterleri

**Durum:** Onaylandı

**Katman 1 — ODE İntegrasyon Toleransı:**

| Parametre | Değer | Kullanım |
|-----------|-------|----------|
| RelTol | 1e-8 | MATLAB ode45 ve Boost.Odeint |
| AbsTol | 1e-10 | MATLAB ode45 ve Boost.Odeint |

MATLAB varsayılan toleransları (RelTol = 1e-3, AbsTol = 1e-6) kesinlikle kullanılmayacak.

**Katman 2 — İnterpolasyon Hatası Toleransı:**

| Büyüklük | Tolerans | Tip |
|----------|---------|-----|
| Pozisyon (ρ, x) | \|ε\| < 1e-4 mm | Mutlak |
| Türev (dρ/ds, dx/ds) | \|ε\| < 1e-6 | Mutlak (boyutsuz) |
| Eğrilik (κ) | \|ε_rel\| < 1e-4 | Bağıl (%0.01) |

Bu değerler GATE-1a'da ölçülen interpolasyon hatası için kabul kriteridir.
Eşik aşılırsa interpolasyon yöntemi reddedilir.

**Katman 3 — Geometrik Sorgu Toleransı:**
- Silindir-dome sınırında ±1e-6 mm tolerans bandı
- Band içindeki noktalar dome tarafına atanır
- Gerekçe: dome tarafına yanlış atama sonuç doğruya yakın, silindir tarafına
  yanlış atama eğrilik kaybına neden olur — güvenli yön dome tarafı

---

### KARAR-12: Fiber/Tow Giriş Parametreleri

**Durum:** Onaylandı  
**Revizyon:** Karar-21 ile BT zorunlu girdiye eklendi.

**Phase-2 zorunlu girdi listesi:**

| Parametre | Birim | Açıklama |
|-----------|-------|----------|
| BW | mm | Tek tow bant genişliği |
| BT | mm | Tek katman bant kalınlığı |
| N_tow | — | Aynı anda sarılan tow sayısı (ayarlanabilir) |
| Fiber_tension | N | Fiber gerginlik kuvveti |
| Winding_type | — | Wet / Dry |

**Türetilen büyüklükler:**
- Efektif bant genişliği: BW_eff = N_tow × BW
- Efektif polar açıklık yarıçapı: RE = r₀ + BW_eff / 2

---

### KARAR-13: Winding Pattern Tanımı — p/q Notasyonu

**Durum:** Onaylandı  
**Seçim:** Hibrit yaklaşım (Opsiyon C)

Her iki yön desteklenir:
- Pattern-driven: Kullanıcı p/q girer, yazılım uyumlu α hesaplar
- Angle-driven: Kullanıcı α girer, yazılım en yakın p/q önerir

L_cyl ince ayar parametresi olarak kullanılabilir — tam pattern kapama için
milimetrik düzeltme.

**Temel tanımlar:**
- Circuit (devre): Bir tam dome-to-dome-to-dome hareketi
- Pattern (desen): p adet devrenin mandrel yüzeyini tam kaplaması
- Angular step: Δφ = 2πq/p
- gcd(p,q) = 1 koşulu zorunlu (S-PAT-02)

---

### KARAR-14: Feed Eye Modeli ve Çarpışma Önleme

**Durum:** Onaylandı  
**Seçim:** Aşamalı strateji — önce bounding box, sonra quasi-numerical zarflama

**Phase-3 zorunlu girdi parametreleri:**

| Parametre | Birim | Açıklama |
|-----------|-------|----------|
| D_eye | mm | Feed eye roller çapı |
| W_eye | mm | Feed eye genişliği (mandrel eksenine paralel) |
| H_eye | mm | Feed eye yüksekliği (radyal yön) |
| L_arm | mm | Feed eye kolu uzunluğu |
| d_min | mm | Minimum izin verilen feed eye-mandrel yüzey mesafesi |

Bu değerler girilmeden Phase-3 kinematik çözücü başlatılamaz.

**Çarpışma kontrolü stratejisi:**
- Phase-3 ilk implementasyon: Opsiyon A (bounding box) — doğrulama ve debug kolaylığı
- Phase-4 öncesi geçiş: Opsiyon B (Koussios quasi-numerical zarflama) — performans
  gereksinimi ortaya çıktığında

**Kapsam sınırı:** Reçine banyosu feed eye mekanizmasının öncesinde ve dışında konumlanıyor.
Kinematik modelde ayrı bileşen olarak modellenmeyecek, çarpışma kontrolü kapsamı dışında.

---

### KARAR-15: Singülarite ve Edge Case Kataloğu

**Durum:** Onaylandı

#### Kategori 1 — Geometrik Singülariteler (Phase-1a/1b)

**S-GEO-01: Polar açıklık noktası (ρ → r₀)**
Dome meridyen profilinde polar açıklık yakınında dρ/ds → 0 ve eğrilik κ → maksimum.
Isotensoid ODE'sinde stiff davranış. Adaptif adım çözücü bu bölgede adım boyutunu
otomatik küçültecek.

**S-GEO-02: Ekvator noktası (ρ = R_eq)**
Silindir-dome geçiş noktası. Meridyen eğim açısı β = 0. C¹ süreklilik sağlanmalı.
İkinci türevde süreksizlik olabilir (C² süreksizliği).

**S-GEO-03: Elipsoidal dome k → 0 limiti**
k çok küçük olduğunda dome düz disk haline gelir — meridyen profili dejenere olur.
Pratik alt sınır: k_min = 0.15 (altında hata fırlatılır).

**S-GEO-04: Hemispherical dome k = 1 özel durumu**
Elipsoidal k → 1 limitinde EllipsoidalProfile ve HemisphericalProfile sınıfları aynı
geometriyi üretmeli. Birebir örtüşme GATE-1a'da doğrulanacak.

#### Kategori 2 — Sarım Yolu Singülariteleri (Phase-2a/2b)

**S-WIND-01: Dome dönüş noktası (ρ = RE)**
Winding açısı α → 90°. dφ/ds → ∞. Numerik integrasyonda özel işlem gerektirir.

**S-WIND-02: Hoop winding α → 90°**
Clairaut sabiti c ≈ R_eq. Fiber dome'a giremez. Hoop winding sadece silindirik bölgede
uygulanabilir — otomatik kontrol zorunlu.

**S-WIND-03: Polar winding α → 0°**
Düşük açılarda çok az açısal ilerleme. Pratik alt sınır var.

#### Kategori 3 — Kinematik Singülariteler (Phase-3)

**S-KIN-01: Determinant singülaritesi**
Lathe-winder kinematik denklemlerinde det[N] = 0 olduğunda çözüm tanımsız.
Phase-3'te tespit ve önleme gerekir.

**S-KIN-02: Dome dönüş noktasında hız değişimi**
Anlık yön değişimi. Mandrel rotasyon / carriage hız oranında hızlı değişim.
Phase-4 trajectory planning'de yumuşatılacak.

#### Kategori 4 — Pattern Singülariteleri (Phase-2b)

**S-PAT-01: Tam kaplama uyumsuzluğu**
Açısal ilerleme irrasyonel olduğunda tam p/q pattern bulunamaz.
L_cyl ince ayarı ile çözülecek.

**S-PAT-02: p/q indirgenemez kesir kontrolü**
gcd(p,q) = 1 koşulu zorunlu. gcd(p,q) > 1 olduğunda fiber p/gcd adet devrede
tam kapatmayı tamamlar ve kalan devreler üst üste biner — kaplama eksik kalır.
Bu kontrol Phase-2b'de pattern validation adımının ilk koşulu olarak uygulanacak.

**Kapsam dışı (makine kontrolcüsü sorumluluğunda):**
- Eksen stroke limitleri (X, Y, C, A)
- Acil durdurma senaryoları
- Malzeme kopma durumu

---

### KARAR-16: MATLAB ↔ C++ Doğrulama İş Akışı — Gate Protokolü

**Durum:** Onaylandı

**Adım 1 — MATLAB referans implementasyonu:**
Matematiksel model önce MATLAB'da implement edilir. Bu implementasyon referans çözüm
(ground truth) statüsündedir. MATLAB kodu repo'da `MATLAB/phaseXX/` altında
versiyon kontrollü tutulur.

**Adım 2 — MATLAB iç doğrulaması:**
Analitik çözümü bilinen test senaryolarıyla karşılaştırma, bilinen singülaritelerde
davranış kontrolü, literatürdeki sayısal örneklerle karşılaştırma. Bu adım geçilmeden
C++ implementasyonuna başlanamaz.

**Adım 3 — C++ implementasyonu:**
Doğrulanmış MATLAB modeli C++'a aktarılır.

**Adım 4 — Cross-validation:**
C++ çıktısı ile MATLAB çıktısı aynı test senaryolarında karşılaştırılır.
Karar-11 toleransları uygulanır. Sonuçlar `docs/` altında belgelenir.

**Standart test senaryoları:**

| ID | Açıklama | R_eq [mm] | r₀ [mm] |
|----|---------|-----------|---------|
| TEST-01 | ASTM Subscale | 73 | 22 |
| TEST-02 | Endüstriyel COPV | 152.4 | 45 |
| TEST-03 | Küçük Açıklık | 150 | 10 |
| TEST-04 | H₂ Aerospace | 200 | 50 |

---

### KARAR-17: Hata Yönetimi Stratejisi

**Durum:** Onaylandı  
**Seçim:** Katmanlı hibrit (Opsiyon C)

- **Üst katman** (kullanıcı girdisi, dosya I/O, konfigürasyon, mandrel oluşturma):
  Exception-based. std::out_of_range, std::invalid_argument ve özel domain exception
  sınıfları kullanılacak.
- **Alt katman** (numerik hesap çekirdeği, interpolasyon sorguları, ODE adımları,
  kinematik tight loops): Return-code based. C++17 std::optional veya std::expected
  kullanılacak. Bu katmanda hiçbir exception fırlatılmayacak.

**Katman sınırı tanımı:** MandrelGeometry ve IMeridianProfile public API'si üst katman.
MeridianLookupTable sorgu metodları ve ODE integrasyon adımları alt katman. Her metodun
hangi katmana ait olduğu header dosyalarında açıkça belgelenecek.

---

### KARAR-18: Python Visualization Katmanı

**Durum:** Onaylandı  
**Seçim:** Kademeli pybind11 binding

**Rol ayrımı:**
- **MATLAB:** Matematiksel doğrulama aracı. Sadece geliştirme sürecinde, son üründe yer almaz.
- **Python:** Son üründe kullanıcıya sunulan görselleştirme katmanı.

**Kütüphane tercihi:**
- **Matplotlib:** 2D analitik çıktılar (winding açısı dağılımı, hız/ivme profilleri, meridyen eğrisi)
- **PyVista:** 3D görselleştirme (mandrel yüzeyi, sarım yolları, process animasyonu)

**Zaman çizelgesi:**
- Phase-1a/2a: MATLAB visualization yeterli
- Phase-1b: Temel pybind11 binding başlar (tamamlanma kriteri değil)
- Her fazın tamamlanmasıyla Python katmanı kademeli genişler

---

### KARAR-19: Build Sistemi ve Proje Yapılandırması

**Durum:** Onaylandı

**Derleme hedefleri:**

| Target | Tip | Açıklama |
|--------|-----|----------|
| filament_core | Statik kütüphane | Tüm numerik hesap çekirdeği |
| filament_python | Shared library | pybind11 binding modülü |
| filament_tests | Executable | Google Test birim testleri |

**Bağımlılıklar:**
- Boost.Odeint — ODE çözücü
- Eigen3 — Lineer cebir (Phase-3 kinematik matris işlemleri)
- pybind11 — Python binding
- Google Test — Test framework
- nlohmann/json — JSON konfigürasyon parser (header-only)

**Compiler ayarları:**
- Standard: C++17
- Optimizasyon: Release -O2, Debug -O0 -g
- Uyarılar: -Wall -Wextra -Wpedantic

**KRİTİK KURAL:** `-ffast-math` hiçbir koşulda kullanılmayacak. IEEE 754 uyumluluğu
zorunludur. Gerekçe: singülarite tespiti (Karar-15) ve tolerans yönetimi (Karar-11)
bu flag ile güvenilir çalışamaz. Bu kural CMakeLists.txt yorumlarına da işlenecek.

---

### KARAR-20: Konfigürasyon Dosya Formatı

**Durum:** Onaylandı  
**Seçim:** JSON (nlohmann/json, header-only)

Konfigürasyon dosyası fazlara göre bölümlendirilecek:

```json
{
  "_units": "Uzunluk: mm, Açı: derece, Hız: mm/s, Kuvvet: N",

  "mandrel": {
    "R_eq": 152.4,
    "r0": 45.0,
    "L_cyl": 300.0,
    "dome_type": "isotensoid",
    "k": null
  },

  "tow": {
    "BW": 6.35,
    "BT": 0.25,
    "N_tow": 1,
    "Fiber_tension": 50.0,
    "Winding_type": "wet"
  },

  "winding_sequence": [
    {"winding_type": "helical", "alpha_deg": 30.0, "N_layers": 2},
    {"winding_type": "hoop", "alpha_deg": 89.0, "N_layers": 1},
    {"winding_type": "helical", "alpha_deg": 30.0, "N_layers": 2},
    {"winding_type": "hoop", "alpha_deg": 89.0, "N_layers": 1}
  ],

  "feed_eye": {
    "D_eye": 0.0,
    "W_eye": 0.0,
    "H_eye": 0.0,
    "L_arm": 0.0,
    "d_min": 0.0
  },

  "trajectory": {
  }
}
```

Birim bilgisi: Karar-8'deki SI-mm sistemi geçerli. Konfigürasyon dosyasında
uzunluklar mm, açılar derece (kullanıcı I/O katmanı), hızlar mm/s.

---

### KARAR-21: Winding Sekans Yönetimi

**Durum:** Onaylandı

Yazılım hem tek sarım sekansını hem de çoklu laminat sekansını destekleyecek.
Sekans kullanıcı tarafından konfigürasyon dosyasında tanımlanır.

**Kalınlık birikimi yaklaşımı:**

- **Kademe 1 (benimsenen):** Silindirik bölgede R_eff = R_eq + n × BT.
  Dome profili değişmez kabul edilir.
- **Kademe 2 (Phase-2a sonrası iyileştirme hedefi):** Dome kalınlık dağılımı
  t(s) = BT / cos(α(s)) ile güncellenir. Phase-2a geodesic solver çıktısı gerektirir.
- **Kademe 3 (kapsam dışı):** Tam iteratif çözüm.

**Cascading etkiler:**
Her katman sonrası R_eff güncellendiğinde:
- sin(α_min_geo) = RE / R_eff yeniden hesaplanır
- R_eff arttıkça izin verilen minimum winding açısı azalır
- Sekans sıralaması fiziksel sonucu etkiler — Phase-2b'de dikkate alınacak

**MandrelGeometry tasarımı:** Immutable + factory pattern. Her sekans için yeni
MandrelGeometry nesnesi oluşturulur. Bellek maliyeti kabul edilebilir (katman sayısı
onlarla sınırlı).

---

## 3. Doğrulama Kapıları (GATE)

### GATE-1a: Phase-1a Geometri Doğrulama Kapısı

**Durum:** ✅ ONAYLI (2026-03-02)

Phase-1a'dan Phase-1b'ye geçiş için tüm koşulların sağlanması zorunludur.

| # | Koşul | Kriter | Sonuç |
|---|-------|--------|-------|
| 1 | İnterpolasyon yöntemi seçimi | Spline vs Chebyshev — analitik türev referansıyla hata analizi | ✅ |
| 2 | Adaptif örnekleme kriteri | Eğrilik tabanlı yoğunlaştırma kuralı belirlenmesi | ✅ |
| 3 | Tablo çözünürlüğü vs hata | Trade-off analizi, Karar-11 toleransları referans | ✅ |
| 4 | C¹ süreklilik doğrulaması | Her dome tipi için ekvator noktasında teğet sürekliliği | ✅ |
| 5 | k=1 örtüşme doğrulaması | EllipsoidalProfile(k=1) = HemisphericalProfile (S-GEO-04) | ✅ |

**GATE-1a-01 Sonuçları:**

Hemispherical, elipsoidal ve isotensoid dome profilleri MATLAB'da doğrulanmıştır.
Toplam **28/28 test PASS** — tüm profil tipleri, tüm test senaryoları (TEST-01 ila
TEST-04), C¹ süreklilik, S-GEO-04 k=1 örtüşme, eğrilik çapraz doğrulama ve uç
durum senaryoları dahil.

**Silindirik bölge kararı:** Silindirik bölge için ayrı MATLAB doğrulama scripti
gereksizdir. ρ = R_eq = sabit, κ_m = 0, β = 0 olan trivial geometri, Phase-1b'de
MandrelGeometry sınıfı içinde analitik olarak tanımlanacaktır. Ayrı bir profil
fonksiyonu veya lookup table gerekmez.

**Struct arayüz yaması:** GATE-1a doğrulama sürecinde profil struct çıktı
arayüzünde eksik alanlar tespit edilmiş ve tüm profil fonksiyonlarına
(hemispherical_dome_profile, ellipsoidal_dome_profile, isotensoid_dome_profile)
aşağıdaki alanlar eklenmiştir:

| Alan | Açıklama |
|------|----------|
| `kappa_eq` | Ekvator meridyen eğriliği κ_m(s=0) [1/mm] |
| `kappa_pol` | Polar açıklık meridyen eğriliği κ_m(s=s_total) [1/mm] |
| `alpha_w` | Ekvator winding açısı α_w = arcsin(r₀/R_eq) [rad] |
| `aspect_r` | Dome aspect ratio (hemispherical: 1.0, elipsoidal: k, isotensoid: h_dome/R_eq) |

Bu alanlar Phase-1b C++ arayüzüne ve Phase-2 geodesic solver girdilerine
doğrudan karşılık gelmektedir. Tüm profil fonksiyonları artık ortak bir üst
küme struct arayüzü üzerinden tutarlı çıktı üretmektedir.

Tolerans kriterleri: Karar-11 Katman 2 değerleri uygulanmıştır.
Test senaryoları: Karar-16 TEST-01 ila TEST-04.

### GATE-1b-01: C++ ODE Altyapısı Doğrulama Kapısı

**Durum:** Bekliyor

| # | Koşul | Kriter |
|---|-------|--------|
| 1 | Aynı test senaryosu MATLAB ode45 ve Boost.Odeint ile çözülecek | TEST-01 ila TEST-04 |
| 2 | Meridyen profil noktaları (ρ(s), x(s)) karşılaştırılacak | Karar-11 Katman 2 toleransları |
| 3 | Birinci türevler ve eğrilik karşılaştırılacak | Karar-11 Katman 2 toleransları |
| 4 | Polar açıklık yakını sapma ayrıca raporlanacak | r → r₀ bölgesi özel inceleme |
| 5 | Maksimum sapma ve RMS sapma belgelenecek | docs/ altında rapor |

### GATE-2a: Phase-2a Geodesic Path Doğrulama Kapısı

**Durum:** Bekliyor  
Phase-2a'da tanımlanacak. Yapı: GATE-1a ile aynı protokol.

### GATE-2b: Phase-2b Pattern Doğrulama Kapısı

**Durum:** Bekliyor  
Phase-2b'da tanımlanacak. İlk koşul: gcd(p,q) = 1 kontrolü (S-PAT-02).

---

## 4. Açık Maddeler

| ID | Açıklama | Çözülecek Faz |
|----|---------|---------------|
| AÇIK-P3-01 | Y-ekseni işaret konvansiyonu (mandrel'e doğru +/-) | Phase-3 öncesi |
| AÇIK-P5-01 | G-code açı birimi (makine kontrolcüsü: derece/radyan) | Phase-5 öncesi |

---

## 5. Phase-1a Sonuçları Özeti

### Doğrulanan Dome Profilleri

| Dome Tipi | Matematik | MATLAB | Doğrulama | Durum |
|-----------|-----------|--------|-----------|-------|
| Hemispherical | ✅ | ✅ | ✅ 28/28 PASS | Tamamlandı |
| Elipsoidal | ✅ | ✅ | ✅ 28/28 PASS | Tamamlandı |
| Isotensoid | ✅ | ✅ | ✅ 28/28 PASS | Tamamlandı |
| Silindirik bölge | — | — | Gereksiz | Phase-1b'de analitik |

### MATLAB Dosyaları (MATLAB/phase1a_geometry/)

| Dosya | İşlev |
|-------|-------|
| `hemispherical_dome_profile.m` | Hemispherical dome meridyen profili üretici |
| `ellipsoidal_dome_profile.m` | Elipsoidal dome meridyen profili üretici |
| `isotensoid_dome_profile.m` | Isotensoid dome meridyen profili üretici (ODE tabanlı) |
| `verify_hemispherical_dome.m` | Hemispherical doğrulama scripti |
| `verify_ellipsoidal_dome.m` | Elipsoidal doğrulama scripti (S-GEO-04 dahil) |
| `verify_isotensoid_dome.m` | Isotensoid doğrulama scripti |

### Kritik Mühendislik Notları (Phase-1a'dan)

- Ekvator eğrilik sıçraması dome tipine bağlıdır: hemispherical κ_m(0) = 1/R_eq,
  elipsoidal κ_m(0) = 1/(R_eq·k²), isotensoid κ_m(0) ODE çözümünden.
  Yassı dome'larda (k < 1) sıçrama daha büyüktür.
- S-GEO-04 doğrulanmıştır: ellipsoidal(k=1) ve hemispherical profilleri
  Karar-11 toleransları içinde birebir örtüşmektedir.
- Tüm profil fonksiyonları ortak struct arayüzü üzerinden tutarlı çıktı üretmektedir
  (kappa_eq, kappa_pol, alpha_w, aspect_r alanları dahil).

---

## 6. Referanslar

- Koussios, S. — "Filament Winding: A Unified Approach" (birincil referans)
- Do Carmo — "Differential Geometry of Curves and Surfaces"
- Numerical Recipes, 3rd Edition
- ASTM D 2585 (COPV test standardı)
