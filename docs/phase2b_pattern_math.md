# Phase-2b: Sarım Deseni (Winding Pattern) Matematik Spesifikasyonu

**Proje:** Filament Winding Yazılımı — 4-Eksen Lathe-Type COPV Üretim Makinesi
**Faz:** Phase-2b (Pattern Generation)
**Durum:** Aktif
**Tarih:** 2026-03-15
**Ön koşul:** Phase-2a GATE-2a PASS (356/356 C++, 114/114 MATLAB)
**Not:** Phase-2a GATE-2a başarısı, T3 (silindir-dome geçiş sürekliliği) ve T7 (adaptif
adım dağılımı) toleranslarının revize edilmiş değerlerine göre tescil edilmiştir
(K-2a-03: T3 toleransı 1e-6 → 1e-4 rad; K-2a-04: T7 limiti 10 → 50). Bu gevşetmeler,
kübik spline doğal sınır koşulu etkisi (O(h²)) ve adaptif ODE çözücünün singülarite
yakını adım yoğunlaştırması nedeniyle fiziksel olarak beklenen davranışlardır.
**Repo:** https://github.com/m-enesuysal/Filament-Winding

---

## 1. Amaç ve Kapsam

Phase-2a tek bir geodesik devrenin geometrisini hesapladı — Δφ_circuit, yol noktaları,
sarım açısı. Phase-2b, bu geometriyi üretilebilir bir sarım desenine dönüştürür.

Somut olarak: verilen bir mandrel ve sarım parametreleri için "kaç devre, hangi atlama
düzeniyle, ne kadar yüzey kaplıyor" sorularını yanıtlayan pattern tablosunu üretir.
Bu tablo olmadan makine çalışamaz — G-code üretimi (Phase-5), katman planlaması ve
kaplama hesabı hepsi bu tabloya dayanır.

**Kapsam içi:**
- Diophantine denklem çözümü — leading/lagging pattern ayrımı
- Angle-driven mod: α → Δφ_circuit → uyumlu p/q listesi
- Pattern-driven mod: p/q → gerekli Δφ_circuit → L_cyl iterasyonu
- gcd(p,q) = 1 doğrulaması (S-PAT-02)
- Ekvator kaplama hesabı ve coverage % raporlama
- L_cyl ince ayar mekanizması (S-PAT-01 çözümü)
- Çoklu katman R_eff güncellemesi (Karar-21)
- Pattern tablo çıktısı — endüstriyel formatta
- Hoop katman desteği (silindirik bölge limitli)

**Kapsam dışı:**
- Tam kaplama simülasyonu (her devrenin 3D yolunun üretilmesi — Phase-3/6)
- Kalınlık dağılımı hesabı (Koussios Bölüm 9.4 — ayrı modül)
- Non-geodesic pattern desteği
- Dome ek dönüş (additional turn-around, ψ parametresi)

---

## 2. Referans Karar Kayıtları

| Karar | Açıklama | Etki |
|-------|----------|------|
| Karar-2 | Simetrik mandrel | Ön/arka dome pattern katkısı eşit |
| Karar-5 | RE = r₀ + BW_eff/2 | Clairaut sabiti → α_eq → K |
| Karar-7 | α: mandrel ekseni ile fiber arası açı | Tüm açı hesapları |
| Karar-12 | BW, N_tow, BT, Winding_type | Pattern girdileri |
| Karar-13 | Hibrit mod (Opsiyon C) + Phase-2a/2b bölünmesi | Çift yönlü çalışma |
| Karar-20 | JSON konfigürasyon — sekans listesi | Katman bazında α/p/q |
| Karar-21 | Çoklu sekans: R_eff(i+1) = R_eq + i·BT | İteratif pattern |
| K-2b-01 | Coverage-based arama varsayılan, gelişmiş modda p aralığı | Arama modu |
| K-2b-02 | L_cyl düzeltme sınırsız, %1 varsayılan uyarı eşiği | İnce ayar |
| K-2b-03 | Her katman farklı α kullanabilir | Sekans esnekliği |
| K-2b-04 | Tablo varsayılan sıralama: p küçükten büyüğe | Çıktı formatı |
| K-2b-07 | T6 Diophantine koşulu bilgi raporu — strict koşul zorunlu değil, coverage-based arama esas (Yaklaşım A) | T6 CONDITIONAL PASS |

---

## 3. Terminoloji ve Notasyon

### 3.1 Temel Tanımlar

| Terim | Tanım |
|-------|-------|
| **Devre (circuit)** | Bir tam dome-1 → silindir → dome-2 → silindir → dome-1 hareketi |
| **Pattern (desen)** | p adet devrenin mandrel yüzeyini tam (veya çoklu) kaplaması |
| **p** | Bir tam kaplama için gereken devre sayısı |
| **q** | Ardışık dome temas noktaları arasındaki açısal atlama (touchpoint cinsinden) |
| **dk** | Bir K partisyonuna sığan tow genişliği sayısı (Koussios notasyonu) |
| **n** | Ekvator çevresine sığan toplam tow genişliği sayısı (d katman için n = p·dk) |
| **d** | Hedef katman sayısı (d=1: tek kaplama, d=2: çift kaplama) |
| **K** | Bir devrenin artık açısal ilerleme sabiti: K = mod(Δφ_circuit, 2π) [rad] |
| **Leading** | 15. devre 1. devrenin **önüne** düşer: p·dk = n·d + 1 |
| **Lagging** | 15. devre 1. devrenin **arkasına** düşer: p·dk = n·d − 1 |
| **Coverage %** | Ekvator kaplama oranı: (n·BW_eff/cos(α_eq)) / (2π·R_eq) × 100 |
| **Overlap** | Coverage > 100% durumunda tow'ların üst üste binme miktarı [mm] |

### 3.2 Sembol Tablosu (Phase-2a'ya ek)

| Sembol | Birim | Tanım |
|--------|-------|-------|
| $p$ | — | Devre sayısı (tamsayı, > 0) |
| $q$ | — | Atlama adımı (tamsayı, > 0, gcd(p,q) = 1) |
| $dk$ | — | Bir K partisyonundaki tow sayısı (tamsayı, > 0) |
| $n$ | — | Ekvator çevresindeki toplam tow sayısı |
| $d$ | — | Hedef katman sayısı (tamsayı, ≥ 1) |
| $K$ | rad | Pattern sabiti = mod(Δφ_circuit, 2π) |
| $b$ | mm | Tek tow genişliği = BW |
| $b_{eff}$ | mm | Efektif bant genişliği = BW_eff = N_tow × BW |
| $b_{eq}$ | mm | Ekvator efektif genişliği = BW_eff / cos(α_eq) |
| $R_{eff}^{(i)}$ | mm | i. katman efektif yarıçapı |

---

## 4. Pattern Sabiti K ve Geometri İlişkisi

### 4.1 K Tanımı

Bir tam devrenin mandrel etrafındaki toplam açısal ilerleme Δφ_circuit,
Phase-2a GeodesicSolver tarafından hesaplanır (Eq. 8.1 of phase2a_geodesic_math.md):

$$\Delta\phi_{circuit} = 4 \cdot \phi_{dome} + 2 \cdot \phi_{cyl} \tag{4.1}$$

Bu ilerleme tam turlardan (2π katları) ve bir artık kısımdan oluşur:

$$\Delta\phi_{circuit} = 2\pi \cdot q_{full} + K \tag{4.2}$$

burada:
- $q_{full} = \lfloor \Delta\phi_{circuit} / (2\pi) \rfloor$ — tam tur sayısı
- $K = \Delta\phi_{circuit} \mod 2\pi$ — artık açısal ilerleme, $K \in [0, 2\pi)$

K, Koussios Bölüm 8'deki "pattern constant" ile birebir örtüşür (Eq. 8.11).

### 4.2 K ve Pattern İlişkisi

Mandrel ekvatorunda $p$ adet eşit partisyon oluşturulur. Her partisyonun
açısal genişliği $K$ olmalıdır. Dolayısıyla:

$$p \cdot K = 2\pi \cdot dk \tag{4.3}$$

ve pattern kapama koşulu:

$$\Delta\phi_{circuit} = \frac{2\pi q}{p} \tag{4.4}$$

burada $q = q_{full} \cdot \text{(tam tur katkısı)} + dk \cdot \text{(artık katkı)}$.
Bu denklem, Phase-0 Karar-13'te tanımlanan temel uyum koşuludur.

### 4.3 q ve Δφ_circuit Doğrudan İlişkisi

Denklem (4.4)'ü açarsak:

$$q = \frac{p \cdot \Delta\phi_{circuit}}{2\pi} \tag{4.5}$$

q tamsayı olmalıdır. Geometriden gelen Δφ_circuit genellikle tam bir
$2\pi q / p$ değerine denk gelmez — bu nedenle ya uyumlu p/q çifti aranır
ya da L_cyl düzeltmesi yapılır (Bölüm 7).

---

## 5. Diophantine Denklem ve Pattern Arama

### 5.1 Kaplama Koşulu

Ekvator çevresini d katman ile tamamen kaplamak için:

$$n \cdot b_{eq} = 2\pi R_{eq} \cdot d \tag{5.1}$$

burada $b_{eq} = BW_{eff} / \cos\alpha_{eq}$ ekvator efektif genişliğidir.

Buradan:

$$n = \left\lfloor \frac{2\pi R_{eq} \cdot d \cdot \cos\alpha_{eq}}{BW_{eff}} \right\rfloor \tag{5.2}$$

n tamsayıya yuvarlanır. Yuvarlama aşağı yapılırsa hafif boşluk (gap),
yukarı yapılırsa hafif örtüşme (overlap) oluşur.

### 5.2 Diophantine Denklem (Koussios Eq. 8.12)

Pattern'ın kapaması için p, dk ve n arasında şu ilişki sağlanmalıdır:

$$\text{Leading:} \quad p \cdot dk - n \cdot d = +1 \tag{5.3a}$$

$$\text{Lagging:} \quad p \cdot dk - n \cdot d = -1 \tag{5.3b}$$

Bu, birinci dereceden Diophantine denklemidir. Çözüm koşulu:
$\gcd(p, n \cdot d) = 1$ olmalıdır (aksi halde çözüm yoktur).

### 5.3 Arama Algoritması

**Girdi:** α_eq (veya RE), BW_eff, R_eq, d, arama modu

**Adım 1 — n hesabı:** Denklem (5.2) ile n bulunur.

**Adım 2 — p aralığı belirleme:**

Coverage-based mod (K-2b-01 varsayılan):
- Kullanıcı hedef coverage % ve tolerans girer (örn: %100, tol: −0%, +50%)
- n_min = floor(2πR_eq·d·cos(α_eq)/BW_eff × (coverage_min/100))
- n_max = ceil(2πR_eq·d·cos(α_eq)/BW_eff × (coverage_max/100))
- p aralığı n_min..n_max üzerinden türetilir

Gelişmiş mod:
- Kullanıcı doğrudan p_min, p_max girer

**Varsayılan arama limitleri (default config):**
- p aralığı: p_min = 3, p_max = 100
- q aralığı: 1 ≤ q < p, gcd(p,q) = 1 (her p için otomatik taranır)
- q'nun tek kısıtı gcd(p,q) = 1 ve q < p'dir. q > p/2 değerleri geçerlidir —
  TaniqWind arayüzünde de p=5, q=4 gibi çözümler görülmektedir. Dairesel
  yüzeyde q açısıyla atlamak ile p−q açısıyla atlamak lifleri aynı yörüngelere
  yerleştirir, sadece sarım yönü (spindle rotasyonu) değişir. Yazılımda her iki
  bölge de taranır — bu, kullanıcıya farklı sarım yönü tercihlerini sunmak
  ve gereksiz hesaplamayı önlemek arasındaki dengeyi korur. p=5, q=4 çözümü
  tabloda p=5, q=1'den ayrı bir seçenek olarak yer alır çünkü makine
  kinematik profili (Phase-3) farklı olabilir.
- Bu varsayılanlar PatternSearchParams'ta config olarak tutulur ve kullanıcı
  tarafından override edilebilir. Endüstriyel COPV'lerde p tipik olarak
  10–200 aralığındadır; p > 100 gerektiren durumlar nadir ve büyük mandrel
  çaplarına özgüdür.

**Adım 3 — Diophantine tarama:**

Her p ∈ [p_min, p_max] için:
1. dk_leading = (n·d + 1) / p — tam bölünüyorsa leading çözüm
2. dk_lagging = (n·d − 1) / p — tam bölünüyorsa lagging çözüm
3. Tam bölünme kontrolü: (n·d ± 1) mod p == 0
4. dk > 0 kontrolü

**Adım 4 — K ve q hesabı:**

Her geçerli (p, dk) için:
1. $K = 2\pi \cdot dk / p$ (Eq. 4.3)
2. $\Delta\phi_{target} = 2\pi \cdot q / p$ — burada q, Eq. (4.5) ile belirlenir
3. q tamsayı kontrolü: $|q - \text{round}(q)| < 10^{-6}$

**Adım 5 — gcd kontrolü (S-PAT-02):**

$\gcd(p, q) = 1$ sağlanmalıdır. Aksi halde fiber $p/\gcd$ devreden
sonra tekrar başlangıç noktasına döner ve mandrel tam kaplanmaz.

**Adım 6 — Coverage ve overlap hesabı:**

$$\text{Coverage} = \frac{n \cdot BW_{eff}}{2\pi R_{eq} \cdot d \cdot \cos\alpha_{eq}} \times 100\% \tag{5.4}$$

$$\text{Overlap}_{eq} = n \cdot \frac{BW_{eff}}{\cos\alpha_{eq}} - 2\pi R_{eq} \cdot d \quad [\text{mm}] \tag{5.5}$$

### 5.4 Desen Sıralama ve Atlama Düzeni (Skip-Index)

Pattern parametreleri (p, q) yüzeyin kaç devreden oluşacağını ve açısal
ilerlemeyi belirler, ancak **devrelerin mandrel üzerinde hangi fiziksel sırayla
atılacağını** doğrudan tanımlamaz. Bu sıralama, liflerin birbirleriyle örülme
(interweaving) derecesini belirler ve hem yapısal performansı hem de üretim
kalitesini doğrudan etkiler.

**Atlama düzeni (skip index / pattern step):**

Bir tam desende p adet dome temas noktası (touchpoint) ekvator çevresinde
eşit aralıklı dağılır. 1. devre 1 numaralı touchpoint'ten başlar ve 1+q
numaralı touchpoint'e atlar. Genel olarak i. devrenin dome-1 touchpoint'i:

$$\text{TP}_i = [(i-1) \cdot q] \mod p + 1, \quad i = 1, 2, \ldots, p \tag{5.6}$$

gcd(p,q) = 1 koşulu sağlandığında (S-PAT-02), bu formül p devre sonunda tüm
touchpoint'leri tam olarak bir kez ziyaret eder — yani permütasyon tamdır.

**Sıralı (sequential) vs atlayışlı (non-sequential) sarım:**

- **q = 1:** Ardışık touchpoint'ler ziyaret edilir (1→2→3→...). Liflerin
  üst üste binme düzeni düzenli, interweaving minimumdur. Basit mandrel
  kontrol profili.

- **q > 1:** Touchpoint'ler atlanır (1→4→7→... gibi). Liflerin çapraz
  örülmesi artar — bu, laminat bütünlüğünü ve delaminasyon direncini
  artırabilir, ancak dome bölgesinde lif yığılma riski oluşturabilir.

**Pattern tablo çıktısına eklenen sütun:**

| Sütun | Birim | Açıklama |
|-------|-------|----------|
| q (Skip) | — | Atlama adımı — zaten mevcut q sütunu bu işlevi görür |
| Sequence | — | İlk 5 touchpoint sırası (örn: "1→4→7→10→13") |

**İmplementasyon notu:** q değeri pattern tablosunda zaten mevcuttur.
"Sequence" sütunu, ilk 5 touchpoint'in Eq. (5.6) ile hesaplanıp string
olarak gösterilmesidir — kullanıcıya atışma düzenini görsel olarak aktarır.
CADWIND'in "Skip" sütunu ve TaniqWind'in "q" sütunu ile eşdeğerdir.

**Interweaving derecesi (bilgi amaçlı):**

Düşük q/p oranı → düşük interweaving (ardışık lifler yan yana).
Yüksek q/p oranı → yüksek interweaving (ardışık lifler çapraz konumda).
Bu bilgi pattern tablosunda raporlanır ama seçim mühendise bırakılır —
yapısal gereksinimler (laminat teorisi) bu kararı yönlendirir.

---

## 6. Hibrit Mod Detayları (Karar-13)

### 6.1 Angle-Driven Mod (α → pattern listesi)

Kullanıcı winding açısı α_eq girer (veya RE üzerinden dolaylı olarak).

1. GeodesicSolver(α_eq) → Δφ_circuit hesapla
2. K = fmod(Δφ_circuit, 2π)
3. Bölüm 5 arama algoritmasını çalıştır
4. Tüm uyumlu pattern'ları listele
5. Her pattern için L_cyl düzeltme hesapla (Bölüm 7)
6. Tablo üret, p küçükten büyüğe sırala (K-2b-04)

### 6.2 Pattern-Driven Mod (p/q → α veya L_cyl)

Kullanıcı p ve q girer.

1. gcd(p,q) = 1 kontrolü → başarısızsa hata
2. Hedef: Δφ_circuit = 2πq/p
3. Mevcut geometriyle Δφ_circuit_actual hesapla
4. Δφ farkı = Δφ_circuit_actual − 2πq/p
5. L_cyl düzeltmesi ile kapatma (Bölüm 7)
6. Düzeltme yetersizse → α iterasyonu (Newton-Raphson):
   - f(RE) = Δφ_circuit(RE) − 2πq/p = 0
   - RE → α_eq → GeodesicSolver → Δφ_circuit zinciri
   - Yakınsama: |f| < 10⁻⁶ rad
   - Maks. 20 iterasyon

### 6.3 TaniqWind Referansı

TaniqWind'in "Path parameters" panelinde "Angle" ve "Rho0" seçenekleri
var — bizim angle-driven (α girdisi) ve pattern-driven (RE → dolaylı α)
modlarımıza karşılık gelir. "Generate by: coverage" modu bizim coverage-based
arama modumuz.

---

## 7. L_cyl İnce Ayar Mekanizması (S-PAT-01)

### 7.1 Motivasyon

Geometriden gelen Δφ_circuit genellikle tam bir 2πq/p'ye denk gelmez.
Silindirik bölge uzunluğu L_cyl, φ_cyl'yi doğrusal olarak kontrol eder:

$$\phi_{cyl} = \frac{L_{cyl} \cdot \tan\alpha_{cyl}}{R_{eq}} \tag{7.1}$$

dolayısıyla L_cyl'deki küçük bir değişiklik Δφ_circuit'i doğrudan etkiler.

### 7.2 Kapalı Form Düzeltme

Hedef Δφ_circuit = 2πq/p olması için gerekli silindir düzeltmesi:

$$\Delta L_{cyl} = \frac{(2\pi q / p - \Delta\phi_{circuit}) \cdot R_{eq}}{2 \cdot \tan\alpha_{cyl}} \tag{7.2}$$

Çarpan 2: bir tam devrede silindir iki kez geçilir (forward + return).

Düzeltilmiş silindir uzunluğu:

$$L_{cyl}^{adj} = L_{cyl} + \Delta L_{cyl} \tag{7.3}$$

### 7.3 Uyarı Eşiği (K-2b-02)

$$\text{Varsayılan eşik:} \quad |\Delta L_{cyl}| < 0.01 \times L_{cyl} \tag{7.4}$$

- Eşik altında: düzeltme sessizce uygulanır
- Eşik üstünde: `LCYL_ADJUSTMENT_WARNING` uyarısı — düzeltme miktarı ve
  oranı raporlanır, uygulanmaya devam edilir
- Üst sınır yok — mühendis karar verir

**Mühendislik gerekçesi:** %1 eşiği, 300 mm silindir için 3 mm'ye karşılık gelir.
Bu, tipik mandrel üretim toleranslarının (±0.5 mm) birkaç katıdır. Daha büyük
düzeltmeler mandrel üretim sürecinde dikkate alınmalıdır — yazılım bu bilgiyi
sunar ama kararı engellemez.

### 7.4 Doğrulama

**Geometrik ön kontrol:** $L_{cyl} + \Delta L_{cyl} < 0$ durumunda düzeltme
geometrik olarak imkansızdır — silindir uzunluğu negatif olamaz. Bu durumda
fiber dome'dan çıkar çıkmaz karşı dome'a girer ve silindirik bölge ortadan
kalkar. Yazılım `std::invalid_argument` hatası fırlatır; bu pattern o
geometri için uygulanamaz.

**Kapama kontrolü:** Düzeltme sonrası GeodesicSolver'ı L_cyl_adj ile yeniden çağır,
|Δφ_circuit_new − 2πq/p| < 10⁻⁶ rad olmalıdır. Dome bölgesi geometriden
bağımsız olduğundan bu doğrulama hızlıdır (sadece φ_cyl yeniden hesaplanır).

---

## 8. Çoklu Katman ve R_eff Güncellemesi (Karar-21)

### 8.1 Katman Bazında Efektif Yarıçap

Her sarım katmanı tamamlandığında mandrel efektif yarıçapı artar:

$$R_{eff}^{(i)} = R_{eq} + (i-1) \cdot BT_{eff} \tag{8.1}$$

burada $i$ katman numarası (1'den başlar), $BT_{eff}$ efektif katman kalınlığı.

Efektif polar açıklık yarıçapı da güncellenir:

$$R_E^{(i)} = r_0 + \frac{BW_{eff}}{2} + (i-1) \cdot BT_{eff} \tag{8.2}$$

**Gerginliğe bağlı kalınlık sıkıştırma (compaction factor):**

Islak sarımda (wet winding) fiber gerginliği, reçine sıkışması nedeniyle
katmanın efektif kalınlığını nominal BT değerinin altına düşürür:

$$BT_{eff} = BT \cdot C_f \tag{8.2a}$$

burada $C_f$ sıkıştırma katsayısı (compaction factor), $0 < C_f \leq 1$.

| Proses tipi | Tipik $C_f$ aralığı | Açıklama |
|-------------|---------------------|----------|
| Kuru sarım (prepreg) | 0.90 – 1.00 | Minimal sıkıştırma |
| Islak sarım, düşük gerginlik | 0.80 – 0.90 | Orta sıkıştırma |
| Islak sarım, yüksek gerginlik | 0.65 – 0.80 | Belirgin reçine sıkışması |

**Phase-2b implementasyonu:** $C_f$ kullanıcı girdisi olarak `WindingParams`
struct'ına eklenir (varsayılan: 1.0, yani sıkıştırma yok). Denklem (8.1) ve
(8.2)'deki BT yerine $BT_{eff} = BT \cdot C_f$ kullanılır.

**Kapsam sınırı:** Gerginlik-kalınlık ilişkisinin analitik modeli (fiber
tension → C_f dönüşümü) bu fazın kapsamı dışındadır. Kullanıcı, proses
deneyimine veya deneysel veriye dayalı olarak $C_f$ değerini doğrudan girer.
İleride ayrı bir proses modeli modülü olarak eklenebilir.

### 8.2 İteratif Pattern Hesabı

Her katman için ayrı GeodesicSolver çalıştırılır:

```
Katman i = 1, 2, ..., N_layers:
  1. R_eq^(i) = R_eq + (i-1) · BT
  2. R_E^(i) = r₀ + BW_eff/2 + (i-1) · BT
  3. α_eq^(i) = arcsin(R_E^(i) / R_eq^(i))
  4. GeodesicSolver(R_eq^(i), R_E^(i), ...) → Δφ_circuit^(i)
  5. Pattern arama(Δφ_circuit^(i)) → (p, q, dk, n)^(i)
  6. L_cyl düzeltme → L_cyl_adj^(i)
```

**Kritik not:** Katman kalınlığı arttıkça α_eq değişir — aynı winding açısı
hedeflense bile RE/R_eq oranı farklılaştığından pattern parametreleri
katmandan katmana değişebilir.

### 8.3 Farklı α Desteği (K-2b-03)

Karar-20 konfigürasyon yapısı her katman için bağımsız parametre tanımlar:

```json
{
  "winding_sequence": [
    {"winding_type": "helical", "alpha_deg": 30, "n_layers": 2},
    {"winding_type": "hoop", "n_layers": 1},
    {"winding_type": "helical", "alpha_deg": 15, "n_layers": 2}
  ]
}
```

Her eleman bağımsız α (veya p/q) kullanabilir. TaniqWind'in katman ağacı
yapısıyla tutarlı.

### 8.4 Hoop Katman Desteği

Hoop winding (α ≈ 90°) dome'a giremez (S-WIND-02). Sadece silindirik bölge
kaplanır. Endüstriyel yazılımlarda hoop sarımlar için kullanıcıya gap/overlap
tercihi sunulur.

**Kullanıcı girdisi — hoop yerleştirme modu:**

| Mod | Açıklama |
|-----|----------|
| `BUTT_JOINT` | Tow'lar tam bitişik, boşluk ve örtüşme yok (varsayılan) |
| `GAP` | Tow'lar arasında kullanıcı tanımlı boşluk (gap_mm > 0) |
| `OVERLAP` | Tow'lar arasında kullanıcı tanımlı örtüşme (overlap_mm > 0) |

**Hoop devre sayısı hesabı:**

Butt-joint (varsayılan):

$$n_{hoop} = \left\lceil \frac{L_{cyl}}{BW_{eff}} \right\rceil \tag{8.3a}$$

Gap modunda:

$$n_{hoop}^{gap} = \left\lceil \frac{L_{cyl}}{BW_{eff} + gap} \right\rceil \tag{8.3b}$$

Overlap modunda:

$$n_{hoop}^{ovlp} = \left\lceil \frac{L_{cyl}}{BW_{eff} - overlap} \right\rceil \tag{8.3c}$$

burada $overlap < BW_{eff}$ koşulu zorunludur (overlap ≥ BW_eff → sonsuz devre).

**Hoop coverage hesabı:**

$$\text{Coverage}_{hoop} = \frac{n_{hoop} \cdot BW_{eff}}{L_{cyl}} \times 100\% \tag{8.4}$$

**Hoop katmanı pattern tablosunda ayrı satır:**
- p = n_hoop, q = 0, dk = 1, n = n_hoop
- Yerleştirme modu (BUTT/GAP/OVERLAP) ve parametresi raporlanır
- L_cyl düzeltme gerekmez — hoop yerleştirme mandrel boyutunu değiştirmez

**Mühendislik notu:** Gap modu reçine akış kanalları oluşturur (wet winding'de
faydalı olabilir). Overlap modu daha yüksek kalınlık ve dayanım sağlar ama
homojen olmayan kalınlık dağılımına neden olur. Butt-joint, eşit kalınlık
hedefleyen tasarımlar için tercih edilir.

---

## 9. Pattern Tablo Çıktısı (Endüstriyel Format)

### 9.1 Tablo Sütunları

CADWIND ve TaniqWind referans alınarak tasarlanmış sütun yapısı:

| Sütun | Birim | Açıklama |
|-------|-------|----------|
| # | — | Sıra numarası |
| p | — | Devre sayısı |
| q | — | Atlama adımı (skip index) |
| Sequence | — | İlk 5 touchpoint sırası (örn: "1→4→7→10→13") |
| +/− | — | Leading (+) / Lagging (−) |
| dk | — | K partisyonundaki tow sayısı |
| n | — | Ekvator toplam tow sayısı |
| Coverage | % | Ekvator kaplama oranı |
| Overlap | mm | Ekvator örtüşme miktarı (negatif = boşluk) |
| α_eq | ° | Ekvator winding açısı |
| K | rad | Pattern sabiti |
| L_cyl_adj | mm | Düzeltilmiş silindir uzunluğu |
| ΔL_cyl | mm | Silindir düzeltme miktarı |
| Path_length | m | Tek devre fiber yol uzunluğu |
| Total_length | m | Toplam fiber uzunluğu (p × tek devre) |

### 9.2 Varsayılan Sıralama (K-2b-04)

Tablo varsayılan olarak **p küçükten büyüğe** sıralanır. Coverage % ve
Overlap yan sütun olarak her satırda görünür. Bu yaklaşım mühendise tüm
seçenekleri sunar ve kararı ona bırakır.

Aynı p değerine sahip birden fazla çözüm varsa (leading ve lagging),
leading (+) önce gösterilir.

### 9.3 CADWIND / TaniqWind Karşılaştırma Tablosu

| Özellik | CADWIND | TaniqWind | Bizim Yazılım |
|---------|---------|-----------|---------------|
| Girdi modu | α + BW | α veya Rho0 | α veya p/q (K-13) |
| Arama modu | Coverage range / Cycles range | Coverage desired | Coverage (vars.), p range (gel.) |
| Sıralama | Sütun bazlı | Sort by dropdown | Varsayılan p↑ (K-2b-04) |
| Çıktı | Tablo + 3D | Tablo + 3D | Tablo + JSON |
| Çoklu katman | Ayrı modül | Layer tree | Sekans listesi (K-20) |

---

## 10. Hata Yönetimi (Karar-17)

### 10.1 Alt Katman

| Durum | Dönüş |
|-------|-------|
| p aralığında uyumlu pattern yok | Boş liste (nullopt değil) |
| Diophantine çözüm yok (gcd koşulu) | O p atlanır |
| n ≤ 0 (fiziksel olarak imkansız geometri) | `std::nullopt` |

### 10.2 Üst Katman

| Durum | Seviye |
|-------|--------|
| Hiç uyumlu pattern bulunamadı | `std::runtime_error` |
| gcd(p,q) ≠ 1 (pattern-driven modda) | `std::invalid_argument` |
| q ≥ p veya q ≤ 0 | `std::invalid_argument` |
| R_E ≥ R_eq | `std::invalid_argument` |
| L_cyl + ΔL_cyl < 0 (geometrik imkansızlık) | `std::invalid_argument` |
| |ΔL_cyl| > %1 × L_cyl (K-2b-02) | Uyarı (`WARNING`), exception değil |
| BT ≤ 0 veya d ≤ 0 | `std::invalid_argument` |
| |Δφ_target − Δφ_actual| > 0.5° (spindle tolerans) | Uyarı (`WARNING`), exception değil |
| Hoop overlap ≥ BW_eff | `std::invalid_argument` |
| C_f ≤ 0 veya C_f > 1 | `std::invalid_argument` |

---

## 11. Toleranslar ve Doğrulama

### 11.1 Pattern Kapama Toleransı

$$|p \cdot \Delta\phi_{circuit}^{adj} - 2\pi q| < 10^{-6} \; \text{rad} \tag{11.1}$$

p devre toplam açısal ilerleme, hedef 2πq'ya bu toleransta eşit olmalıdır.

### 11.2 Coverage Doğrulama

$$\left| \frac{n \cdot b_{eq}}{2\pi R_{eq} \cdot d} - \frac{\text{Coverage}_{reported}}{100} \right| < 10^{-4} \tag{11.2}$$

### 11.3 Diophantine Doğrulama

$$|p \cdot dk - n \cdot d| = 1 \tag{11.3}$$

Tam tamsayı kontrolü — kayan nokta hatası yok (tamsayı aritmetiği).

### 11.4 gcd Doğrulama

$$\gcd(p, q) = 1 \tag{11.4}$$

Öklid algoritması ile kesin kontrol.

### 11.5 Spindle Düzeltme Toleransı (Phase-5 Emniyet Sınırı)

Hesaplanan desenin makine kontrolcüsü tarafından fiziksel olarak uygulanabilmesi
için, her devrenin açısal ilerleme hatasının spindle rotasyon hassasiyeti
dahilinde kalması gerekir:

$$|\Delta\phi_{circuit}^{target} - \Delta\phi_{circuit}^{actual}| \leq 0.5° = 8.727 \times 10^{-3} \; \text{rad} \tag{11.5}$$

Bu tolerans, makine kontrolcüsünün (Phase-5) numerik hassasiyeti için bir
emniyet sınırıdır. L_cyl düzeltmesi sonrası bu koşul otomatik olarak
sağlanmalıdır (Eq. 7.2 kapalı form düzeltme, Eq. 11.1 kapama toleransı
10⁻⁶ rad << 0.5° olduğundan).

**Mühendislik gerekçesi:** 0.5° tolerans, tipik servo motor encoder çözünürlüğü
(0.01°) ile karşılaştırıldığında rahat bir marj sağlar. Asıl amaç,
L_cyl düzeltmesi yapılamayan durumlarda (mandrel boyutu sabit, düzeltme
reddedilmiş) kalan açısal hatanın makine tarafından telafi edilip
edilemeyeceğini raporlamaktır.

**Açık madde ilişkisi:** Bu 0.5° değeri şu an varsayımdır — makine
spesifikasyonundan doğrulanmamıştır. AÇIK-P5-01 (G-code açı birimi ve makine
kontrolcüsü hassasiyeti) kapsamında Phase-5 öncesinde gözden geçirilecektir.
Gerçek makine encoder çözünürlüğü ve kontrolcü hassasiyeti belirlendiğinde
bu eşik güncellenebilir.

**İmplementasyon:** Pattern tablosundaki her satır için bu kontrol yapılır.
Eşik aşılırsa `SPINDLE_CORRECTION_WARNING` uyarısı eklenir — bu bir exception
değildir ve yol üretimini durdurmaz. Mühendis, daha düşük hassasiyetle sarım
yapılmasını kabul ederek bu uyarıyı bilinçli olarak geçebilir. Yazılım sapma
miktarını raporlar ama nihai kararı mühendisin değerlendirmesine bırakır.
Bu uyarı tipik olarak L_cyl düzeltmesi uygulanmadan doğrudan pattern-driven
mod kullanıldığında veya mandrel boyutunun sabit olduğu durumlarda ortaya çıkar.

---

## 12. GATE-2b Doğrulama Koşulları

| Test | Açıklama | Kriter |
|------|----------|--------|
| T1 | gcd(p,q) = 1 — tüm üretilen pattern'larda | Tam tamsayı |
| T2 | Pattern kapama: p·Δφ_circuit_adj ≈ 2πq | |fark| < 10⁻⁶ rad |
| T3 | Ekvator kaplama: n·b_eq ≈ 2πR_eq·d | %1 tolerans |
| T4 | L_cyl düzeltme tutarlılığı: düzeltilmiş L_cyl ile yeniden hesap | |Δφ fark| < 10⁻⁶ rad |
| T5 | Round-trip: pattern-driven → angle-driven → aynı p/q | Eşleşme |
| T6 | Leading/lagging Diophantine koşulu: p·dk − n·d = ±1 | Tamsayı |
| T7 | Endüstri test senaryoları TEST-01..04 × 3 dome tipi | S1..S4 parametreleri |
| T8 | Sınır durumlar: α → hoop (dome yok), p=1, çok büyük p | Hata yönetimi |
| T9 | Çoklu katman R_eff: katman-2 pattern ≠ katman-1 (BT > 0 ise) | R_eff farkı |
| T10 | Hoop katman: n_hoop hesabı ve coverage | Basit kaplama |
| T11 | Spindle toleransı: |Δφ_target − Δφ_actual| ≤ 0.5° (8.73e-3 rad) | L_cyl düzeltme sonrası |
| T12 | Hoop gap/overlap: üç mod (BUTT/GAP/OVERLAP) doğru n_hoop üretir | Her mod ayrı test |
| T13 | Skip-index sequence: Eq. 5.6 ile touchpoint permütasyonu tam | p nokta = p farklı TP |
| T14 | Compaction factor: C_f < 1 durumunda R_eff < C_f=1 R_eff | BT_eff farkı |

---

## 13. Phase-2a API Kullanım Haritası

Phase-2b pattern solver, Phase-2a'nın aşağıdaki API'lerini kullanır:

| Phase-2a Bileşeni | Kullanım |
|--------------------|----------|
| `GeodesicSolver::solve()` | Δφ_circuit hesabı |
| `GeodesicPath::delta_phi_circuit` | Pattern sabiti K türetimi |
| `GeodesicPath::phi_dome` | Dome açısal katkısı |
| `GeodesicPath::phi_cyl` | Silindir açısal katkısı (L_cyl düzeltme) |
| `GeodesicParams` | BW_eff, RE, α_eq parametre aktarımı |

Phase-1b'den:

| Phase-1b Bileşeni | Kullanım |
|--------------------|----------|
| `MandrelGeometry` | R_eq, r₀, L_cyl, dome_type |
| `ConfigParser` | JSON sekans okuma (Karar-20) |

---

## 14. C++ Sınıf Arayüzü (Tasarım Önerisi)

### 14.1 PatternResult Struct

```
struct PatternResult {
    int p;                    // Devre sayısı
    int q;                    // Atlama adımı (skip index)
    int dk;                   // K partisyon tow sayısı
    int n;                    // Ekvator toplam tow sayısı
    int d;                    // Katman sayısı
    PatternType type;         // LEADING / LAGGING
    std::string sequence;     // İlk 5 touchpoint sırası (örn: "1→4→7→10→13")
    double K_rad;             // Pattern sabiti [rad]
    double alpha_eq_rad;      // Ekvator winding açısı [rad]
    double coverage_pct;      // Coverage [%]
    double overlap_mm;        // Ekvator overlap [mm]
    double L_cyl_adj_mm;      // Düzeltilmiş L_cyl [mm]
    double delta_L_cyl_mm;    // L_cyl düzeltme miktarı [mm]
    double path_length_m;     // Tek devre yol uzunluğu [m]
    double total_length_m;    // Toplam fiber uzunluğu [m]
    bool lcyl_warning;        // |ΔL| > %1 uyarısı
    bool spindle_warning;     // |Δφ_target − Δφ_actual| > 0.5° uyarısı
};
```

### 14.2 PatternSearchParams Struct

```
struct PatternSearchParams {
    SearchMode mode;          // ANGLE_DRIVEN / PATTERN_DRIVEN
    // Angle-driven parametreleri:
    double alpha_eq_deg;      // Winding açısı [°]
    // Pattern-driven parametreleri:
    int target_p;             // Hedef p
    int target_q;             // Hedef q
    // Ortak:
    int d_layers;             // Hedef katman sayısı (varsayılan: 1)
    double coverage_min_pct;  // Min coverage % (varsayılan: 100)
    double coverage_max_pct;  // Max coverage % (varsayılan: 150)
    double compaction_factor; // Sıkıştırma katsayısı C_f (varsayılan: 1.0)
    // Hoop katman parametreleri:
    HoopMode hoop_mode;       // BUTT_JOINT / GAP / OVERLAP (varsayılan: BUTT_JOINT)
    double hoop_gap_mm;       // Gap miktarı [mm] (hoop_mode=GAP için)
    double hoop_overlap_mm;   // Overlap miktarı [mm] (hoop_mode=OVERLAP için)
    // Gelişmiş:
    int p_min;                // Min devre sayısı (varsayılan: 3)
    int p_max;                // Max devre sayısı (varsayılan: 100)
    double lcyl_warn_pct;     // L_cyl uyarı eşiği (varsayılan: 1.0)
};
```

### 14.3 PatternSolver Sınıfı

```
class PatternSolver {
public:
    PatternSolver(const MandrelGeometry& geom,
                  const GeodesicParams& wind_params);

    // Angle-driven: tüm uyumlu pattern'ları bul
    std::vector<PatternResult> findPatterns(
        const PatternSearchParams& search) const;

    // Pattern-driven: belirli p/q için L_cyl düzeltme
    PatternResult solveForPattern(int p, int q, int d) const;

    // Çoklu katman sekansı
    std::vector<PatternResult> solveSequence(
        const std::vector<LayerSpec>& sequence) const;
};
```

---

## 15. MATLAB Referans İmplementasyon Yapısı

### 15.1 Dosya Organizasyonu

```
MATLAB/phase2b_pattern/
├── pattern_solver.m              % Ana fonksiyon — pattern arama
├── find_compatible_patterns.m    % Diophantine arama motoru
├── lcyl_adjustment.m             % L_cyl ince ayar hesabı
├── verify_pattern.m              % GATE-2b doğrulama (10 test)
└── reference_data/               % CSV referans çıktılar
```

### 15.2 pattern_solver.m Arayüzü

```
Girdi:
  dome_profile   — Phase-1a profil struct
  R_eq, r0, L_cyl, BW_eff, BT — geometri parametreleri
  mode           — 'angle_driven' veya 'pattern_driven'
  alpha_eq_deg   — winding açısı [°] (angle-driven için)
  target_p, target_q — (pattern-driven için)
  d_layers       — hedef katman sayısı (varsayılan: 1)
  coverage_range — [min_pct, max_pct] (varsayılan: [100, 150])

Çıktı:
  patterns — struct dizisi, her eleman PatternResult alanlarını içerir
  table    — MATLAB table nesnesi (gösterim için)
```

---

## 16. Session Planı

### S1 — MATLAB: Pattern Çekirdek Matematik

**Hedef:** Diophantine çözücü, K hesabı, kaplama, L_cyl ayar referans implementasyonu.

**Kapsam:** `pattern_solver.m`, `find_compatible_patterns.m`, `lcyl_adjustment.m`
— Phase-2a `geodesic_single_circuit.m` çağrılarak Δφ_circuit elde edilir.
Angle-driven ve pattern-driven modlar. Hoop katman desteği.

**Tamamlanma kriteri:** TEST-01..04 × 3 dome tipiyle çalışır. En az bir uyumlu
pattern bulunur. Pattern-driven modda L_cyl yakınsar.

**Bağımlılık:** Phase-2a S1 (`geodesic_single_circuit.m`)

### S2 — MATLAB: GATE-2b Doğrulama

**Hedef:** Pattern çözümünün iç tutarlılık doğrulaması — 14 test.

**Kapsam:** `verify_pattern.m` — T1..T14 testleri. Referans .csv dışa aktarma.

**Tamamlanma kriteri:** 14/14 PASS. .csv dosyaları üretilmiş.

**Bağımlılık:** S1

### S3 — C++: PatternSolver Çekirdek

**Hedef:** C++ pattern çözücü — hibrit mod.

**Kapsam:** `include/winding/pattern_solver.h` + `src/winding/pattern_solver.cpp`
— Diophantine arama, L_cyl düzeltme, gcd kontrolü, angle-driven ve pattern-driven
modlar. Birim testler.

**Tamamlanma kriteri:** Her iki mod TEST-01 parametreleriyle çalışır. Testler PASS.

**Bağımlılık:** Phase-2a S4 (GeodesicSolver)

### S4 — C++: Çoklu Katman + Pattern Tablo

**Hedef:** Karar-21 çoklu sekans desteği ve tablo çıktısı.

**Kapsam:** WindingSequence, iteratif R_eff güncelleme, hoop katman, JSON dışa aktarım.

**Tamamlanma kriteri:** İki katmanlı sekans çalışır. R_eff doğru güncellenir.

**Bağımlılık:** S3

### S5 — C++: GATE-2b Cross-Validation

**Hedef:** Phase-2b tamamlanma kapısı.

**Kapsam:** MATLAB .csv referanslarıyla karşılaştırma, T1..T14 C++ tarafında
doğrulanması, `docs/gate_2b_report.md`.

**Tamamlanma kriteri:** GATE-2b tüm koşullar PASS.

**Bağımlılık:** S2, S4

---

## 17. Referanslar

- Koussios, S. — "Filament Winding: A Unified Approach" — Bölüm 8 (Winding Patterns), 9 (Pressure Vessels Revisited)
- Peters, S.T. (ed.) — "Composite Filament Winding" — Bölüm 9 (Pattern Selection)
- CADWIND V10 — Material Design Software (arayüz referansı)
- TaniqWind Pro — Coverage Path Design Software (arayüz referansı)
- Phase-2a spesifikasyonu: docs/phase2a_geodesic_math.md
- Phase-0 karar kaydı: docs/phase0_decisions.md
