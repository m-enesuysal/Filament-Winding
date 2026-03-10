# Phase-2a: Geodesic Sarım Yolu Matematik Spesifikasyonu

**Proje:** Filament Winding Yazılımı — 4-Eksen Lathe-Type COPV Üretim Makinesi  
**Faz:** Phase-2a (Geodesic Path Matematiği)  
**Durum:** Aktif  
**Tarih:** 2026-03-08  
**Repo:** https://github.com/m-enesuysal/Filament-Winding

---

## 1. Amaç ve Kapsam

Bu doküman, Phase-2a geodesic sarım yolu hesaplama modülünün matematiksel
temellerini tanımlar. Hedef: tek devre (single circuit) geodesic fiber yolunun
tam mandrel üzerinde üretilmesi — dome-1 dönüşünden dome-2 dönüşüne ve geri.

**Kapsam içi:**
- Clairaut ilişkisi ve winding açısı dağılımı
- Dome bölgesinde paralel açı (φ) entegrasyonu
- Silindirik bölgede helikal ilerleme
- Dome dönüş noktası (turnaround) singülarite yönetimi
- Tek devre toplam açısal ilerleme (Δφ_circuit) hesabı
- Hoop winding kısıtı (S-WIND-02)
- Konveksite kontrolü ve lif köprülemesi riski (S-WIND-04)
- Geodeziklikten sapma emniyet marjı (S-WIND-05)
- Adaptif ODE çıktısından uniform ızgara yeniden örnekleme stratejisi
- Tek katman, sabit R_eq

**Kapsam dışı:**
- Sonlu bant genişliği etkisi (finite bandwidth, Koussios Bölüm 9 ε parametresi)
- Çoklu sekans / R_eff güncelleme (Phase-2b'ye ertelendi)
- Pattern generation (Phase-2b)
- Non-geodesic yollar

---

## 2. Referans Karar Kayıtları

| Karar | Açıklama | Etki |
|-------|----------|------|
| Karar-2 | Simetrik mandrel — ön ve arka dome özdeş | φ_dome1 = φ_dome2 |
| Karar-4 | Silindirik (ρ, φ, x) birincil koordinat sistemi | Tüm formülasyon |
| Karar-5 | RE = r₀ + BW_eff/2, efektif polar açıklık yarıçapı | Clairaut sabiti |
| Karar-7 | α: mandrel ekseni ile fiber arası açı (Konvansiyon A) | sin(α) formülasyonu |
| Karar-9 | Pre-computed lookup table, kübik spline interpolasyon | Dome metrik sorgusu |
| Karar-11 | Üç katmanlı tolerans sistemi | Numerik doğrulama |
| Karar-12 | BW_eff = N_tow × BW | RE hesabı |
| Karar-13 | Phase-2 = Phase-2a + Phase-2b | Faz bölünmesi |
| Karar-17 | Alt katman: std::optional, üst katman: exception | Hata yönetimi |

---

## 3. Koordinat Sistemi ve Notasyon

### 3.1 Koordinat Tanımları

Mandrel ekseni $x$ yönündedir. Silindirik koordinat sistemi $(\\rho, \\phi, x)$:

- $\\rho$: Mandrel eksenine dik mesafe (yarıçap) [mm]
- $\\phi$: Paralel açı (azimuthal açı, mandrel ekseni etrafında) [rad]
- $x$: Mandrel ekseni boyunca konum [mm]

### 3.2 Yay Uzunluğu Parametrizasyonu

Phase-1b'den miras alınan global yay uzunluğu sistemi:

$$s_{global} \\in [0,\\; s_{dome} + L_{cyl} + s_{dome}]$$

- $s \\in [0, s_{dome}]$: Dome-1 bölgesi (ekvatordan kutba doğru)
- $s \\in [s_{dome}, s_{dome} + L_{cyl}]$: Silindirik bölge
- $s \\in [s_{dome} + L_{cyl}, 2 s_{dome} + L_{cyl}]$: Dome-2 bölgesi

### 3.3 Sembol Tablosu

| Sembol | Birim | Tanım |
|--------|-------|-------|
| $\\alpha$ | rad | Winding açısı — mandrel ekseni ile fiber arası (Karar-7) |
| $\\rho(s)$ | mm | Meridyen profil yarıçapı (MeridianLookupTable'dan) |
| $\\phi$ | rad | Paralel açı (azimuthal konum) |
| $\\beta(s)$ | rad | Meridyen eğim açısı |
| $\\kappa_m(s)$ | 1/mm | Meridyen eğriliği |
| $R_{eq}$ | mm | Ekvator yarıçapı |
| $r_0$ | mm | Fiziksel polar açıklık (boss) yarıçapı |
| $R_E$ | mm | Efektif polar açıklık yarıçapı = $r_0 + BW_{eff}/2$ |
| $c$ | mm | Clairaut sabiti = $R_E$ |
| $E(\\nu)$ | mm² | Birinci temel form katsayısı — paralel metrik = $\\rho(\\nu)^2$ |
| $G(\\nu)$ | mm² | Birinci temel form katsayısı — meridyen metrik |
| $L_{cyl}$ | mm | Silindirik bölge uzunluğu |
| $\\Delta\\phi_{circuit}$ | rad | Tek devredeki toplam açısal ilerleme |

---

## 4. Clairaut İlişkisi — Geodesic Winding Açısı

### 4.1 Temel Denklem

Bir dönel yüzey (shell of revolution) üzerindeki geodesic eğri, Clairaut
ilişkisiyle tanımlanır:

$$\\rho(s) \\cdot \\sin\\alpha(s) = c = R_E \\tag{4.1}$$

Bu denklem, Euler-Lagrange varyasyonel prensibinden türetilir ve fiziksel
olarak "geodesic yolda, yerel yarıçap ile winding açısının sinüsünün çarpımı
sabittir" anlamına gelir.

### 4.2 Winding Açısı Dağılımı

Clairaut ilişkisinden winding açısı doğrudan çıkar:

$$\\alpha(s) = \\arcsin\\!\\left(\\frac{R_E}{\\rho(s)}\\right) \\tag{4.2}$$

**Geçerlilik koşulu:** $\\rho(s) \\geq R_E$ her noktada sağlanmalıdır. $\\rho < R_E$
durumunda arcsin argümanı 1'i aşar — bu fiziksel olarak fiberin o yarıçapa
ulaşamayacağı anlamına gelir.

### 4.3 Özel Durumlar

**Ekvator noktası** ($\\rho = R_{eq}$):

$$\\alpha_{eq} = \\arcsin\\!\\left(\\frac{R_E}{R_{eq}}\\right) \\tag{4.3}$$

Bu, silindirik bölgedeki sabit helikal winding açısıdır.

**Dome dönüş noktası** ($\\rho = R_E$):

$$\\alpha_{turn} = \\arcsin(1) = \\frac{\\pi}{2} \\tag{4.4}$$

Fiber tam çevresel yönde (hoop yönünde) ilerler ve yön değiştirir.

**Silindirik bölge** ($\\rho = R_{eq} = \\text{sabit}$):

$$\\alpha_{cyl} = \\alpha_{eq} = \\text{sabit} \\tag{4.5}$$

Clairaut otomatik olarak sabit açılı helikal ilerleme üretir.

### 4.4 Hoop Winding Kısıtı (S-WIND-02)

Hoop winding'de $\\alpha \\approx 90°$, dolayısıyla Clairaut sabiti:

$$c = R_{eq} \\cdot \\sin(90°) = R_{eq}$$

Bu durumda $c = R_{eq} > R_E$ olduğundan dome dönüş noktası $\\rho = c = R_{eq}$
ekvator noktasının kendisidir — fiber dome bölgesine giremez.

**İmplementasyon kuralı:** Sarım tipi `HOOP` seçildiğinde dome yolu üretilmez.
Hoop winding yalnızca silindirik bölgede uygulanır. Bu kontrol geodesic
solver giriş noktasında zorunlu olarak yapılmalıdır.

---

## 5. Dome Bölgesinde ODE Sistemi

### 5.1 Paralel Açı Diferansiyel Denklemi

Dome bölgesinde fiber yolunun paralel açı ($\\phi$) ilerlemesini bulmak için
birinci temel form katsayılarını kullanırız. Dönel yüzey metrikleri:

$$E(\\nu) = \\rho(\\nu)^2, \\quad F(\\nu) = 0, \\quad G(\\nu) = \\left(\\frac{d\\rho}{d\\nu}\\right)^2 + \\left(\\frac{dx}{d\\nu}\\right)^2 \\tag{5.1}$$

burada $\\nu$ meridyen profil parametresidir. $F = 0$ dönel simetriyi ifade eder.

Bir eğrinin dönel yüzey üzerindeki paralel açı ilerleme hızı (Koussios Eq. 2.39):

$$\\frac{d\\phi}{d\\nu} = \\tan\\alpha(\\nu) \\cdot \\frac{\\sqrt{G(\\nu)}}{\\sqrt{E(\\nu)}} \\tag{5.2}$$

### 5.2 Yay Uzunluğu Parametrizasyonu

MeridianLookupTable $s$ (yay uzunluğu) parametrizasyonunda çalıştığından,
$\\nu = s$ alındığında $G(s) = 1$ olur (yay uzunluğu parametrizasyonunun doğal
özelliği: $ds^2 = d\\rho^2 + dx^2$ ve $\\sqrt{G} \\cdot d\\nu = ds$).

Bu durumda paralel açı denklemi sadeleşir:

$$\\frac{d\\phi}{ds} = \\frac{\\tan\\alpha(s)}{\\rho(s)} \\tag{5.3}$$

Clairaut ilişkisinden $\\tan\\alpha$'yı açalım:

$$\\sin\\alpha = \\frac{R_E}{\\rho}, \\quad \\cos\\alpha = \\sqrt{1 - \\frac{R_E^2}{\\rho^2}} = \\frac{\\sqrt{\\rho^2 - R_E^2}}{\\rho}$$

$$\\tan\\alpha = \\frac{\\sin\\alpha}{\\cos\\alpha} = \\frac{R_E}{\\sqrt{\\rho^2 - R_E^2}} \\tag{5.4}$$

Denklem (5.3)'e yerleştirerek:

$$\\boxed{\\frac{d\\phi}{ds} = \\frac{R_E}{\\rho(s) \\cdot \\sqrt{\\rho(s)^2 - R_E^2}}} \\tag{5.5}$$

Bu, Phase-2a'nın **çekirdek ODE denklemidir**. Sağ taraf yalnızca $\\rho(s)$
fonksiyonuna bağlıdır ve $\\rho(s)$ MeridianLookupTable'dan kübik spline
interpolasyonla elde edilir.

### 5.3 Entegrasyon Yönü ve Başlangıç Koşulu

Entegrasyon ekvatordan ($s = 0$) polar açıklığa ($s = s_{dome}$) doğru yapılır.

**Başlangıç koşulu:**

$$\\phi(s = 0) = \\phi_0 \\tag{5.6}$$

burada $\\phi_0$ devrenin başlangıç paralel açısıdır. İlk devre için $\\phi_0 = 0$
alınabilir.

### 5.4 Entegral Formu

Denklem (5.5) ayrıca kapalı entegral formunda yazılabilir:

$$\\phi_{dome} = \\int_0^{s_{turn}} \\frac{R_E}{\\rho(s) \\cdot \\sqrt{\\rho(s)^2 - R_E^2}} \\, ds \\tag{5.7}$$

burada $s_{turn}$, dome dönüş noktasına karşılık gelen yay uzunluğu değeridir
($\\rho(s_{turn}) = R_E$). Bu integral, dönüş noktasında integrand'ın sonsuza
ıraksaması nedeniyle **uygunsuz (improper) integraldir** — Bölüm 6'daki
singülarite yönetimi zorunludur.

---

## 6. Dome Dönüş Noktası — Singülarite Yönetimi (S-WIND-01)

### 6.1 Singülaritenin Kaynağı

$\\rho \\to R_E$ durumunda denklem (5.5)'teki payda sıfıra yaklaşır:

$$\\rho^2 - R_E^2 \\to 0 \\quad \\Rightarrow \\quad \\frac{d\\phi}{ds} \\to +\\infty \\tag{6.1}$$

Bu, fiberin dönüş noktasında çevresel yönde çok hızlı ilerlediğini ifade eder
— fiziksel olarak geçerli, numerik olarak sorunlu.

### 6.2 Çözüm Stratejisi: Analitik Limit Tamamlama

**Onaylanmış karar:** Numerik entegrasyon $\\rho = R_E + \\varepsilon$'da durdurulur,
kalan $\\phi$ katkısı analitik olarak hesaplanır.

**Adım 1 — Numerik entegrasyon:** Denklem (5.5), $s \\in [0, s_\\varepsilon]$ aralığında
sayısal olarak integre edilir. $s_\\varepsilon$, $\\rho(s_\\varepsilon) = R_E + \\varepsilon$
koşulunu sağlayan yay uzunluğu değeridir.

**Adım 2 — $s_\\varepsilon$ tespiti:** MeridianLookupTable üzerinde ters sorgu ile bulunur:

$$\\rho(s_\\varepsilon) = R_E + \\varepsilon \\tag{6.2}$$

**Adım 3 — Analitik tamamlama:** Dönüş noktası yakınında $\\rho(s)$ profili
lokal olarak doğrusallaştırılır. $\\rho \\approx R_E$ bölgesinde meridyen eğim açısı
$\\beta$ küçük olduğundan (fiberin polar açıklığa yakın dönüş noktasında meridyen
neredeyse dikey):

$$\\rho(s) \\approx R_E + |\\dot{\\rho}(s_{turn})| \\cdot (s_{turn} - s) \\tag{6.3}$$

burada $\\dot{\\rho} = d\\rho/ds$ MeridianLookupTable'dan interpolasyonla elde
edilir. $|\\dot{\\rho}_{turn}|$ dönüş noktasındaki meridyen yarıçap değişim hızıdır.

Kalan $\\phi$ katkısı için $u = \\rho - R_E$ değişken dönüşümü uygulayarak:

$$\\Delta\\phi_{tail} = \\int_{s_\\varepsilon}^{s_{turn}} \\frac{R_E}{\\rho \\sqrt{\\rho^2 - R_E^2}} \\, ds \\tag{6.4}$$

Doğrusallaştırma altında $\\rho \\approx R_E + |\\dot{\\rho}| \\cdot (s_{turn} - s)$ ve
$\\rho^2 - R_E^2 \\approx 2 R_E \\cdot (\\rho - R_E)$ yazılırsa:

$$\\Delta\\phi_{tail} \\approx \\int_0^{\\varepsilon} \\frac{R_E}{(R_E + u)\\sqrt{2 R_E \\cdot u}} \\cdot \\frac{du}{|\\dot{\\rho}_{turn}|} \\tag{6.5}$$

$R_E + u \\approx R_E$ yaklaşımıyla:

$$\\boxed{\\Delta\\phi_{tail} \\approx \\frac{1}{|\\dot{\\rho}_{turn}|} \\cdot \\sqrt{\\frac{2\\varepsilon}{R_E}}} \\tag{6.6}$$

Bu kapalı form ifadesi, singüler bölgedeki $\\phi$ katkısını doğrudan hesaplar.

### 6.3 ε Parametresi Seçimi

$$\\varepsilon = 10 \\times \\text{tol}_{interp} = 10 \\times 10^{-4} \\; \\text{mm} = 10^{-3} \\; \\text{mm} \\tag{6.7}$$

Gerekçe: Karar-11 Katman 2 interpolasyon toleransının 10 katı. Yeterince
büyük ki adaptif solver verimsiz küçük adımlara zorlanmasın, yeterince küçük
ki doğrusallaştırma hatası ihmal edilebilir olsun.

**MATLAB doğrulama gereksinimi:** $\\varepsilon \\in \\{10^{-4}, 10^{-3}, 10^{-2}, 10^{-1}\\}$
aralığında $\\phi_{dome}$ yakınsama testi yapılmalı. Referans değer en küçük
$\\varepsilon$ ile elde edilen sonuçtur.

### 6.4 $\\dot{\\rho}_{turn}$ Hesabı

MeridianLookupTable, $d\\rho/ds$ türevini kübik spline interpolasyonla sağlar.
Dönüş noktasında $\\rho = R_E$ olduğu yerde bu değer sorgulanır.

**Özel durum kontrolü:** İzotensoid dome'da dönüş noktası yakınında $|\\dot{\\rho}|$
çok küçük olabilir (S-GEO-01 stiff davranış). $|\\dot{\\rho}_{turn}| < 10^{-8}$
durumunda ikinci dereceden yaklaşıma geçilmeli veya hata döndürülmeli.

---

## 7. Silindirik Bölgede Geodesic Yol

### 7.1 Sabit Açılı Helikal İlerleme

Silindirik bölgede $\\rho = R_{eq} = \\text{sabit}$, dolayısıyla Clairaut:

$$\\alpha_{cyl} = \\arcsin\\!\\left(\\frac{R_E}{R_{eq}}\\right) = \\text{sabit} \\tag{7.1}$$

Bu, geodesic yolun silindirde sabit pitch heliksi olduğu anlamına gelir.

### 7.2 Paralel Açı İlerleme

Silindirik bölgede $\\rho = R_{eq}$, dolayısıyla denklem (5.3):

$$\\frac{d\\phi}{dx} = \\frac{\\tan\\alpha_{cyl}}{R_{eq}} \\tag{7.2}$$

burada $x$ mandrel ekseni boyunca konumdur. Silindirik bölgede $ds = dx$ (meridyen
düz olduğundan $d\\rho = 0$, $ds = \\sqrt{d\\rho^2 + dx^2} = dx$).

Toplam silindirik $\\phi$ ilerleme:

$$\\phi_{cyl} = \\frac{L_{cyl} \\cdot \\tan\\alpha_{cyl}}{R_{eq}} \\tag{7.3}$$

### 7.3 Aksiyel İlerleme

Silindirik bölgede fiberin birim çevresel ilerleme başına aksiyel ilerleme
miktarı sabit olan bir heliks oluşur. Heliks adımı (pitch):

$$p_{helix} = \\frac{2\\pi R_{eq}}{\\tan\\alpha_{cyl}} \\tag{7.4}$$

---

## 8. Tek Devre Toplam Açısal İlerleme

### 8.1 Devre Tanımı

Bir tam devre (circuit): fiber dome-1 dönüş noktasından başlar → ekvator-1
üzerinden silindirik bölgeye girer → silindir boyunca ilerler → ekvator-2
üzerinden dome-2'ye girer → dome-2 dönüş noktasına ulaşır → geri döner →
dome-2'den çıkar → silindir boyunca geri ilerler → dome-1'e girer → dome-1
dönüş noktasına ulaşır.

Simetrik mandrel varsayımıyla (Karar-2) dome-1 ve dome-2 katkıları özdeştir.

### 8.2 Toplam Açısal İlerleme Formülü

$$\\boxed{\\Delta\\phi_{circuit} = 4 \\cdot \\phi_{dome} + 2 \\cdot \\phi_{cyl}} \\tag{8.1}$$

burada:

- $\\phi_{dome}$: Tek dome yarısı (ekvatordan dönüş noktasına) paralel açı ilerleme
- Çarpan 4: İki dome × iki geçiş (gidiş + dönüş), her geçişte aynı $\\phi_{dome}$
- $\\phi_{cyl}$: Tek silindir geçişi paralel açı ilerleme ($L_{cyl}$ boyunca)
- Çarpan 2: İki silindir geçişi (gidiş + dönüş)

### 8.3 Detaylı Bileşen Hesabı

$$\\phi_{dome} = \\underbrace{\\int_0^{s_\\varepsilon} \\frac{R_E}{\\rho(s) \\sqrt{\\rho(s)^2 - R_E^2}} \\, ds}_{\\text{numerik}} + \\underbrace{\\Delta\\phi_{tail}}_{\\text{analitik (Eq. 6.6)}} \\tag{8.2}$$

$$\\phi_{cyl} = \\frac{L_{cyl} \\cdot \\tan\\alpha_{cyl}}{R_{eq}} \\tag{8.3}$$

### 8.4 İşaret Konvansiyonu

Forward stroke (dome-1 → dome-2): $\\phi$ artar ($d\\phi > 0$).

Return stroke (dome-2 → dome-1): $\\phi$ yine artar (fiber aynı yönde sarar,
geri hareket mandrel eksenindedir, çevresel ilerleme her zaman pozitif).

Bu nedenle $\\Delta\\phi_{circuit}$ her zaman pozitiftir ve denklem (8.1) doğrudan
toplanır.

---

## 9. Tam Fiber Yolu Veri Yapısı

### 9.1 Çıktı Veri Seti

Tek devre geodesic yolu, eşit yay uzunluğu artışlarıyla örneklenen bir
nokta dizisi olarak üretilir:

| Alan | Birim | Açıklama |
|------|-------|----------|
| $s_i$ | mm | Global yay uzunluğu (mandrel boyunca kümülatif) |
| $\\rho_i$ | mm | Yarıçap (MandrelGeometry'den) |
| $x_i$ | mm | Aksiyel konum (MandrelGeometry'den) |
| $\\phi_i$ | rad | Paralel açı (entegrasyondan) |
| $\\alpha_i$ | rad | Winding açısı (Clairaut'tan, Eq. 4.2) |
| $k_{n,i}$ | 1/mm | Normal eğrilik (Eq. 9.2, S-WIND-04 bridging kontrolü) |

### 9.2 3D Kartezyen Dönüşüm

Görselleştirme ve Phase-3 kinematik girdisi için:

$$X_i = \\rho_i \\cos\\phi_i, \\quad Y_i = \\rho_i \\sin\\phi_i, \\quad Z_i = x_i \\tag{9.1}$$

### 9.3 Yol Sürekliliği Doğrulama Kriterleri

Üretilen yol aşağıdaki koşulları sağlamalıdır (system prompt gereksinimleri):

1. **Tek sürekli eğri:** Dome'dan dome'a süreksizlik yok
2. **Geodesic koşulu:** Her noktada $\\rho_i \\sin\\alpha_i = R_E$ (Karar-11 Katman 2 tol.)
3. **Polar açıklık teğetliği:** Dönüş noktasında $\\alpha \\to \\pi/2$
4. **Simetri:** Forward ve return stroke'lar mandrel ekseni etrafında simetrik
5. **Öz-kesişme yok:** Tek stroke içinde yol kendini kesmez
6. **Dome geçiş pürüzsüzlüğü:** Silindir-dome birleşim noktasında $\\alpha$ süreksizliği yok
7. **Konveksite koşulu (S-WIND-04):** Her noktada normal eğrilik $k_n > 0$ olmalıdır;
   $k_n \\leq 0$ durumunda fiber mandrel yüzeyinden ayrılma riski (bridging) vardır ve
   yazılım `BRIDGING_RISK` uyarısı fırlatmalıdır. Bkz. Bölüm 9.4.
8. **Geodeziklikten sapma emniyet marjı (S-WIND-05):** Sayısal hassasiyet hataları
   veya yüzey düzensizlikleri nedeniyle hesaplanan yolun ideal geodezikten sapması
   ($|\\rho \\sin\\alpha - c|$), mevcut sürtünme kapasitesi $\\mu$ limitleri içinde
   kalmalıdır. Bkz. Bölüm 9.5.

### 9.4 Konveksite Kontrolü ve Lif Köprülemesi (S-WIND-04)

Geodesic yollar yalnızca konveks (dışbükey) yüzeylerde kararlıdır. Konkav
(içbükey) bölgelerde fiber mandrel yüzeyine temas kaybeder ve havada köprü
oluşturur (fiber bridging). Bu durum hem yapısal bütünlüğü hem de sarım
kalitesini bozar.

**Normal eğrilik hesabı:** Fiber yolunun normal eğriliği, Euler eğrilik
denklemiyle (Koussios Eq. 2.15) hesaplanır:

$$k_n(s) = \\kappa_m(s) \\cos^2\\alpha(s) + \\kappa_p(s) \\sin^2\\alpha(s) \\tag{9.2}$$

burada $\\kappa_m(s)$ meridyen eğriliği (MeridianLookupTable'dan), $\\kappa_p(s)$
paralel eğriliktir:

$$\\kappa_p(s) = \\frac{\\sin\\beta(s)}{\\rho(s)} = \\frac{-d\\rho/ds}{\\rho(s) \\cdot \\sqrt{1 - (d\\rho/ds)^2}} \\cdot \\frac{1}{\\rho(s)} \\tag{9.3}$$

Her iki eğrilik de Phase-1b MeridianLookupTable'dan mevcut verilerle
($\\rho$, $d\\rho/ds$, $\\kappa_m$) türetilebilir.

**İmplementasyon kuralı:** Yol üretimi sırasında her nokta için $k_n$ hesaplanır.
$k_n \\leq 0$ tespit edildiğinde:
- Yol üretimi durdurulmaz (sadece uyarı)
- İlgili noktanın $s$ değeri ve $k_n$ büyüklüğü kayıt altına alınır
- Çıktı struct'ına `bridging_risk_indices` alanı eklenir
- Kullanıcı seviyesinde `BRIDGING_RISK` uyarı mesajı fırlatılır

**Pratik not:** Desteklenen dome geometrileri (hemispherical, ellipsoidal,
isotensoid) ve silindirik bölge için $k_n > 0$ koşulu teorik olarak her zaman
sağlanır — çünkü bu yüzeyler her yerde konvekstir veya düzdür ($k_n = 0$
silindirik bölgede, $\\alpha \\neq 0$ olduğu sürece). Bununla birlikte, sayısal
hatalardan kaynaklanan marjinal $k_n \\leq 0$ değerleri kontrol edilmeli ve ileride
desteklenebilecek özel geometrilere (örn. konkav geçiş bölgeleri) karşı
kernel savunmalı olmalıdır.

### 9.5 Geodeziklikten Sapma Emniyet Marjı (S-WIND-05)

Teorik olarak geodesic yolda sürtünme gereksinimi sıfırdır ($k_g = 0$).
Ancak pratikte sayısal hassasiyet hataları, mandrel yüzey pürüzlülüğü ve
tow yerleştirme mekanizmasının hassasiyeti nedeniyle fiber yolu ideal
geodezikten sapabilir. Bu sapma, fiberin kaymasına (slippage) neden olur.

Kayma eğilimi (slippage tendency, Koussios Eq. 6.17):

$$|\\lambda| = \\left|\\frac{k_g}{k_n}\\right| \\tag{9.4}$$

Fiberin kaymaması için $|\\lambda| < \\mu$ koşulu sağlanmalıdır; $\\mu$ mandrel
yüzeyindeki statik sürtünme katsayısıdır.

Hesaplanan geodesic yolda $k_g = 0$ olsa da, sayısal Clairaut sapması
($\\delta c = |\\rho \\sin\\alpha - R_E|$) eşdeğer bir yapay $k_g$ üretir. Bu sapmanın
üretilebilirlik emniyet marjı içinde kalması kontrol edilmelidir:

$$\\frac{\\delta c}{\\rho \\cdot k_n} < \\mu_{available} \\tag{9.5}$$

**Tipik $\\mu$ değerleri** (Peters Tablo 9.7):
- Islak sarım (wet winding), karbon/epoksi: $\\mu \\leq 0.11$
- Kuru sarım (prepreg), karbon/epoksi: $\\mu \\leq 0.20$
- Islak sarım, cam/epoksi: $\\mu \\leq 0.18$

**İmplementasyon kuralı:** Winding_type parametresine (Karar-12) göre $\\mu_{max}$
otomatik atanır. GATE-2a'da raporlanır; eşik aşımı uyarı üretir, red etmez.

---

## 10. Dome Tipi Bazında Özel Formüller

### 10.1 Hemispherical Dome

Yarıçap profili: $\\rho(\\theta) = R_{eq} \\cos\\theta$, $\\theta \\in [0, \\theta_{turn}]$

Dönüş noktası açısı: $\\theta_{turn} = \\arccos(R_E / R_{eq})$

Paralel açı ilerleme (Koussios Eq. 5.24 ile uyumlu, k = 1 küre durumu):

$$\\phi_{dome}^{hemi} = \\int_0^{\\theta_{turn}} \\frac{R_E}{R_{eq} \\cos\\theta \\cdot \\sqrt{R_{eq}^2 \\cos^2\\theta - R_E^2}} \\cdot R_{eq} \\, d\\theta \\tag{10.1}$$

Küre üzerinde geodesic bir büyük daire olduğundan, $\\phi_{dome}^{hemi}$ yarım
dönüş açısı $\\pi/2$'dir — bu, sphere ile geodesic kaplama yapılamamasının
nedenidir (Koussios Bölüm 5.4).

**Doğrulama notu:** $R_E / R_{eq} = \\sin\\alpha_{eq}$ küçük olduğunda ($\\alpha_{eq}
\\lesssim 15°$, polar winding bölgesi) $\\phi_{dome}^{hemi} \\approx \\pi/2$ olur.
$R_E / R_{eq} \\to 1$ durumunda ($\\alpha_{eq} \\to 90°$, hoop limiti) dönüş noktası
ekvator'a yaklaşır ve $\\phi_{dome} \\to 0$.

### 10.2 Ellipsoidal Dome

Yarıçap profili: $\\rho(\\theta) = R_{eq} \\cos\\theta$, $x(\\theta) = k \\cdot R_{eq} \\sin\\theta$

burada $k$ aspect ratio parametresidir (Karar-5).

Paralel açı ilerleme Koussios Eq. 5.23-5.26 ile verilir. Genel elipsoid
için analitik çözüm karmaşıktır; sayısal entegrasyon tercih edilir.

Koussios'un turn-around angle yaklaşım polinomu (Eq. 5.25-5.26)
kullanılabilir ama **sadece doğrulama referansı olarak** — implementasyon
sayısal entegrasyon yapmalı, çünkü bu polinom $0.05 \\leq \\varrho \\leq 0.85$
ve $0.21 \\leq k$ aralığında %99.95 doğruluk sağlar ama aralık dışında geçersizdir.

### 10.3 Isotensoid Dome

Profil, Koussios $q$-$r$ parametrizasyonuyla tanımlıdır. Karar-2 uyarınca
simetrik mandrel, $r = 0$.

Paralel açı ilerleme, Koussios Eq. 4.22'den eliptik integral olarak verilir:

$$\\phi_{dome}^{iso}(q, r) = \\frac{\\sqrt{2q - rq + 1}}{\\sqrt{(2q+1)q}} \\cdot \\Pi\\!\\left(\\frac{q}{(2q+1)(1-rq)}; \\frac{1}{q+1}\\right) \\tag{10.2}$$

burada $\\Pi(n; m)$ üçüncü tür tamamlanmamış eliptik integraldir.

Tam devre turn-around açısı (Koussios Eq. 4.28):

$$\\phi_{total}^{iso} = 4 \\cdot \\phi_{dome}^{iso} + \\frac{2 \\cdot s \\cdot \\tan\\alpha_{eq}}{Y_{eq}} \\tag{10.3}$$

burada $s$ boyutsuz silindir uzunluğu, $Y_{eq}$ boyutsuz ekvator yarıçapı.

**İmplementasyon notu:** Eliptik integral çözüm doğrulama referansıdır.
Sayısal entegrasyon (Eq. 5.5) tüm dome tipleri için ortak yoldur ve tercih
edilen implementasyon yöntemidir.

---

## 11. Numerik Entegrasyon Stratejisi

### 11.1 Çözücü Seçimi

Denklem (5.5) birinci dereceden sıradan diferansiyel denklem (ODE):

$$\\frac{d\\phi}{ds} = f(s) = \\frac{R_E}{\\rho(s) \\cdot \\sqrt{\\rho(s)^2 - R_E^2}}$$

Sağ taraf $\\phi$'ye bağlı değildir — bu bir **kuadratür (quadrature) problemi**dir.
$\\phi$ çözüm değişkeni olarak doğrudan entegre edilir; adaptif adım kontrolü
gerekmez çünkü stiffness yoktur.

Ancak dönüş noktası yakınında $f(s) \\to \\infty$ olduğundan, adaptif adım
kontrollü bir çözücü (Karar-6: Boost.Odeint RK45 Dormand-Prince) pratik
avantaj sağlar.

### 11.2 Alternatif: Trapez/Simpson Kuadratürü

$f(s)$ yalnızca bilinen $\\rho(s)$ fonksiyonuna bağlı olduğundan, sabit adımlı
veya adaptif kuadratür yöntemleri de uygulanabilir. Avantajı: Boost.Odeint
bağımlılığı olmadan çalışır, dezavantajı: adaptif adım kontrolünü manuel
yönetmek gerekir.

**Tavsiye:** Phase-1b'de zaten altyapısı hazır olan Boost.Odeint kullanılsın.
Dome yakını adaptif adım doğal olarak küçülür.

### 11.3 Entegrasyon Sınırları

- **Alt sınır:** $s = 0$ (ekvator noktası)
- **Üst sınır:** $s = s_\\varepsilon$ ($\\rho(s_\\varepsilon) = R_E + \\varepsilon$)

$s_\\varepsilon$ bulma: MeridianLookupTable üzerinde ters sorgu (inverse query).
Tablo $s \\to \\rho$ yönünde çalıştığından, $\\rho \\to s$ ters sorgusu gerekir.
Monoton azalan $\\rho(s)$ varsayımı altında (dome bölgesinde her zaman geçerli)
binary search uygulanabilir.

### 11.4 Adım Boyutu Önerisi

İlk adım tahmini:

$$h_0 = \\frac{s_{dome}}{N_{steps}}, \\quad N_{steps} \\geq 500 \\tag{11.1}$$

Adaptif çözücü bu adımı otomatik olarak dome dönüş noktası yakınında
küçültecektir. Karar-11 Katman 1 ODE toleransları uygulanır:

$$\\text{atol} = 10^{-10}, \\quad \\text{rtol} = 10^{-8} \\tag{11.2}$$

---

## 12. Tam Devre Algoritması — Adım Adım

### 12.1 Giriş Parametreleri

| Parametre | Kaynak |
|-----------|--------|
| `MandrelGeometry` nesnesi | Phase-1b çıktısı |
| $R_{eq}$, $r_0$, $L_{cyl}$, dome tipi | Konfigürasyon (Karar-20) |
| $BW$, $N_{tow}$ | Kullanıcı girdisi (Karar-12) |
| Sarım tipi (helical/polar/hoop) | Kullanıcı girdisi |

### 12.2 Ön Hesap

1. $BW_{eff} = N_{tow} \\times BW$
2. $R_E = r_0 + BW_{eff} / 2$
3. $\\alpha_{eq} = \\arcsin(R_E / R_{eq})$
4. **Hoop kontrolü:** $\\alpha_{eq} > 85°$ ise dome yolu üretilmez (S-WIND-02)
5. **Polar uyarı:** $\\alpha_{eq} < 5°$ ise kullanıcıya bilgi (S-WIND-03, red yok)

### 12.3 Dome-1 Yolu (Ekvator → Dönüş Noktası)

1. MeridianLookupTable'dan $s_\\varepsilon$ bul: $\\rho(s_\\varepsilon) = R_E + \\varepsilon$
2. ODE entegrasyonu: $\\phi(s)$ hesapla, $s \\in [0, s_\\varepsilon]$, Eq. (5.5)
3. Analitik tamamlama: $\\Delta\\phi_{tail}$ hesapla, Eq. (6.6)
4. $\\phi_{dome} = \\phi(s_\\varepsilon) + \\Delta\\phi_{tail}$
5. Eşit $\\Delta s$ aralıklarla nokta dizisi oluştur: $\\{s_i, \\rho_i, x_i, \\phi_i, \\alpha_i\\}$

### 12.4 Silindir Geçişi (Ekvator-1 → Ekvator-2)

1. $\\alpha_{cyl} = \\alpha_{eq}$ (sabit)
2. $\\phi$ doğrusal artar: $\\phi(x) = \\phi_{dome} + x \\cdot \\tan\\alpha_{cyl} / R_{eq}$
3. $N_{cyl}$ adet eşit aralıklı nokta üret

### 12.5 Dome-2 Yolu (Ekvator → Dönüş Noktası)

Karar-2 simetri: dome-1 ile aynı $\\phi_{dome}$ katkısı. Yay uzunluğu global
koordinatta devam eder. $\\phi$ birikimli olarak artar.

### 12.6 Return Stroke (Dome-2 Dönüş → Dome-1 Dönüş)

Return stroke, forward stroke'un ayna görüntüsüdür:
- Dome-2 → ekvator-2: $\\phi$ artar ($+\\phi_{dome}$)
- Ekvator-2 → ekvator-1: $\\phi$ artar ($+\\phi_{cyl}$)
- Ekvator-1 → dome-1 dönüş: $\\phi$ artar ($+\\phi_{dome}$)

### 12.7 Sonuç

$$\\Delta\\phi_{circuit} = 4 \\cdot \\phi_{dome} + 2 \\cdot \\phi_{cyl} \\tag{12.1}$$

### 12.8 Adaptif ODE Çıktısından Uniform Izgara Yeniden Örnekleme

Boost.Odeint adaptif adım çözücüsü, dome dönüş noktası yakınında adımları
otomatik küçültür, silindirik bölgede ise büyütür. Bu durum düzensiz aralıklı
$\\{s_j, \\phi_j\\}$ çıktısı üretir. Ancak Phase-3 kinematik çözücüsü ve Phase-4
trajectory planner, uniform $\\Delta s$ adımlarında örneklenmiş yol verisi bekler
(Koussios Bölüm 10.2, eşit yay uzunluğu artışları prensibi).

**Yeniden örnekleme stratejisi:**

**Adım 1 — Ham ODE çıktısı:** Adaptif çözücü $\\{s_j, \\phi_j\\}_{j=0}^{M}$
üretir; $M$ ve $s_j$ aralıkları çözücünün iç kararıdır.

**Adım 2 — Monoton kübik interpolasyon:** Ham $\\phi(s)$ verileri üzerine
`pchip` (piecewise cubic Hermite interpolating polynomial) uygulayarak
$\\tilde{\\phi}(s)$ sürekli fonksiyonu oluşturulur. `pchip` tercih sebebi:

- Monotonluk koruması: $\\phi(s)$ dome bölgesinde kesinlikle monoton artar;
  standart kübik spline bu garantiyi vermez ve aşma (overshoot) üretebilir
- $C^1$ süreklilik: Phase-3 kinematik denklemlerinin türev sürekliliği
  gereksinimiyle uyumlu
- Koussios'un Bölüm 10.1'deki "interpolating polynomial" önerisinin
  pratik karşılığı

**Adım 3 — Uniform ızgara:** Ana meridyen profilinin uniform yay uzunluğu
noktaları $\\{s_i\\}_{i=0}^{N}$ ($\\Delta s = s_{total} / N$, $N$ = `N_points`
parametresi) zaten MeridianLookupTable'da mevcuttur. $\\tilde{\\phi}(s)$
interpolasyonu bu noktalarda değerlendirilir:

$$\\phi_i = \\tilde{\\phi}(s_i), \\quad i = 0, 1, \\ldots, N \\tag{12.2}$$

**Adım 4 — Türetilen büyüklükler:** Uniform $s_i$ noktalarında $\\rho_i$,
$x_i$ doğrudan MeridianLookupTable'dan sorgulanır; $\\alpha_i$ Clairaut'tan
(Eq. 4.2) hesaplanır. Nihai çıktı $\\{s_i, \\rho_i, x_i, \\phi_i, \\alpha_i\\}$
tamamen uniform aralıklıdır.

**Silindirik bölge istisnası:** Silindirik bölgede $\\phi(x)$ zaten doğrusaldır
(Eq. 7.2), yeniden örnekleme gerekmez — doğrudan analitik formülden
hesaplanır.

**Doğrulama:** Yeniden örnekleme sonrası Clairaut invariant kontrolü
(Eq. 14.1) tekrar uygulanır. Interpolasyon hatası Karar-11 Katman 2
toleransı ($10^{-4}$ mm) içinde kalmalıdır.

---

## 13. Hata Yönetimi (Karar-17)

### 13.1 Alt Katman (Geometri Sorgu)

`std::optional` dönüş tipi, exception yok.

| Durum | Dönüş |
|-------|-------|
| $\\rho(s) < R_E$ (arcsin argümanı > 1) | `std::nullopt` |
| MeridianLookupTable aralık dışı sorgu | `std::nullopt` |
| $\|\\dot{\\rho}_{turn}\| < 10^{-8}$ (dejenere dönüş) | `std::nullopt` |
| Ters sorgu yakınsamama | `std::nullopt` |
| $k_n \\leq 0$ tespit (S-WIND-04) | Değer döner + `bridging_risk` flag set |

### 13.2 Üst Katman (Kullanıcı Arayüzü)

Exception fırlatma.

| Durum | Exception |
|-------|-----------|
| Hoop winding + dome yolu talebi | `std::invalid_argument` |
| $R_E \\geq R_{eq}$ (fiziksel olarak imkansız) | `std::invalid_argument` |
| $BW_{eff} \\leq 0$ veya $r_0 \\leq 0$ | `std::invalid_argument` |
| Dome entegrasyonu başarısız | `std::runtime_error` |
| $k_n \\leq 0$ tespiti (S-WIND-04) | Uyarı (warning), exception değil |
| $\\delta c / (\\rho \\cdot k_n) > \\mu_{max}$ (S-WIND-05) | Uyarı (warning), exception değil |

---

## 14. Toleranslar ve Doğrulama (Karar-11)

### 14.1 ODE Entegrasyon Toleransı (Katman 1)

$$\\text{atol} = 10^{-10}, \\quad \\text{rtol} = 10^{-8}$$

### 14.2 Clairaut Doğrulama Toleransı (Katman 2)

Her üretilen nokta için:

$$|\\rho_i \\sin\\alpha_i - R_E| < 10^{-4} \\; \\text{mm} \\tag{14.1}$$

Bu, geodesic koşulunun karşılandığının çalışma-zamanı doğrulamasıdır.

### 14.3 Silindir-Dome Geçiş Toleransı (Katman 3)

Geçiş noktasında winding açısı süreksizliği kontrolü:

$$|\\alpha_{dome}(s = 0) - \\alpha_{cyl}| < 10^{-6} \\; \\text{rad} \\tag{14.2}$$

Bu, $\\rho_{dome}(s=0) = R_{eq}$ ve Clairaut'un her iki bölgede de aynı sabiti
kullanmasından otomatik olarak sağlanır — numerik doğrulama amaçlıdır.

---

## 15. MATLAB Referans İmplementasyon Yapısı

### 15.1 Dosya Organizasyonu

```
MATLAB/phase2a_winding/
├── geodesic_single_circuit.m     % Ana fonksiyon — tek devre yolu üretimi
├── verify_geodesic_path.m        % Doğrulama scripti (GATE-2a)
└── utils/
    └── dome_phi_integration.m    % Dome φ entegrasyonu (adaptif kuadratür)
```

### 15.2 geodesic_single_circuit.m Arayüzü

```
Girdi:
  dome_profile   — Phase-1a profil struct (hemispherical/ellipsoidal/isotensoid)
  R_eq           — ekvator yarıçapı [mm]
  r0             — polar açıklık yarıçapı [mm]
  L_cyl          — silindir uzunluğu [mm]
  BW_eff         — efektif bant genişliği [mm]
  N_points       — yol nokta sayısı (varsayılan: 2000)

Çıktı struct:
  .s             — global yay uzunluğu [mm] (N_points × 1)
  .rho           — yarıçap [mm]
  .x             — aksiyel konum [mm]
  .phi           — paralel açı [rad]
  .alpha          — winding açısı [rad]
  .kn            — normal eğrilik [1/mm] (S-WIND-04 bridging kontrolü)
  .lambda_max    — maksimum kayma eğilimi oranı δc/(ρ·kn) (S-WIND-05)
  .delta_phi_circuit — tek devre toplam açısal ilerleme [rad]
  .phi_dome      — tek dome yarısı φ katkısı [rad]
  .phi_cyl       — tek silindir geçişi φ katkısı [rad]
  .bridging_risk_indices — kn ≤ 0 olan nokta indeksleri (boş = risk yok)
```

### 15.3 Doğrulama Testleri (GATE-2a)

| Test | Açıklama | Kriter |
|------|----------|--------|
| T1 | Clairaut invariant: $\\rho \\sin\\alpha = R_E$ her noktada | max sapma $< 10^{-4}$ mm |
| T2 | Dome dönüş noktasında $\\alpha \\to \\pi/2$ | $|\\alpha_{turn} - \\pi/2| < 10^{-3}$ rad |
| T3 | Silindir-dome geçiş sürekliliği | $|\\Delta\\alpha| < 10^{-6}$ rad |
| T4 | Hemispherical dome: $\\phi_{dome} \\approx \\pi/2$ (küre özelliği) | %1 tolerans |
| T5 | Epsilon yakınsama: $\\phi_{dome}(\\varepsilon)$ kararlılığı | 4 farklı $\\varepsilon$ |
| T6 | Endüstri test senaryoları (TEST-01..04) ile çapraz doğrulama | Karar-16 |
| T7 | Yol 3D sürekliliği: ardışık noktalar arası mesafe kontrolü | Sıçrama yok |
| T8 | Hoop winding kısıt kontrolü | Dome yolu üretilmemeli |
| T9 | Konveksite kontrolü (S-WIND-04): tüm yol boyunca $k_n > 0$ | Bridging bölgesi yok |
| T10 | Kayma emniyet marjı (S-WIND-05): $\\delta c / (\\rho \\cdot k_n) < \\mu_{max}$ | Winding_type'a göre $\\mu$ |
| T11 | Uniform yeniden örnekleme sonrası Clairaut invariant | max sapma $< 10^{-4}$ mm |

---

## 16. Phase-1b API Kullanım Haritası

Phase-2a geodesic solver, Phase-1b'nin aşağıdaki API'lerini kullanır:

| Phase-1b Bileşeni | Kullanım |
|--------------------|----------|
| `MandrelGeometry::point(s_global)` | $\\rho(s)$, $x(s)$ sorgusu |
| `MandrelGeometry::isOnDome1/Cylinder/Dome2(s)` | Bölge tespiti |
| `MeridianLookupTable::query(s)` | $\\rho$, $d\\rho/ds$, $dx/ds$, $\\kappa_m$ |
| `MeridianLookupTable::inverseLookup(\\rho)` | $s_\\varepsilon$ tespiti ($\\rho \\to s$) |
| `MandrelGeometry::sTotal()` | Toplam yay uzunluğu |
| `MandrelGeometry::sDome()` | Dome yay uzunluğu |
| `MandrelGeometry::cylinderLength()` | $L_{cyl}$ |
| `MandrelGeometry::equatorRadius()` | $R_{eq}$ |

**Türetilen büyüklük — paralel eğrilik:** $\\kappa_p$ doğrudan MeridianLookupTable'da
depolanmaz ancak mevcut verilerden türetilebilir (Eq. 9.3). Phase-2a solver'ı
$\\kappa_p(s) = -\\dot{\\rho}(s) / (\\rho(s) \\cdot \\sqrt{1 - \\dot{\\rho}(s)^2})$ formülüyle
her noktada hesaplar. Bu hesap S-WIND-04 bridging kontrolü (Eq. 9.2) için
zorunludur.

**Kritik not:** `inverseLookup(rho)` fonksiyonu Phase-1b'de henüz implement
edilmemiş olabilir. Bu, Phase-2a implementasyonunun ilk adımlarından biri
olarak MeridianLookupTable'a eklenmelidir.

**İmplementasyon detayları:**

Dome bölgesinde $\\rho(s)$ monoton azalır ($d\\rho/ds < 0$), dolayısıyla
$\\rho \\to s$ ters eşleme tekil değerlidir ve binary search ile $O(\\log n)$
karmaşıklıkta çözülebilir. Ancak aşağıdaki sınır durumları özel işlem
gerektirir:

(a) **Ekvator noktası** ($\\rho = R_{eq}$, $s = 0$): Bu noktada $d\\rho/ds = 0$
(teğet yatay, meridyen eğim açısı $\\beta = 0$). Silindir-dome geçişinde
$\\rho(s)$'ın türevi sıfır olduğundan, binary search bu bölgede yavaş
yakınsayabilir. Çözüm: $|\\rho - R_{eq}| < $ Karar-11 Katman 3 toleransı
($10^{-6}$ mm) durumunda doğrudan $s = 0$ döndürülmeli.

(b) **Polar açıklık noktası** ($\\rho = R_E$, $s \\approx s_{dome}$): Dönüş
noktası yakınında $d\\rho/ds$ çok küçük olabilir (isotensoid'de S-GEO-01
stiff davranış). $|\\rho - R_E| < \\varepsilon$ durumunda $s_\\varepsilon$ doğrudan
döndürülmeli — bu, Bölüm 6.2'deki analitik tamamlama stratejisiyle
tutarlıdır.

(c) **Aralık dışı sorgular:** $\\rho > R_{eq}$ veya $\\rho < r_0$ durumunda
`std::nullopt` döndürülmeli (Karar-17 alt katman kuralı).

(d) **Silindirik bölge:** $\\rho = R_{eq}$ sabit olduğundan ters sorgu
anlamsızdır. `inverseLookup` yalnızca dome bölgesinde çağrılmalıdır;
çağıran kodun bölge kontrolü yapması gerekir (`isOnDome1/Dome2`).

**Yakınsama kriteri:** Binary search, $|\\rho(s_{mid}) - \\rho_{target}| < 10^{-6}$
mm (Karar-11 Katman 3) toleransına ulaşınca sonlandırılır. Maksimum
iterasyon sayısı $\\lceil \\log_2(N_{table}) \\rceil + 5$ ile sınırlandırılır.

---

## 17. Referanslar

- Koussios, S. — "Filament Winding: A Unified Approach" — Bölüm 2, 4, 5, 6, 8, 9, 10
- Peters, S.T. (ed.) — "Composite Filament Winding" — Bölüm 4, 5, 9
- Vargas-Rojas & Collombet — "Non-geodesic filament winding" (CAD 181, 2025) — Appendix B
- Do Carmo — "Differential Geometry of Curves and Surfaces" — Bölüm 4
- Phase-0 karar kaydı: docs/phase0_decisions.md
- Phase-1b GATE raporu: docs/gate_1b_01_report.md
