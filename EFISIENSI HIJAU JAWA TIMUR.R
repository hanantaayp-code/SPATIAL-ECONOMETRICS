# ==============================================================================
# ------------------------- INFORMASI PROYEK & PENELITIAN ----------------------
# ==============================================================================
# JUDUL       : Limpahan Spasial Efisiensi Ekonomi Hijau di Jawa Timur, Indonesia: 
#               Peran Kualitas Kelembagaan dalam Memoderasi Tekanan Sektoral
# TARGET      : Publikasi Jurnal Wilayah dan Lingkungan (SINTA 2)
# TIM PENULIS : 1. Hananta Angger
#               2. Rinto Martinus
#               3. Himawan Nugraha
# PEMBIMBING  : Dr. agr. Deden Dinar Iskandar
# DESKRIPSI   : Script replikasi untuk estimasi model Data Panel (Fixed/Random Effect)
#               dan Ekonometrika Spasial (SAR/SEM/SDM) serta pemetaan LISA.
# INPUT       : 1. Shapefile Jawa Timur (.shp) peta wilayah
#               2. Data Excel (.xlsx) berisi variabel ekonomi/sosial per tahun
# OUTPUT      : 1. Model statistik
#               2. Uji autokorelasi spasial (MORAN I), Peta Kluster Local Autokorelasi (LISA)
#               2. Estimasi Panel Non Spasial
#               4. Estimasi Panel Spasial Fixed Effect Model "Individual"
#               5. Uji Diasnostik dan Goodnes of fit
# NOTE        : Script by Rinto Martinus, edited by Nanta on December 14, 2025
# ==============================================================================

# ==============================================================================
# --- SPESIFIKASI MODEL DASAR ---
# ==============================================================================
# SBM_Efficiency ~ SAKIP + AGRN + SXA + LMS + POP + URB + INT + LNINF
# VARIABEL DEPENDEN  : NILAI EFISIENSI HIJAU HASIL OLAH MENGGUNAKAN SBM DEA
# VARIABEL INDEPENDEN: 1. SAKIP : PROKSI KUALITAS INSTITUSI
#                      2. AGRN  : PROKSI KEGIATAN EKONOMI SEKTOR AGRIKULTUR
#                      3. SXA
#                      4. LMS   : RATA-RATA LAMA SEKOLAH
#                      5. POP   : POPULASI
#                      6. URB   : PRESENTASE PENDUDUK KOTA
#                      7. INT   : 
#                      8. LNINF : INFRASTRUKTUR DIPROKSI DENGAN PANJANG JALAN MANTAP
# =================================================================================

# ==============================================================================
# --- 1. INSTALL & LOAD LIBRARY ---
# ==============================================================================
# Jika belum install, jalankan:
# install.packages(c("sf", "spdep", "splm", "dplyr", "plm", 
#                    "readxl", "tmap", "ggplot2", "lmtest"))

library(sf)       # Manipulasi data spasial (peta)
library(spdep)    # Matriks pembobot & Uji Moran
library(splm)     # Model panel spasial (SAR/SEM/SDM)
library(dplyr)    # Data wrangling (filter, join, dll)
library(plm)      # Model panel standar (OLS/FE/RE)
library(readxl)   # Baca Excel
library(tmap)     # Visualisasi Peta
library(ggplot2)  # Visualisasi Grafik
library(lmtest)   # Uji diagnostik tambahan, Menguji jenis Autokorelasi ( Y, X, OR EROR)

# ==============================================================================
# BAGIAN I: PERSIAPAN DATA (DATA PREPARATION)
# ==============================================================================

# --- 2. Load Data ---
# [PENTING] Sesuaikan path di bawah ini dengan lokasi file di komputer Anda
path_shp  <- "D:/Lomba/Jawa Timur/Jawa_Timur_ADMIN_BPS.shp" 
path_data <- "DATA FIX.xlsx"

shapefile <- st_read(path_shp, quiet = TRUE)
data_FIX  <- read_excel(path_data)

# --- 3. Cleaning & Standardisasi Nama Wilayah ---
# Mengubah nama kabupaten menjadi huruf kecil dan menghapus spasi berlebih
shapefile$Kabupaten <- trimws(tolower(shapefile$Kabupaten))
data_FIX$Kabupaten  <- trimws(tolower(data_FIX$Kabupaten))

# --- 4. Membangun Struktur Data Panel pada Shapefile ---
# Shapefile asli bersifat cross-section (1 tahun). Kita perlu duplikasi baris
# agar sesuai dengan jumlah tahun pengamatan (5 Tahun: 2018-2022).
n_tahun <- 5
periode <- 2018:2022

shp_panel <- shapefile[rep(1:nrow(shapefile), each = n_tahun), ]
shp_panel$Tahun <- rep(periode, times = nrow(shapefile))

# --- 5. Menggabungkan (Merge) Peta dengan Data Excel ---
shp_merged <- left_join(shp_panel, data_FIX, by = c("Kabupaten", "Tahun"))

# --- 6. Hapus Duplikat & Validasi ---
shp_merged <- shp_merged %>% 
distinct(Kabupaten, Tahun, .keep_all = TRUE)

# Buat dataframe non-spasial untuk analisis panel standar
panel_df <- st_drop_geometry(shp_merged)

# Hapus kolom yang varians-nya 0 (jika ada)
varians  <- sapply(panel_df[, sapply(panel_df, is.numeric)], var, na.rm=TRUE)
panel_df <- panel_df[, c("Kabupaten", "Tahun", names(varians[varians > 0]))]

# Set index panel (Individu = Kabupaten, Waktu = Tahun)
pdata <- pdata.frame(panel_df, index = c("Kabupaten", "Tahun"))

# --- 7. Definisi Model (FORMULA) ---
# Y = SBM_Efficiency
# X = SAKIP, AGRN, SXA, LMS, POP, URB, INT, LNINF
model_formula <- SBM_Efficiency ~ SAKIP + AGRN + SXA + LMS + POP + URB + INT + LNINF

# ==============================================================================
# BAGIAN II: ANALISIS DATA PANEL STANDAR (NON-SPASIAL)
# ==============================================================================
# Bagian ini untuk menentukan apakah model terbaik adalah OLS, FE, atau RE
# sebelum masuk ke pemodelan spasial.

# --- A. Estimasi Model Standar ---

# 1. Common Effect Model (CEM) / Pooled OLS
pooled_model <- plm(model_formula, data = pdata, model = "pooling")

# 2. Fixed Effect Model (FEM)
fe_model <- plm(model_formula, data = pdata, model = "within")

# 3. Random Effect Model (REM)
re_model <- plm(model_formula, data = pdata, model = "random")

# --- B. Uji Pemilihan Model (Chow & Hausman) ---

# 1. Uji Chow (Membandingkan CEM vs FEM)
# H0: CEM lebih baik; H1: FEM lebih baik
# Kita lakukan manual F-test karena struktur data spasial kadang kompleks
pool_lm <- lm(model_formula, data = panel_df)
fixed_lm <- lm(model_formula, data = panel_df, factor(Kabupaten)) # Dummy wilayah
chow_test <- anova(pool_lm, fixed_lm)

print("--- HASIL UJI CHOW (CEM vs FEM) ---")
print(chow_test) 
# Jika Pr(>F) < 0.05, maka Tolak H0 -> Pilih FEM.

# 2. Uji Hausman (Membandingkan FEM vs REM)
# H0: REM lebih baik (konsisten); H1: FEM lebih baik
hausman_test <- phtest(fe_model, re_model)

print("--- HASIL UJI HAUSMAN (FEM vs REM) ---")
print(hausman_test)
# Jika p-value < 0.05, maka Tolak H0 -> Pilih FEM (Fixed Effect).
# (Biasanya data level kabupaten/provinsi cenderung ke Fixed Effect).

# ==============================================================================
# BAGIAN III : MEMBUAT MATRIKS PEMBOBOT SPASIAL (SPATIAL WEIGHTS)
# ==============================================================================

# Kita butuh peta referensi (hanya 1 tahun) untuk membuat hubungan antar wilayah
shp_ref <- shp_merged %>% filter(Tahun == 2018) %>% st_make_valid()

# --- Opsi 1: Pembobotan Queen (Berdasarkan persinggungan batas) ---
nb_queen    <- poly2nb(shp_ref)
listw_queen <- nb2listw(nb_queen, style = "W")
# Cek hasil Queen
# summary(listw_queen)
# moran.test(as.numeric(shp_ref$SBM_Efficiency), listw_queen)

# --- Opsi 2: Pembobotan Kernel (Berdasarkan jarak Gaussian) ---
coords   <- st_coordinates(st_centroid(shp_ref))
dist_mat <- as.matrix(dist(coords))
bandwidth <- mean(dist_mat) # Bandwidth rata-rata jarak

gaussian_kernel <- function(d, bw) { exp(-(d^2) / (2 * bw^2)) }
W_kernel <- gaussian_kernel(dist_mat, bandwidth)
diag(W_kernel) <- 0 # Wilayah tidak mempengaruhi dirinya sendiri
W_kernel_norm  <- W_kernel / rowSums(W_kernel) # Normalisasi baris
W_kernel_norm[is.na(W_kernel_norm)] <- 0
listw_kernel <- mat2listw(W_kernel_norm, style = "W")

# --- Opsi 3: Pembobotan KNN (k-Nearest Neighbors) ---
# (Ini yang akan digunakan utama dalam model di bawah)
k_neighbors <- 4 # Jumlah tetangga terdekat
nb_knn      <- knn2nb(knearneigh(coords, k = k_neighbors))
listw_knn   <- nb2listw(nb_knn, style = "W")


# ==============================================================================
# BAGIAN IV: UJI AUTOKORELASI SPASIAL GLOBAL MENGGUNAKAN MORAN I
# ==============================================================================
# Uji Global Moran's I dengan QUEEN
moran_queen_result <- moran.test(as.numeric(shp_ref$SBM_Efficiency), listw_queen)
print("--- Hasil Global Moran's I (QUEEN) ---")
print(moran_queen_result)

# Uji Global Moran's I dengan KERNEL
moran_kernel_result <- moran.test(as.numeric(shp_ref$SBM_Efficiency), listw_kernel)
print("--- Hasil Global Moran's I (KERNEL) ---")
print(moran_kernel_result)

# Uji Global Moran's I dengan KNN
moran_knn_result <- moran.test(as.numeric(shp_ref$SBM_Efficiency), listw_knn)
print("--- Hasil Global Moran's I (KNN) ---")
print(moran_knn_result)

# ==============================================================================
# BAGIAN V: ANALISIS LOKAL (LISA / Local Moran's I) & MAPPING
# ==============================================================================
# Loop otomatis untuk tahun 2018-2022 agar script tidak panjang & repetitif

tmap_mode("view") # Set mode peta interaktif

# Fungsi untuk menentukan Kuadran LISA
get_quadrant <- function(val, lag, p_val) {
  mean_val <- mean(val, na.rm=TRUE)
  mean_lag <- mean(lag, na.rm=TRUE)
  
  quad <- ifelse(p_val > 0.15, "Not significant",
          ifelse(val >= mean_val & lag >= mean_lag, "High-High",
          ifelse(val <= mean_val & lag <= mean_lag, "Low-Low",
          ifelse(val >= mean_val & lag <= mean_lag, "High-Low",
          ifelse(val <= mean_val & lag >= mean_lag, "Low-High", "Not significant")))))
  return(quad)
}

# Loop Proses LISA per Tahun
for (thn in 2018:2022) {
  
  cat(paste0("\n--- Memproses LISA Tahun ", thn, " ---\n"))
  
  # 1. Filter Data & Validasi Geometri
  shp_curr <- shp_merged %>% filter(Tahun == thn) %>% st_make_valid()
  
  # 2. Ambil Nilai Variabel & Standarisasi
  nilai_curr <- as.numeric(scale(shp_curr$SBM_Efficiency))
  
  # 3. Hitung Local Moran
  lm_res <- localmoran(nilai_curr, listw_knn)
  
  # 4. Masukkan Hasil ke Shapefile
  shp_curr$Ii     <- lm_res[, 1]
  shp_curr$P.Ii   <- lm_res[, 5]
  
  # 5. Hitung Lag Spatial & Tentukan Kuadran
  lag_curr <- lag.listw(listw_knn, nilai_curr)
  shp_curr$quad <- get_quadrant(nilai_curr, lag_curr, shp_curr$P.Ii)
  
  # 6. Plot Peta LISA (Tmap)
  peta <- tm_shape(shp_curr) +
    tm_fill("quad",
            palette = c("High-High" = "red", "Low-Low" = "blue", 
                        "High-Low" = "yellow", "Low-High" = "green", 
                        "Not significant" = "grey90"),
            title = paste("LISA Quadrant", thn)) +
    tm_borders(alpha = 0.5) +
    tm_layout(legend.outside = TRUE)
  
  print(peta) # Tampilkan peta
  
  # 7. Plot Scatter Moran (Ggplot2)
  df_scatter <- data.frame(Nilai = nilai_curr, Lag = lag_curr, Quad = shp_curr$quad)
  
  plot_scatter <- ggplot(df_scatter, aes(x = Nilai, y = Lag, color = Quad)) +
    geom_point() +
    geom_smooth(method = "lm", se = FALSE, color = "black", linetype="dashed") +
    geom_hline(yintercept = 0) + geom_vline(xintercept = 0) +
    labs(title = paste("Moran Scatter Plot -", thn),
         x = "SBM Efficiency (Standardized)", y = "Spatial Lag") +
    scale_color_manual(values = c("High-High" = "red", "Low-Low" = "blue", 
                                  "High-Low" = "yellow", "Low-High" = "green", 
                                  "Not significant" = "grey")) +
    theme_minimal()
  
  print(plot_scatter) # Tampilkan scatter plot
}


# ==============================================================================
# BAGIAN VII: UJI DIAGNOSTIK SPASIAL (LAGRANGE MULTIPLIER TEST)
# ==============================================================================
# Tujuannya: Menentukan apakah spasial efek ada di Lag (SAR) atau Error (SEM).
# Hipotesis Nol (H0): Tidak ada autokorelasi spasial.
# Jika p-value < 0.05, maka tolak H0 (artinya ada pengaruh spasial).

cat("\n--- MELAKUKAN UJI LAGRANGE MULTIPLIER (LM TEST) ---\n")

# 1. LM Test untuk Spatial Lag (SAR)
lm_lag <- slmtest(model_formula, data = pdata, listw = listw_knn, 
                  model = "within", test = "lml")
print("--- 1. LM Test for Spatial Lag (SAR) ---")
print(lm_lag)

# 2. LM Test untuk Spatial Error (SEM)
lm_err <- slmtest(model_formula, data = pdata, listw = listw_knn, 
                  model = "within", test = "lme")
print("--- 2. LM Test for Spatial Error (SEM) ---")
print(lm_err)

# --- UJI ROBUST (Dilakukan jika KEDUANYA Signifikan) ---
# Jika LM Lag DAN LM Error sama-sama signifikan (p < 0.05),
# maka kita lihat versi Robust-nya (RLM).

# 3. Robust LM Test untuk Spatial Lag
rlm_lag <- slmtest(model_formula, data = pdata, listw = listw_knn, 
                   model = "within", test = "rlml")
print("--- 3. Robust LM Test for Spatial Lag ---")
print(rlm_lag)

# 4. Robust LM Test untuk Spatial Error
rlm_err <- slmtest(model_formula, data = pdata, listw = listw_knn, 
                   model = "within", test = "rlme")
print("--- 4. Robust LM Test for Spatial Error ---")
print(rlm_err)


# ==============================================================================
# PANDUAN KEPUTUSAN (RULE OF THUMB) - ANSELIN (1988)
# ==============================================================================
# Catat hasil p-value di atas, lalu ikuti aturan ini untuk Publikasi:
#
# KONDISI A: Hanya LM Lag yang signifikan -> Gunakan Model SAR.
# KONDISI B: Hanya LM Error yang signifikan -> Gunakan Model SEM.
# KONDISI C: Keduanya signifikan -> Lihat Robust LM (RLM).
#            - Jika RLM Lag sig & RLM Error tidak -> SAR.
#            - Jika RLM Error sig & RLM Lag tidak -> SEM.
#            - Jika keduanya masih signifikan -> SDM (Spatial Durbin Model) 
#              seringkali menjadi pilihan teraman karena mencakup keduanya.
# ==============================================================================

# ==============================================================================
# BAGIAN VII: MODEL EKONOMETRIKA SPASIAL (SAR, SEM, SDM)
# ==============================================================================
# Menggunakan pembobot KNN (listw_knn) dan asumsi Fixed Effect ("individual")

# 1. Spatial Autoregressive Model (SAR)
# Lag pada variabel dependen (Y dipengaruhi Y tetangga)
model_sar <- spml(
  formula = model_formula,
  data = pdata,
  listw = listw_knn,
  spatial.error = "none", 
  lag = TRUE, 
  model = "within", # Menggunakan Fixed Effect (sesuaikan jika hasil Hausman beda)
  effect = "individual"
)
print("--- Ringkasan Model SAR ---")
summary(model_sar)

# 2. Spatial Error Model (SEM)
# Error saling berkorelasi antar wilayah
model_sem <- spml(
  formula = model_formula,
  data = pdata,
  listw = listw_knn,
  spatial.error = "b", 
  lag = FALSE, 
  model = "within",
  effect = "individual"
)
print("--- Ringkasan Model SEM ---")
summary(model_sem)

# 3. Spatial Durbin Model (SDM)
# Y dipengaruhi Y tetangga DAN X tetangga
model_sdm <- spml(
  formula = model_formula,
  data = pdata,
  listw = listw_knn,
  spatial.error = "none", 
  lag = TRUE, 
  spatial.lag.x = TRUE, # Lag pada variabel independen
  model = "within",
  effect = "individual"
)
print("--- Ringkasan Model SDM ---")
summary(model_sdm)

# =================================================================================
# BAGIAN VIII: GOODNES OF FIT ---
# =================================================================================

# ---. Pemilihan Model Terbaik (AIC Comparison) ---

# Fungsi hitung AIC manual (karena spml tidak output AIC langsung)
# AIC = -2 * LogLik + 2 * (Jumlah Parameter)
get_aic <- function(model, k_extra=0) {
  ll <- model$logLik
  k  <- length(model$coefficients) + k_extra
  return(-2 * ll + 2 * k)
}

aic_vals <- data.frame(
  Model = c("SAR", "SEM", "SDM"),
  LogLikelihood = c(model_sar$logLik, model_sem$logLik, model_sdm$logLik),
  AIC = c(get_aic(model_sar, 1), get_aic(model_sem, 1), get_aic(model_sdm, 1))
)

print("--- PERBANDINGAN MODEL (Pilih AIC Terkecil) ---")
print(aic_vals)

# ==============================================================================
# BAGIAN IX: DEKOMPOSISI DAMPAK (DIRECT, INDIRECT, TOTAL) - WAJIB UNTUK SDM
# ==============================================================================
# Karena model spasial (terutama SDM) mengandung feedback loop, koefisien
# tidak bisa diinterpretasi langsung. Harus dihitung marginal effect-nya.
# Bagian ini akan menghitung rata-rata dampak.

cat("\n--- MENGHITUNG DAMPAK LANGSUNG & SPILLOVER (KHUSUS SDM) ---\n")

# Cek apakah SDM model terbaik? Jika ya, jalankan kode ini:
if(aic_vals$AIC[3] <= aic_vals$AIC[1] & aic_vals$AIC[3] <= aic_vals$AIC[2]) {
  
  # Hitung Impacts
  # R = 100 simulasi (bisa dinaikkan jadi 500-1000 untuk hasil final publikasi)
  impacts_sdm <- impacts(model_sdm, listw = listw_knn, R = 100, time = n_tahun) 
  
  print("--- Direct, Indirect (Spillover), & Total Effects ---")
  print(summary(impacts_sdm, zstats=TRUE, short=TRUE))
  
  cat("NOTE: 
  - Direct: Efek perubahan X di suatu kab terhadap Y di kab itu sendiri (termasuk feedback).
  - Indirect: Efek perubahan X di suatu kab terhadap Y di kab TETANGGA (Spillover).
  - Total: Jumlah keduanya.\n")
  
} else {
  print("SDM bukan model dengan AIC terendah. Silakan sesuaikan kode impacts untuk SAR jika SAR terpilih.")
}

#====================================================================================
print("---PUJI TUHAN OLAH DATA SELESAI ---")
print("--- SCRIPT SELESAI: SIAP UNTUK PUBLIKASI ---")
#===================================================================================

