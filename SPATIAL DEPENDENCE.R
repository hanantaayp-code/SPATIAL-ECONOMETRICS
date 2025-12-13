#===============================================================================
#   INTTERREGIONAL SPILLOVER EFFCT AND LOCAL DETERMINAN OF FOOD SECURITY IN EASTERN EAST INDONESIA
#   SKRIPSI HANANTA ANGGER YUGA PRAWIRA 
#   PART #SPATIAL PANEL DATA ANALYSIS (TWO-WAYS FIXED EFFECTS)
#
#   TUJUAN          : Estimasi Model Spasial (SAR, SEM, SDM) dengan Twoways Effect
#   DATA            : Y = IKP (Level), X = Transformed (Log Level, Log(X+1))
#   METODE ESTIMASI : Spatial Panel Fixed Effects (Twoways)
#   MATRIKS W       : K-Nearest Neighbors (KNN) Optimasi Moran's I
#   
#   LOG PERUBAHAN   :
#   1. STATISTIKA DESKRIPTIF
#   2. PANEL PEROSEDUR
#   3. AUTOKORELASI SPASIAL (GLOBAL AND LOCAL)
#   4. LM TEST AND ROBAS LM TEST
#   5. SPASIAL TWO WAYS EFFECT
#
#===============================================================================

# 0.0. BERSIH-BERSIH MEMORY
rm(list = ls()) 
gc() 
cat("\014") # Clear console

#===============================================================================
# 0.1. LOAD PUSTAKA (PACKAGES)
#===============================================================================
message("--- TAHAP 0.1: Memuat Library ---")

# --- GRUP 1: SPASIAL & PETA ---
library(sf)             # Manipulasi data vektor spasial
library(tmap)           # Visualisasi peta tematik
library(classInt)       # Klasifikasi interval kelas peta
library(spdep)          # Dependensi spasial (Moran's I, Matriks W)
library(splm)           # Spatial Panel Linear Models (SAR, SEM, SDM)

# --- GRUP 2: DATA WRANGLING ---
library(tidyverse)      # Koleksi paket data science (dplyr, ggplot2, tidyr)
library(janitor)        # Membersihkan nama variabel
library(readxl)         # Baca Excel
library(writexl)        # Tulis Excel
library(fuzzyjoin)      # Join data yang namanya agak beda dikit

# --- GRUP 3: EKONOMETRIKA DASAR ---
library(plm)            # Panel Data Linear Models (Standard)
library(lmtest)         # Uji diagnostik (BP test, dll)
library(car)            # VIF test
library(sandwich)       # Robust standard errors
library(broom)          # Merapikan output model

# --- GRUP 4: VISUALISASI ---
library(corrplot)       # Plot korelasi matriks
library(GGally)         # Scatterplot matrix

#===============================================================================
# 0.2. FUNGSI PEMBANTU (HELPER FUNCTIONS)
#===============================================================================
message("--- TAHAP 0.2: Definisi Fungsi Helper ---")

# Normalisasi teks (uppercase + hapus spasi berlebih) untuk key join
normalize <- function(x) {
  x |> toupper() |> stringr::str_squish()
}

# Memperbaiki geometri yang error/invalid
fix_geom <- function(g) {
  sf::sf_use_s2(FALSE) # Matikan spherical geometry sementara
  if ("st_make_valid" %in% getNamespaceExports("sf")) {
    g_fixed <- sf::st_make_valid(g)
  } else {
    g_fixed <- sf::st_buffer(g, 0)
  }
  sf::sf_use_s2(TRUE)
  return(g_fixed)
}

# Membuat Matriks W berbasis KNN dan output diagnostik
build_knn_w <- function(coords, k) {
  nb  <- knn2nb(knearneigh(coords, k = k))
  lw  <- nb2listw(nb, style = "W", zero.policy = TRUE)
  deg <- sapply(nb, length)
  
  tib <- tibble(
    k = k,
    components = n.comp.nb(nb)$nc, # Cek apakah ada pulau terpisah
    min_deg = min(deg),
    mean_deg = mean(deg),
    max_deg = max(deg)
  )
  list(nb = nb, lw = lw, diag = tib)
}

#===============================================================================
# 0.3. SETTING DIREKTORI KERJA
#===============================================================================
# Ganti path ini sesuai lokasi folder skripsi di laptopmu
setwd("~/SKRIPSI") 
message("Direktori kerja saat ini: ", getwd())

#===============================================================================
# 1. DATA PREPARATION & CLEANING
#===============================================================================
message("--- TAHAP 1: Persiapan Data Panel & SHP ---")

# --- 1.A. Load & Fix Shapefile (Peta) ---
dir.create("data", showWarnings = FALSE)

# Cek apakah file gpkg sudah ada, jika belum buat dari SHP mentah
if (!file.exists("data/petaindotimur.gpkg")) {
    # Asumsi file SHP mentah ada di folder root
    shp_full <- st_read("INDONESIA.shp", quiet=TRUE) |> clean_names()
    daftar185 <- read_csv("NAMAKABKOT2.csv", show_col_types=FALSE) |> clean_names()
    
    # Matching nama
    shp_full <- shp_full |> mutate(nm_shp = normalize(kab_kota))
    daftar185 <- daftar185 |> mutate(nm_list = normalize(kabkot))
    
    matched_exact <- daftar185 |>
      inner_join(shp_full |> st_drop_geometry(), by = c("nm_list" = "nm_shp"))
    
    shp_185 <- shp_full |>
      semi_join(matched_exact, by = c("nm_shp" = "nm_list"))
    
    # Simpan hasil bersih
    sf::st_write(shp_185, "data/petaindotimur.gpkg", delete_dsn = TRUE, quiet = TRUE)
    message("File GPKG baru berhasil dibuat.")
}

# Load Peta Bersih
shp_185 <- st_read("data/petaindotimur.gpkg", quiet = TRUE) |> 
  clean_names() |> 
  mutate(nm_shp = normalize(kab_kota)) |>
  fix_geom() # Fix geometri

# --- 1.B. Load Excel Data Panel ---
xls_path  <- "DATA SET SKRIPSI.xlsx"
panel_raw <- read_xlsx(xls_path, sheet = 1) |> clean_names()

# Rename variabel agar coding lebih mudah
panel <- panel_raw |>
  select(-any_of("id")) |> 
  rename(
    tahun        = tahun,
    Kab_kota_raw = kab_kota, 
    ikp          = ikp,
    fiskal       = fiskal,
    bpk          = bpk,
    ctsr         = ctsr, 
    ln_poktan    = ln_poktan,
    ln_penyuluh  = ln_penyuluh,
    tpt          = tpt,
    akses        = akses,
    ln_kud       = ln_kud,
    rs           = rs,
    apm          = apm,
    rain         = rain,
    ln_lahan     = ln_lahan
  ) |>
  mutate(
    nm_list = normalize(kab_kota_data), # Gunakan kolom kab_kota_data sebagai key
    tahun   = as.integer(tahun)
  )

# --- 1.C. Join Peta + Data (Panel SF) ---
panel_sf <- shp_185 |>
  inner_join(panel, by = c("nm_shp" = "nm_list")) |>
  arrange(kab_kota, tahun) 

# Cek Validasi Baris (Harus 185 kab x 5 tahun = 925)
if (nrow(panel_sf) != 925) stop("ERROR: Jumlah baris data tidak 925. Cek proses join!")

# --- 1.D. Transformasi Variabel (Logaritma) ---
panel_df <- panel_sf |> st_drop_geometry()

panel_df_transformed <- panel_df %>%
  mutate(pop = as.numeric(pop)) %>%
  mutate(
    log_grdp    = log(grdp), 
    log_pop     = log(pop), 
    log_rice    = log1p(rice),    # log(x+1) untuk data ada nol
    log_nonrice = log1p(nonrice)
  ) %>%
  select(-Kab_kota_raw)

# Finalisasi Data Frame
panel_df <- panel_df_transformed |>
  mutate(id = as.factor(kab_kota)) |> 
  arrange(id, tahun) 

# Update Panel SF (Data Spasial + Atribut Baru)
panel_sf <- shp_185 %>%
  inner_join(panel_df, by = c("kab_kota", "nm_shp")) %>%
  arrange(kab_kota, tahun)

saveRDS(panel_sf, "data/panel_sf_eastindo.rds")
saveRDS(panel_df, "data/panel_df_eastindo.rds")
message("Data Panel berhasil disimpan (RDS).")

#===============================================================================
# 2. PEMBENTUKAN MATRIKS PEMBOBOT (SPATIAL WEIGHTS / W)
#===============================================================================
message("--- TAHAP 2: Seleksi Matriks W (KNN) ---")

# --- 2.A. Generate Kandidat W (K=4 sampai 12) ---
# Kita cari mana K yang memberikan Moran's I paling tinggi/stabil
sf::sf_use_s2(FALSE)
coords <- st_coordinates(st_point_on_surface(shp_185))
sf::sf_use_s2(TRUE)

K_list <- 4:12
moran_results <- list()

for (k in K_list) {
  # Buat W
  res <- build_knn_w(coords, k)
  
  # Uji Moran Global per tahun untuk K ini
  mor_k <- lapply(unique(panel_df$tahun), function(yy) {
    d_sub <- panel_df |> filter(tahun == yy)
    # Pastikan urutan data sama dengan urutan W
    mi <- moran.test(d_sub$ikp, res$lw, zero.policy = TRUE)
    tibble(tahun = yy, moran_i = mi$estimate[[1]], k = k)
  }) |> bind_rows()
  
  moran_results[[as.character(k)]] <- mor_k
  
  # Simpan objek W
  saveRDS(res$lw, file = sprintf("data/W_knn_k%d_listw.rds", k))
}

# --- 2.B. Pilih W Terbaik ---
summary_w <- bind_rows(moran_results) |>
  group_by(k) |>
  summarise(avg_moran = mean(moran_i)) |>
  arrange(desc(avg_moran))

print(summary_w)

# Misal kita pilih K=10 berdasarkan hasil (sesuaikan jika hasil run berbeda)
k_final <- 10 
message(sprintf("Memilih K=%d sebagai W Final.", k_final))

lw_final <- readRDS(sprintf("data/W_knn_k%d_listw.rds", k_final))
saveRDS(lw_final, "data/W_final_listw.rds")

#===============================================================================
# 3. ANALISIS DIAGNOSTIK & PEMILIHAN MODEL (NON-SPASIAL)
#===============================================================================
message("--- TAHAP 3: Uji Diagnostik Panel Standar ---")

# Setup Data Panel untuk Paket 'plm'
pdata <- plm::pdata.frame(panel_df, index = c("kab_kota", "tahun"))

# Formula Utama
# Sesuaikan variabel independen di sini
vars_x <- c("log_rice","log_nonrice","log_pop","log_grdp","hdi","pov","fiskal","bpk", "ctsr")
form_eq <- as.formula(paste("ikp ~", paste(vars_x, collapse = " + ")))

# --- 3.A. Estimasi Model Standar (Pool, FE, RE) ---
# Pool (OLS Biasa)
mod_pool <- plm(form_eq, data = pdata, model = "pooling")
# Fixed Effect (Within) - Twoways
mod_fe   <- plm(form_eq, data = pdata, model = "within", effect = "twoways")
# Random Effect
mod_re   <- plm(form_eq, data = pdata, model = "random", effect = "twoways")

# --- 3.B. Uji Pemilihan Model ---
message("1. Uji Chow (Pool vs FE)")
print(pFtest(mod_fe, mod_pool)) 
# H0: Pool (Restricted) lebih baik. H1: FE lebih baik.
# Jika p < 0.05 -> Pilih FE.

message("2. Uji Hausman (FE vs RE)")
print(phtest(mod_fe, mod_re))
# H0: RE konsisten (RE lebih efisien). H1: RE tidak konsisten (Pilih FE).
# Jika p < 0.05 -> Pilih FE.

# --- 3.C. Uji Asumsi Klasik (Pada model Pool/LM) ---
lm_obj <- lm(form_eq, data = panel_df)
message("3. Uji Multikolinearitas (VIF)")
print(vif(lm_obj))

message("4. Uji Heteroskedastisitas (Breusch-Pagan)")
print(bptest(lm_obj))

# --- 3.D. Uji Dependensi Spasial (LM Test) ---
# Untuk melihat apakah butuh model spasial atau cukup panel biasa
message("5. Uji LM Lagrange Multiplier (Spatial Dependence)")
# lml = Lag, lme = Error
print(slmtest(mod_pool, listw = lw_final, test = "lml"))
print(slmtest(mod_pool, listw = lw_final, test = "lme"))
print(slmtest(mod_pool, listw = lw_final, test = "rlml")) # Robust
print(slmtest(mod_pool, listw = lw_final, test = "rlme")) # Robust

#===============================================================================
# 4. ESTIMASI MODEL SPASIAL PANEL (SAR, SEM, SDM) - TWOWAYS
#===============================================================================
message("--- TAHAP 4: Estimasi Model Spasial (TWOWAYS FE) ---")
# CATATAN: 'effect = "twoways"' mengontrol efek spesifik individu DAN waktu.

# --- 4.A. Spatial Autoregressive Model (SAR) ---
# Y = rho*W*Y + X*beta + u
message("Estimasi SAR Twoways...")
model_sar_tw <- spml(
  formula       = form_eq,
  data          = panel_df,
  index         = c("kab_kota", "tahun"),
  listw         = lw_final,
  model         = "within",     # Fixed Effect
  effect        = "twoways",    # Individu + Waktu
  lag           = TRUE,         # Ada lag Y (rho)
  spatial.error = "none"        # Tidak ada error spasial
)
summary(model_sar_tw)

# --- 4.B. Spatial Error Model (SEM) ---
# Y = X*beta + u; u = lambda*W*u + e
message("Estimasi SEM Twoways...")
model_sem_tw <- spml(
  formula       = form_eq,
  data          = panel_df,
  index         = c("kab_kota", "tahun"),
  listw         = lw_final,
  model         = "within",
  effect        = "twoways",
  lag           = FALSE,        # Tidak ada lag Y
  spatial.error = "b"           # Baltagi error dependence (SEM)
)
summary(model_sem_tw)

# --- 4.C. Spatial Durbin Model (SDM) ---
# Y = rho*W*Y + X*beta + W*X*theta + u
# Trik: spml tidak punya opsi langsung SDM otomatis, kita harus buat W*X manual
# atau gunakan trik memasukkan lag X ke formula.

message("Estimasi SDM Twoways (Persiapan Lag X)...")

# Buat Lag X (WX)
vars_x_names <- all.vars(form_eq)[-1]
# Fungsi untuk buat lag per tahun (agar tidak nyebrang tahun/wilayah sembarangan)
panel_sdm_ready <- panel_df |>
  group_by(tahun) |>
  mutate(across(all_of(vars_x_names), 
                ~ lag.listw(lw_final, ., zero.policy = TRUE), 
                .names = "lag.{.col}")) |>
  ungroup()

# Buat formula SDM: Y ~ X + lag.X
lag_vars <- paste0("lag.", vars_x_names)
form_sdm <- as.formula(paste("ikp ~", paste(c(vars_x_names, lag_vars), collapse = " + ")))

model_sdm_tw <- spml(
  formula       = form_sdm,
  data          = panel_sdm_ready,
  index         = c("kab_kota", "tahun"),
  listw         = lw_final,
  model         = "within",
  effect        = "twoways",
  lag           = TRUE,         # SDM punya rho (lag Y)
  spatial.error = "none"
)
summary(model_sdm_tw)

#===============================================================================
# 5. PEMILIHAN MODEL TERBAIK (AIC/BIC MANUAL)
#===============================================================================
message("--- TAHAP 5: Perbandingan Model (AIC/BIC) ---")

# Fungsi hitung LogLik, AIC, BIC yang aman untuk objek splm
get_criteria <- function(model, name) {
  # Coba ambil LogLik
  ll <- tryCatch(model$logLik, error = function(e) NA)
  if(is.null(ll) || is.na(ll)) {
    # Fallback jika struktur objek berbeda
    ll <- tryCatch(summary(model)$logLik, error = function(e) NA)
  }
  
  # Hitung parameter (k)
  # Jumlah koefisien + 1 (sigma) + 1 (rho/lambda jika ada)
  k <- length(coef(model)) + 1 
  n <- nrow(model$model) # Jumlah observasi
  
  aic <- -2 * ll + 2 * k
  bic <- -2 * ll + log(n) * k
  
  return(tibble(Model = name, LogLik = as.numeric(ll), AIC = aic, BIC = bic))
}

# Bandingkan
res_sar <- get_criteria(model_sar_tw, "SAR Twoways")
res_sem <- get_criteria(model_sem_tw, "SEM Twoways")
res_sdm <- get_criteria(model_sdm_tw, "SDM Twoways")

final_comparison <- bind_rows(res_sar, res_sem, res_sdm) |>
  arrange(AIC) # Urutkan dari AIC terkecil (Terbaik)

print("TABEL PERBANDINGAN MODEL FINAL:")
print(final_comparison)

# Simpan hasil
write_xlsx(final_comparison, "out_estimasi/HASIL_FINAL_AIC_BIC.xlsx")

# --- KETERANGAN TAMBAHAN UNTUK SKRIPSI ---
cat("\n")
cat("=== CATATAN PENTING ===\n")
cat("1. Model yang digunakan adalah Fixed Effect Two-ways (Within).\n")
cat("2. Two-ways mengontrol efek spesifik wilayah (invididu) DAN efek spesifik waktu (tahun).\n")
cat("3. Jika AIC SDM paling kecil, gunakan SDM. Namun perhatikan signifikansi variabel lag X (W*X).\n")
cat("4. Jangan lupa interpretasi 'Direct' dan 'Indirect' effects jika menggunakan SAR/SDM.\n")

message("✅ PROSES SELESAI SEMPURNA. ALHAMDULILLAH.")
