#===============================================================================
#   PROJECT BARU: SPATIAL PANEL (TWOWAYS FE) - PROVINSI
#   Y & X MENGIKUTI SPREADSHEET: hasil_final_id_dan_nama_provinsi.xlsx
#   MATRIKS W: K-Nearest Neighbors (KNN)
#===============================================================================

# 0.0. BERSIH-BERSIH MEMORY
rm(list = ls())
gc()
cat("\014")

#===============================================================================
# 0.1. LOAD PUSTAKA
#===============================================================================
message("--- TAHAP 0.1: Memuat Library ---")

library(sf)
library(tmap)
library(classInt)
library(spdep)
library(splm)

library(tidyverse)
library(janitor)
library(readxl)
library(writexl)

library(plm)
library(lmtest)
library(car)
library(sandwich)
library(broom)

#===============================================================================
# 0.2. HELPER FUNCTIONS
#===============================================================================
message("--- TAHAP 0.2: Definisi Fungsi Helper ---")

normalize <- function(x) {
  x |> as.character() |> toupper() |> stringr::str_squish()
}

fix_geom <- function(g) {
  sf::sf_use_s2(FALSE)
  if ("st_make_valid" %in% getNamespaceExports("sf")) {
    g_fixed <- sf::st_make_valid(g)
  } else {
    g_fixed <- sf::st_buffer(g, 0)
  }
  sf::sf_use_s2(TRUE)
  return(g_fixed)
}

build_knn_w <- function(coords, k) {
  nb  <- knn2nb(knearneigh(coords, k = k))
  lw  <- nb2listw(nb, style = "W", zero.policy = TRUE)
  deg <- sapply(nb, length)

  tib <- tibble(
    k = k,
    components = n.comp.nb(nb)$nc,
    min_deg = min(deg),
    mean_deg = mean(deg),
    max_deg = max(deg)
  )
  list(nb = nb, lw = lw, diag = tib)
}

#===============================================================================
# 0.3. SET WORKDIR
#===============================================================================
setwd("~/SKRIPSI")  # sesuaikan
dir.create("data", showWarnings = FALSE)
dir.create("out_estimasi", showWarnings = FALSE)
message("Direktori kerja saat ini: ", getwd())

#===============================================================================
# 1. LOAD DATA PANEL (EXCEL) - PROJECT BARU
#===============================================================================
message("--- TAHAP 1: Load Data Panel (Provinsi) ---")

xls_path  <- "hasil_final_id_dan_nama_provinsi.xlsx"
panel_raw <- read_xlsx(xls_path, sheet = 1) |> clean_names()

# Standarisasi nama kolom kunci agar konsisten
# (clean_names() akan bikin jadi: prov_id, provinsi, tahun, dst.)
panel <- panel_raw |>
  rename(
    id    = prov_id,
    prov  = provinsi,
    tahun = tahun
  ) |>
  mutate(
    id = as.integer(id),
    prov_key = normalize(prov),
    tahun = as.integer(tahun)
  ) |>
  arrange(id, tahun)

# --- cek dimensi panel (harus balanced: 34 prov x 9 tahun = 306) ---
n_id    <- n_distinct(panel$id)
n_tahun <- n_distinct(panel$tahun)
message(sprintf("Panel check: %d provinsi x %d tahun = %d observasi",
                n_id, n_tahun, nrow(panel)))

#===============================================================================
# 1B. PILIH Y DAN X SESUAI SPREADSHEET
#===============================================================================
message("--- TAHAP 1B: Definisi Y dan X (mengikuti spreadsheet) ---")

# Kandidat Y: kolom yang diawali y_ atau ln_y_
y_candidates <- names(panel) |> stringr::str_subset("^(y_|ln_y_)")
x_candidates <- names(panel) |> stringr::str_subset("^(x_|ln_x_)")

message("Kandidat Y (dari spreadsheet):")
print(y_candidates)

message("Kandidat X (dari spreadsheet):")
print(x_candidates)

# ---- DEFAULT (bisa kamu ganti kalau mau) ----
# Y utama: log deforestasi GFW (lebih stabil)
y_var <- if ("ln_y_deforestasi_gfw" %in% names(panel)) "ln_y_deforestasi_gfw" else
  if ("ln_y_deforestasi_pos" %in% names(panel)) "ln_y_deforestasi_pos" else
    y_candidates[1]

# X utama: pilih yang umum dan menghindari kolinearitas "total vs sektoral"
# (kamu punya ln_x_pdrb_total dan juga pdrb sektoral; biasanya pilih salah satu grup)
x_vars <- c(
  "ln_x_pma_pc",
  "ln_x_pdrb_total_pc",
  "x_ipm",
  "x_iklh",
  "x_iktl",
  "ln_x_kepadatan_penduduk",
  "ln_x_luaswilayah_km2"
)

# drop yang tidak ada (biar aman kalau kolom beda)
x_vars <- x_vars[x_vars %in% names(panel)]

message("Y dipakai:")
print(y_var)

message("X dipakai:")
print(x_vars)

if (length(x_vars) == 0) stop("ERROR: Tidak ada X yang cocok. Cek nama kolom di spreadsheet.")

#===============================================================================
# 1C. LOAD SHAPEFILE PROVINSI + JOIN
#===============================================================================
message("--- TAHAP 1C: Load Shapefile Provinsi + Join ---")

# NOTE:
# Kamu perlu shapefile level PROVINSI.
# Misal: data/INDO_PROVINSI.shp atau data/indonesia_prov.gpkg
# Di bawah ini aku bikin fleksibel: coba baca GPKG dulu, kalau tidak ada baca SHP.

prov_gpkg <- "data/peta_provinsi.gpkg"
prov_shp  <- "data/PROVINSI_INDONESIA.shp"  # <-- ganti sesuai file kamu

if (file.exists(prov_gpkg)) {
  shp_prov <- st_read(prov_gpkg, quiet = TRUE) |> clean_names()
} else if (file.exists(prov_shp)) {
  shp_prov <- st_read(prov_shp, quiet = TRUE) |> clean_names()
} else {
  stop("ERROR: Shapefile/GPKG provinsi tidak ditemukan. Simpan di data/ lalu set path prov_shp/prov_gpkg.")
}

# Cari kolom nama provinsi di shapefile (umum: provinsi, nama_prov, nm_prov, dll.)
# Kamu bisa print(names(shp_prov)) untuk cek.
name_candidates <- c("provinsi", "nama_prov", "nm_prov", "province", "wadmpr", "propinsi")
name_col <- name_candidates[name_candidates %in% names(shp_prov)][1]
if (is.na(name_col)) stop("ERROR: Tidak ketemu kolom nama provinsi di shapefile. Cek names(shp_prov).")

shp_prov <- shp_prov |>
  mutate(prov_key = normalize(.data[[name_col]])) |>
  fix_geom()

panel_sf <- shp_prov |>
  inner_join(panel, by = "prov_key") |>
  arrange(id, tahun)

# validasi join (harus sama dengan nrow(panel))
if (nrow(panel_sf) != nrow(panel)) {
  message(sprintf("WARNING: Baris setelah join = %d, sebelum join = %d", nrow(panel_sf), nrow(panel)))
  message("Cek penulisan nama provinsi di shapefile vs spreadsheet (normalisasi sudah dipakai).")
}

saveRDS(panel_sf, "data/panel_sf_prov.rds")
saveRDS(st_drop_geometry(panel_sf), "data/panel_df_prov.rds")

#===============================================================================
# 2. PEMBENTUKAN MATRIKS W (KNN)
#===============================================================================
message("--- TAHAP 2: Seleksi Matriks W (KNN) ---")

# IMPORTANT:
# KNN harus berdasarkan 1 geometri per provinsi (bukan per tahun).
# Jadi kita ambil unik provinsi dari shapefile, urutannya harus konsisten dengan panel.

shp_prov_unique <- shp_prov |>
  arrange(prov_key) |>
  distinct(prov_key, .keep_all = TRUE)

sf::sf_use_s2(FALSE)
coords <- st_coordinates(st_point_on_surface(shp_prov_unique))
sf::sf_use_s2(TRUE)

K_list <- 3:10  # provinsi lebih sedikit, K terlalu besar bisa “padat”
moran_results <- list()

# bikin data cross-section untuk moran: ambil 1 tahun (atau rata-rata) biar konsisten
# Cara 1: pakai tahun pertama yang tersedia
tahun_ref <- sort(unique(panel$tahun))[1]
cs_ref <- panel |>
  filter(tahun == tahun_ref) |>
  mutate(prov_key = normalize(prov)) |>
  arrange(prov_key)

# pastikan urutan sama dengan shp_prov_unique
cs_ref <- cs_ref |>
  semi_join(shp_prov_unique |> st_drop_geometry(), by = "prov_key") |>
  arrange(prov_key)

if (nrow(cs_ref) != nrow(shp_prov_unique)) {
  message("WARNING: cross-section tahun_ref tidak sama jumlahnya dengan shapefile unique provinsi.")
  message("Pastikan semua provinsi ada pada tahun_ref atau pakai rata-rata semua tahun.")
}

for (k in K_list) {
  res <- build_knn_w(coords, k)

  # Moran untuk Y pada cross-section (tahun_ref)
  mi <- moran.test(cs_ref[[y_var]], res$lw, zero.policy = TRUE)
  moran_results[[as.character(k)]] <- tibble(
    k = k,
    moran_i = mi$estimate[[1]],
    p_value = mi$p.value,
    components = res$diag$components
  )

  saveRDS(res$lw, file = sprintf("data/W_knn_k%d_listw.rds", k))
}

summary_w <- bind_rows(moran_results) |>
  arrange(desc(moran_i))

print(summary_w)

# pilih K final: default ambil Moran I tertinggi yang components=1 (kalau ada)
pick <- summary_w |>
  filter(components == 1) |>
  slice_max(order_by = moran_i, n = 1)

k_final <- if (nrow(pick) == 1) pick$k else summary_w$k[1]
message(sprintf("Memilih K=%d sebagai W Final.", k_final))

lw_final <- readRDS(sprintf("data/W_knn_k%d_listw.rds", k_final))
saveRDS(lw_final, "data/W_final_listw.rds")

#===============================================================================
# 3. PANEL DIAGNOSTIC (NON-SPASIAL)
#===============================================================================
message("--- TAHAP 3: Uji Diagnostik Panel Standar ---")

panel_df <- panel_sf |> st_drop_geometry()

pdata <- plm::pdata.frame(panel_df, index = c("id", "tahun"))

form_eq <- as.formula(paste(y_var, "~", paste(x_vars, collapse = " + ")))
message("Formula utama:")
print(form_eq)

mod_pool <- plm(form_eq, data = pdata, model = "pooling")
mod_fe   <- plm(form_eq, data = pdata, model = "within", effect = "twoways")
mod_re   <- plm(form_eq, data = pdata, model = "random", effect = "twoways")

message("1) Uji Chow (Pool vs FE)")
print(pFtest(mod_fe, mod_pool))

message("2) Uji Hausman (FE vs RE)")
print(phtest(mod_fe, mod_re))

lm_obj <- lm(form_eq, data = panel_df)

message("3) VIF (Multikolinearitas)")
print(vif(lm_obj))

message("4) Breusch-Pagan (Heteroskedastisitas)")
print(bptest(lm_obj))

message("5) LM test dependensi spasial (berbasis model pool)")
print(slmtest(mod_pool, listw = lw_final, test = "lml"))
print(slmtest(mod_pool, listw = lw_final, test = "lme"))
print(slmtest(mod_pool, listw = lw_final, test = "rlml"))
print(slmtest(mod_pool, listw = lw_final, test = "rlme"))

#===============================================================================
# 4. SPATIAL PANEL (SAR / SEM / SDM) - TWOWAYS FE
#===============================================================================
message("--- TAHAP 4: Estimasi Model Spasial (TWOWAYS FE) ---")

message("Estimasi SAR Twoways...")
model_sar_tw <- spml(
  formula       = form_eq,
  data          = panel_df,
  index         = c("id", "tahun"),
  listw         = lw_final,
  model         = "within",
  effect        = "twoways",
  lag           = TRUE,
  spatial.error = "none"
)
print(summary(model_sar_tw))

message("Estimasi SEM Twoways...")
model_sem_tw <- spml(
  formula       = form_eq,
  data          = panel_df,
  index         = c("id", "tahun"),
  listw         = lw_final,
  model         = "within",
  effect        = "twoways",
  lag           = FALSE,
  spatial.error = "b"
)
print(summary(model_sem_tw))

message("Estimasi SDM Twoways (buat WX)...")

vars_x_names <- all.vars(form_eq)[-1]

panel_sdm_ready <- panel_df |>
  group_by(tahun) |>
  mutate(across(all_of(vars_x_names),
                ~ lag.listw(lw_final, ., zero.policy = TRUE),
                .names = "lag.{.col}")) |>
  ungroup()

lag_vars <- paste0("lag.", vars_x_names)
form_sdm <- as.formula(paste(y_var, "~", paste(c(vars_x_names, lag_vars), collapse = " + ")))

model_sdm_tw <- spml(
  formula       = form_sdm,
  data          = panel_sdm_ready,
  index         = c("id", "tahun"),
  listw         = lw_final,
  model         = "within",
  effect        = "twoways",
  lag           = TRUE,
  spatial.error = "none"
)
print(summary(model_sdm_tw))

#===============================================================================
# 5. PERBANDINGAN MODEL (AIC/BIC)
#===============================================================================
message("--- TAHAP 5: Perbandingan Model (AIC/BIC) ---")

get_criteria <- function(model, name) {
  ll <- tryCatch(model$logLik, error = function(e) NA)
  if (is.null(ll) || is.na(ll)) ll <- tryCatch(summary(model)$logLik, error = function(e) NA)

  k <- length(coef(model)) + 1
  n <- nrow(model$model)

  aic <- -2 * ll + 2 * k
  bic <- -2 * ll + log(n) * k

  tibble(Model = name, LogLik = as.numeric(ll), AIC = aic, BIC = bic)
}

res_sar <- get_criteria(model_sar_tw, "SAR Twoways")
res_sem <- get_criteria(model_sem_tw, "SEM Twoways")
res_sdm <- get_criteria(model_sdm_tw, "SDM Twoways")

final_comparison <- bind_rows(res_sar, res_sem, res_sdm) |>
  arrange(AIC)

print(final_comparison)
write_xlsx(final_comparison, "out_estimasi/HASIL_FINAL_AIC_BIC_PROV.xlsx")

message("✅ PROSES SELESAI.")
