#===============================================================================
# SPATIAL PANEL DATA ANALYSIS (TWO-WAYS FIXED EFFECTS)
# Interregional Spillover and Food Security - Eastern Indonesia
# Hananta Angger Yuga Prawira
#===============================================================================

#===============================================================================
# 0. CLEAN UP
#===============================================================================
rm(list = ls()); gc(); cat("\014")

#===============================================================================
# 0.1 PACKAGES
#===============================================================================
message("--- TAHAP 0.1: Load packages ---")

library(sf)
library(spdep)
library(splm)

library(tidyverse)
library(janitor)
library(readxl)
library(writexl)

library(plm)
library(lmtest)
library(car)

#===============================================================================
# 0.2 HELPERS
#===============================================================================
message("--- TAHAP 0.2: Helper functions ---")

normalize <- function(x) x |> toupper() |> stringr::str_squish()

fix_geom <- function(g) {
  sf::sf_use_s2(FALSE)
  g_fixed <- if ("st_make_valid" %in% getNamespaceExports("sf")) {
    sf::st_make_valid(g)
  } else sf::st_buffer(g, 0)
  sf::sf_use_s2(TRUE)
  g_fixed
}

build_knn_w <- function(coords, k) {
  nb  <- knn2nb(knearneigh(coords, k = k))
  lw  <- nb2listw(nb, style = "W", zero.policy = TRUE)
  deg <- sapply(nb, length)
  diag <- tibble(
    k = k,
    components = n.comp.nb(nb)$nc,
    min_deg = min(deg),
    mean_deg = mean(deg),
    max_deg = max(deg)
  )
  list(nb = nb, lw = lw, diag = diag)
}

moran_by_year <- function(df, value_col, lw, years, order_col = "order_id") {
  lapply(years, function(yy) {
    d <- df |> filter(tahun == yy) |> arrange(.data[[order_col]])
    mi <- moran.test(d[[value_col]], lw, zero.policy = TRUE)
    tibble(
      tahun   = yy,
      moran_i = unname(mi$estimate[["Moran I statistic"]]),
      p_value = mi$p.value
    )
  }) |> bind_rows()
}

# Elhorst/LeSage-Pace impacts for SDM:
# S = (I - rho W)^(-1)
# Direct  = mean(diag(S)*(beta) + diag(SW)*(theta)) per variable
# Indirect= mean(rowSums(S)*beta + rowSums(SW)*theta - Direct)
sdm_impacts_elhorst <- function(W, rho, beta, theta, n_sims = 0, seed = 1) {
  # W: dense matrix n x n
  # beta: named vector (X)
  # theta: named vector (WX) - same X names ideally
  set.seed(seed)
  n <- nrow(W)
  I <- diag(n)
  
  S <- solve(I - rho * W)         # (I - rho W)^(-1)
  SW <- S %*% W
  
  # diag / row sums
  Sd  <- diag(S)
  SWd <- diag(SW)
  Sr  <- rowSums(S)
  SWr <- rowSums(SW)
  
  vars <- intersect(names(beta), names(theta))
  out <- lapply(vars, function(v) {
    b <- beta[[v]]
    t <- theta[[v]]
    direct  <- mean(Sd * b + SWd * t)
    total   <- mean(Sr * b + SWr * t)
    indirect <- total - direct
    tibble(variable = v, direct = direct, indirect = indirect, total = total)
  }) |> bind_rows()
  
  out
}

#===============================================================================
# 0.3 WORKDIR
#===============================================================================
setwd("D:\\SKRIPSI 26 JANUARI")
dir.create("data", showWarnings = FALSE)
dir.create("out_estimasi", showWarnings = FALSE)
dir.create("out_maps", showWarnings = FALSE)

message("WD: ", getwd())

#===============================================================================
# 1. LOAD DATA (MAP + PANEL) & BUILD PANEL_DF
#===============================================================================
message("--- TAHAP 1: Load & prepare data ---")

# 1.A Load GPKG
gpkg_path <- "petaindotimur.gpkg"
stopifnot(file.exists(gpkg_path))

shp_185 <- st_read(gpkg_path, quiet = TRUE) |>
  clean_names() |>
  fix_geom()

# >>> pastikan kolom nama wilayah ada: kab_kota (ubah jika beda)
stopifnot("kab_kota" %in% names(shp_185))

shp_185 <- shp_185 |> mutate(nm_shp = normalize(kab_kota))

# 1.B Load Excel panel
xls_path <- "DATA SET SKRIPSI.xlsx"
stopifnot(file.exists(xls_path))

panel_raw <- read_xlsx(xls_path, sheet = 1) |> clean_names()

# Rename sesuai data kamu
panel <- panel_raw |>
  select(-any_of("id")) |>
  rename(
    tahun        = tahun,
    kab_kota_raw = kab_kota,
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
    ln_hujan     = ln_hujan,
    ln_land      = ln_land
  ) |>
  mutate(
    tahun = as.integer(tahun),
    nm_list = normalize(kab_kota_raw)
  )

# 1.C Join
panel_sf <- shp_185 |>
  inner_join(panel, by = c("nm_shp" = "nm_list")) |>
  arrange(nm_shp, tahun)

if (nrow(panel_sf) != 925) stop("ERROR: nrow(panel_sf) != 925. Cek join & nama wilayah.")

panel_df <- panel_sf |>
  st_drop_geometry() |>
  mutate(
    id = as.factor(nm_shp),
    tahun = as.integer(tahun)
  ) |>
  arrange(id, tahun)

message("N unit = ", length(unique(panel_df$id)),
        " | T = ", length(unique(panel_df$tahun)),
        " | n = ", nrow(panel_df))

#===============================================================================
# 2. BUILD W CANDIDATES + GLOBAL MORAN (IKP) + SELECT K BY RESIDUAL MORAN
#===============================================================================
message("--- TAHAP 2: Build W candidates + Global Moran + Select K ---")

# 2.0 Order alignment (order_id from shp)
shp_185 <- shp_185 |> mutate(order_id = row_number())

panel_df <- panel_df |>
  left_join(shp_185 |> st_drop_geometry() |> select(nm_shp, order_id), by = "nm_shp")

if (any(is.na(panel_df$order_id))) stop("ERROR: order_id NA. nm_shp mismatch!")

# 2.A Coordinates
sf::sf_use_s2(FALSE)
coords <- st_coordinates(st_point_on_surface(shp_185))
sf::sf_use_s2(TRUE)

years <- sort(unique(panel_df$tahun))
K_list <- 1:10

moran_y_results <- list()
w_diag_results  <- list()

for (k in K_list) {
  res <- build_knn_w(coords, k)
  w_diag_results[[as.character(k)]] <- res$diag
  
  mor_k <- moran_by_year(panel_df, "ikp", res$lw, years, "order_id") |>
    mutate(k = k)
  
  moran_y_results[[as.character(k)]] <- mor_k
  saveRDS(res$lw, sprintf("data/W_knn_k%d_listw.rds", k))
}

summary_moran_y <- bind_rows(moran_y_results) |>
  group_by(k) |>
  summarise(avg_moran = mean(moran_i), avg_p = mean(p_value), .groups="drop") |>
  arrange(desc(avg_moran))

w_diag_tbl <- bind_rows(w_diag_results) |> arrange(k)

print(summary_moran_y)
print(w_diag_tbl)

write_xlsx(summary_moran_y, "out_estimasi/summary_global_moran_IKP_byK.xlsx")
write_xlsx(w_diag_tbl,      "out_estimasi/diag_W_byK.xlsx")

# 2.B Select K final by Moran residual FE two-ways
vars_x <- c("fiskal","bpk","ctsr","ln_poktan","ln_penyuluh","akses","tpt",
            "ln_kud","rs","apm","ln_hujan","ln_land")
form_eq <- as.formula(paste("ikp ~", paste(vars_x, collapse = " + ")))

pdata <- plm::pdata.frame(panel_df, index = c("id","tahun"))
mod_fe_tw <- plm::plm(form_eq, data = pdata, model = "within", effect = "twoways")

moran_resid_results <- list()
for (k in K_list) {
  lw_k <- readRDS(sprintf("data/W_knn_k%d_listw.rds", k))
  tmp <- panel_df |> mutate(resid_fe = as.numeric(residuals(mod_fe_tw)))
  mor_res_k <- moran_by_year(tmp, "resid_fe", lw_k, years, "order_id") |> mutate(k = k)
  moran_resid_results[[as.character(k)]] <- mor_res_k
}

summary_moran_resid <- bind_rows(moran_resid_results) |>
  group_by(k) |>
  summarise(avg_moran_resid = mean(moran_i), avg_p_resid = mean(p_value), .groups="drop") |>
  arrange(abs(avg_moran_resid))

print(summary_moran_resid)
write_xlsx(summary_moran_resid, "out_estimasi/summary_moran_RESID_byK.xlsx")

k_final <- summary_moran_resid$k[1]
message("✅ K FINAL = ", k_final)

lw_final <- readRDS(sprintf("data/W_knn_k%d_listw.rds", k_final))
saveRDS(lw_final, "data/W_final_listw.rds")

#===============================================================================
# 3. LOCAL AUTOCORRELATION (LISA) + MAPPING (ggplot2)
#===============================================================================
message("--- TAHAP 3: LISA (Local Moran) + mapping ---")

# Pilih 1 tahun untuk peta LISA (misal tahun terakhir)
year_lisa <- max(years)

# Ambil IKP tahun tersebut dan urutkan sesuai shp(order_id)
d_y <- panel_df |>
  filter(tahun == year_lisa) |>
  arrange(order_id)

# Local Moran
lmoran <- localmoran(d_y$ikp, lw_final, zero.policy = TRUE)

lisa_df <- d_y |>
  mutate(
    Ii    = lmoran[, "Ii"],
    EIi   = lmoran[, "E.Ii"],
    VarIi = lmoran[, "Var.Ii"],
    ZIi   = lmoran[, "Z.Ii"],
    pIi   = lmoran[, "Pr(z != E(Ii))"]
  )

# Klasifikasi cluster (HH/LL/HL/LH) memakai standardized values
z_y  <- as.numeric(scale(lisa_df$ikp))
wzy  <- as.numeric(lag.listw(lw_final, z_y, zero.policy = TRUE))  # spatial lag of z

lisa_df <- lisa_df |>
  mutate(
    z_y  = z_y,
    wz_y = wzy,
    sig  = pIi <= 0.05,
    quadrant = case_when(
      sig & z_y >= 0 & wz_y >= 0 ~ "HH (High-High)",
      sig & z_y <  0 & wz_y <  0 ~ "LL (Low-Low)",
      sig & z_y >= 0 & wz_y <  0 ~ "HL (High-Low)",
      sig & z_y <  0 & wz_y >= 0 ~ "LH (Low-High)",
      TRUE ~ "Not significant"
    )
  )

# Gabungkan ke shapefile untuk plot
# (shp_185 urut asli order_id, jadi join by nm_shp aman)
shp_lisa <- shp_185 |>
  left_join(
    lisa_df |> select(nm_shp, ikp, Ii, ZIi, pIi, quadrant),
    by = "nm_shp"
  )

# Plot peta cluster LISA (tanpa tmap)
p_lisa <- ggplot(shp_lisa) +
  geom_sf(aes(fill = quadrant), color = NA) +
  labs(
    title = paste0("LISA Cluster Map IKP (", year_lisa, ")"),
    fill = "Cluster"
  ) +
  theme_minimal()

ggsave(filename = sprintf("out_maps/LISA_IKP_%d.png", year_lisa),
       plot = p_lisa, width = 10, height = 7, dpi = 300)

# Simpan tabel LISA
write_xlsx(
  lisa_df |> select(nm_shp, tahun, ikp, Ii, ZIi, pIi, quadrant),
  sprintf("out_estimasi/LISA_table_IKP_%d.xlsx", year_lisa)
)

message("✅ LISA map & table saved.")

#===============================================================================
# 4. LM TEST (Pooled baseline) - spatial dependence
#===============================================================================
message("--- TAHAP 4: LM Test (Pooled baseline) ---")

mod_pool <- plm::plm(form_eq, data = pdata, model = "pooling")

print(slmtest(mod_pool, listw = lw_final, test = "lml"))
print(slmtest(mod_pool, listw = lw_final, test = "lme"))
print(slmtest(mod_pool, listw = lw_final, test = "rlml"))
print(slmtest(mod_pool, listw = lw_final, test = "rlme"))

#===============================================================================
# 5. SPATIAL PANEL FE TWO-WAYS: SAR / SEM / SDM
#===============================================================================
message("--- TAHAP 5: Spatial panel FE two-ways (SAR/SEM/SDM) ---")

# SAR
model_sar_tw <- spml(
  formula       = form_eq,
  data          = panel_df,
  index         = c("id","tahun"),
  listw         = lw_final,
  model         = "within",
  effect        = "twoways",
  lag           = TRUE,
  spatial.error = "none"
)
print(summary(model_sar_tw))

# SEM
model_sem_tw <- spml(
  formula       = form_eq,
  data          = panel_df,
  index         = c("id","tahun"),
  listw         = lw_final,
  model         = "within",
  effect        = "twoways",
  lag           = FALSE,
  spatial.error = "b"
)
print(summary(model_sem_tw))

# SDM (manual WX per tahun)
vars_x_names <- all.vars(form_eq)[-1]
panel_sdm_ready <- panel_df |>
  group_by(tahun) |>
  mutate(across(all_of(vars_x_names),
                ~ lag.listw(lw_final, ., zero.policy = TRUE),
                .names = "lag.{.col}")) |>
  ungroup()

lag_vars <- paste0("lag.", vars_x_names)
form_sdm <- as.formula(paste("ikp ~", paste(c(vars_x_names, lag_vars), collapse = " + ")))

model_sdm_tw <- spml(
  formula       = form_sdm,
  data          = panel_sdm_ready,
  index         = c("id","tahun"),
  listw         = lw_final,
  model         = "within",
  effect        = "twoways",
  lag           = TRUE,
  spatial.error = "none"
)
print(summary(model_sdm_tw))

#===============================================================================
# 6. DIRECT / INDIRECT / TOTAL EFFECTS (ELHORST-STYLE) for SDM
#===============================================================================
message("--- TAHAP 6: Direct/Indirect/Total effects (Elhorst) ---")

stopifnot(exists("lw_final"), exists("model_sdm_tw"))

# 6.1 W matrix NxN
Wmat <- spdep::nb2mat(lw_final$neighbours, style = "W", zero.policy = TRUE)
print(dim(Wmat))
stopifnot(nrow(Wmat) == 185, ncol(Wmat) == 185)

# 6.2 Ambil rho (spatial autoregressive coefficient) - SUPER ROBUST
rho_hat <- numeric(0)

sm <- summary(model_sdm_tw)

# Fallback A: kalau versi object nyimpen langsung
candA <- c(
  tryCatch(as.numeric(model_sdm_tw$lambda), error=function(e) NA_real_),
  tryCatch(as.numeric(model_sdm_tw$rho),    error=function(e) NA_real_),
  tryCatch(as.numeric(sm$lambda),           error=function(e) NA_real_),
  tryCatch(as.numeric(sm$rho),              error=function(e) NA_real_)
)
candA <- candA[is.finite(candA)]
if (length(candA) > 0) rho_hat <- candA[1]

# Fallback B: kalau summary punya table koef spasial
if (length(rho_hat) == 0) {
  # coba beberapa nama slot yang sering muncul lintas versi
  candB <- c(
    tryCatch(as.numeric(sm$spat.coef["lambda","Estimate"]), error=function(e) NA_real_),
    tryCatch(as.numeric(sm$spat.coef["rho","Estimate"]),    error=function(e) NA_real_),
    tryCatch(as.numeric(sm$spat.coef["lambda","estimate"]), error=function(e) NA_real_),
    tryCatch(as.numeric(sm$spat.coef["rho","estimate"]),    error=function(e) NA_real_)
  )
  candB <- candB[is.finite(candB)]
  if (length(candB) > 0) rho_hat <- candB[1]
}

# Fallback C: parse dari teks summary (PASTI ADA di output)
if (length(rho_hat) == 0) {
  txt <- paste(capture.output(print(sm)), collapse = "\n")
  # cari pola "lambda 0.176071" atau "rho 0.176071"
  m <- regexpr("(lambda|rho)\\s+(-?\\d+\\.\\d+)", txt, perl=TRUE)
  if (m[1] != -1) {
    hit <- regmatches(txt, m)
    rho_hat <- as.numeric(sub("^(lambda|rho)\\s+", "", hit))
  }
}

# Fallback D: kalau masih kosong, isi manual dari output summary kamu
if (length(rho_hat) == 0 || !is.finite(rho_hat)) {
  message("⚠️ rho tidak ketemu otomatis. Isi manual dari output SDM (baris 'lambda').")
  # GANTI sesuai output kamu:
  rho_hat <- 0.176071
}

print(rho_hat)
stopifnot(length(rho_hat) == 1, is.finite(rho_hat))

# 6.3 Ambil beta dan theta
coefs <- coef(model_sdm_tw)

vars_x <- c("fiskal","bpk","ctsr","ln_poktan","ln_penyuluh",
            "akses","tpt","ln_kud","rs","apm","ln_hujan","ln_land")

beta_hat <- coefs[vars_x]
theta_hat <- coefs[paste0("lag.", vars_x)]
names(theta_hat) <- vars_x

if (any(is.na(beta_hat)))  stop("ERROR: beta_hat NA. Cek vars_x vs names(coef(model_sdm_tw)).")
if (any(is.na(theta_hat))) stop("ERROR: theta_hat NA. Cek apakah SDM punya lag.<var> lengkap.")

# 6.4 Fungsi impacts (Elhorst)
sdm_impacts_elhorst <- function(W, rho, beta, theta) {
  n <- nrow(W)
  I <- diag(n)
  
  S  <- solve(I - rho * W)
  SW <- S %*% W
  
  Sd  <- diag(S)
  SWd <- diag(SW)
  Sr  <- rowSums(S)
  SWr <- rowSums(SW)
  
  vars <- names(beta)
  
  dplyr::bind_rows(lapply(vars, function(v) {
    direct   <- mean(Sd  * beta[v] + SWd * theta[v])
    total    <- mean(Sr  * beta[v] + SWr * theta[v])
    indirect <- total - direct
    
    tibble::tibble(
      variable = v,
      direct   = direct,
      indirect = indirect,
      total    = total
    )
  }))
}

# 6.5 Hitung impacts
imp_tbl <- sdm_impacts_elhorst(Wmat, rho_hat, beta_hat, theta_hat) |>
  dplyr::arrange(dplyr::desc(abs(indirect)))

print(imp_tbl)

# 6.6 Save
dir.create("out_estimasi", showWarnings = FALSE)
writexl::write_xlsx(imp_tbl, "out_estimasi/SDM_direct_indirect_total_effects.xlsx")

message("✅ Impacts SDM (direct/indirect/total) tersimpan.")
#===============================================================================
# 6. DIRECT / INDIRECT / TOTAL EFFECTS (ELHORST-STYLE) for SDM
#===============================================================================
message("--- TAHAP 6: Direct/Indirect/Total effects (Elhorst) ---")

stopifnot(exists("lw_final"), exists("model_sdm_tw"))

# 6.1 W matrix NxN
Wmat <- spdep::nb2mat(lw_final$neighbours, style = "W", zero.policy = TRUE)
print(dim(Wmat))
stopifnot(nrow(Wmat) == 185, ncol(Wmat) == 185)

# 6.2 Ambil rho (spatial autoregressive coefficient) - SUPER ROBUST
rho_hat <- numeric(0)

sm <- summary(model_sdm_tw)

# Fallback A: kalau versi object nyimpen langsung
candA <- c(
  tryCatch(as.numeric(model_sdm_tw$lambda), error=function(e) NA_real_),
  tryCatch(as.numeric(model_sdm_tw$rho),    error=function(e) NA_real_),
  tryCatch(as.numeric(sm$lambda),           error=function(e) NA_real_),
  tryCatch(as.numeric(sm$rho),              error=function(e) NA_real_)
)
candA <- candA[is.finite(candA)]
if (length(candA) > 0) rho_hat <- candA[1]

# Fallback B: kalau summary punya table koef spasial
if (length(rho_hat) == 0) {
  # coba beberapa nama slot yang sering muncul lintas versi
  candB <- c(
    tryCatch(as.numeric(sm$spat.coef["lambda","Estimate"]), error=function(e) NA_real_),
    tryCatch(as.numeric(sm$spat.coef["rho","Estimate"]),    error=function(e) NA_real_),
    tryCatch(as.numeric(sm$spat.coef["lambda","estimate"]), error=function(e) NA_real_),
    tryCatch(as.numeric(sm$spat.coef["rho","estimate"]),    error=function(e) NA_real_)
  )
  candB <- candB[is.finite(candB)]
  if (length(candB) > 0) rho_hat <- candB[1]
}

# Fallback C: parse dari teks summary (PASTI ADA di output)
if (length(rho_hat) == 0) {
  txt <- paste(capture.output(print(sm)), collapse = "\n")
  # cari pola "lambda 0.176071" atau "rho 0.176071"
  m <- regexpr("(lambda|rho)\\s+(-?\\d+\\.\\d+)", txt, perl=TRUE)
  if (m[1] != -1) {
    hit <- regmatches(txt, m)
    rho_hat <- as.numeric(sub("^(lambda|rho)\\s+", "", hit))
  }
}

# Fallback D: kalau masih kosong, isi manual dari output summary kamu
if (length(rho_hat) == 0 || !is.finite(rho_hat)) {
  message("⚠️ rho tidak ketemu otomatis. Isi manual dari output SDM (baris 'lambda').")
  # GANTI sesuai output kamu:
  rho_hat <- 0.176071
}

print(rho_hat)
stopifnot(length(rho_hat) == 1, is.finite(rho_hat))

# 6.3 Ambil beta dan theta
coefs <- coef(model_sdm_tw)

vars_x <- c("fiskal","bpk","ctsr","ln_poktan","ln_penyuluh",
            "akses","tpt","ln_kud","rs","apm","ln_hujan","ln_land")

beta_hat <- coefs[vars_x]
theta_hat <- coefs[paste0("lag.", vars_x)]
names(theta_hat) <- vars_x

if (any(is.na(beta_hat)))  stop("ERROR: beta_hat NA. Cek vars_x vs names(coef(model_sdm_tw)).")
if (any(is.na(theta_hat))) stop("ERROR: theta_hat NA. Cek apakah SDM punya lag.<var> lengkap.")

# 6.4 Fungsi impacts (Elhorst)
sdm_impacts_elhorst <- function(W, rho, beta, theta) {
  n <- nrow(W)
  I <- diag(n)

  S  <- solve(I - rho * W)
  SW <- S %*% W

  Sd  <- diag(S)
  SWd <- diag(SW)
  Sr  <- rowSums(S)
  SWr <- rowSums(SW)

  vars <- names(beta)

  dplyr::bind_rows(lapply(vars, function(v) {
    direct   <- mean(Sd  * beta[v] + SWd * theta[v])
    total    <- mean(Sr  * beta[v] + SWr * theta[v])
    indirect <- total - direct

    tibble::tibble(
      variable = v,
      direct   = direct,
      indirect = indirect,
      total    = total
    )
  }))
}

# 6.5 Hitung impacts
imp_tbl <- sdm_impacts_elhorst(Wmat, rho_hat, beta_hat, theta_hat) |>
  dplyr::arrange(dplyr::desc(abs(indirect)))

print(imp_tbl)

# 6.6 Save
dir.create("out_estimasi", showWarnings = FALSE)
writexl::write_xlsx(imp_tbl, "out_estimasi/SDM_direct_indirect_total_effects.xlsx")

message("✅ Impacts SDM (direct/indirect/total) tersimpan.")

message("--- TAHAP 7: AIC, BIC, dan Pseudo R² (SAR / SEM / SDM) ---")

# 7.1 Fungsi umum hitung AIC & BIC manual
calc_aic_bic <- function(model, k, n) {
  ll <- as.numeric(model$logLik)
  if (is.null(ll) || is.na(ll)) {
    ll <- as.numeric(summary(model)$logLik)
  }
  AIC <- -2 * ll + 2 * k
  BIC <- -2 * ll + log(n) * k
  c(AIC = AIC, BIC = BIC)
}

# 7.2 Fungsi Pseudo R² (Elhorst-style)
calc_pseudo_r2 <- function(model, y) {
  u <- as.numeric(residuals(model))
  1 - sum(u^2, na.rm = TRUE) /
    sum((y - mean(y, na.rm = TRUE))^2, na.rm = TRUE)
}

#===============================================================================
# 7.3 Komponen dasar
#===============================================================================
n_obs <- nrow(panel_df)     # jumlah observasi
y_var <- panel_df$ikp       # variabel dependen

# variabel X (harus sama dengan form_eq)
vars_x <- c("fiskal","bpk","ctsr","ln_poktan","ln_penyuluh",
            "akses","tpt","ln_kud","rs","apm","ln_hujan","ln_land")

#===============================================================================
# 7.4 Hitung jumlah parameter (k)
#===============================================================================
# SAR: beta + rho + sigma^2
k_sar <- length(vars_x) + 2

# SEM: beta + lambda + sigma^2
k_sem <- length(vars_x) + 2

# SDM: beta + theta + rho + sigma^2
k_sdm <- length(vars_x) * 2 + 2

#===============================================================================
# 7.5 AIC & BIC
#===============================================================================
aicbic_sar <- calc_aic_bic(model_sar_tw, k_sar, n_obs)
aicbic_sem <- calc_aic_bic(model_sem_tw, k_sem, n_obs)
aicbic_sdm <- calc_aic_bic(model_sdm_tw, k_sdm, n_obs)

#===============================================================================
# 7.6 Pseudo R²
#===============================================================================
r2_sar <- calc_pseudo_r2(model_sar_tw, y_var)
r2_sem <- calc_pseudo_r2(model_sem_tw, y_var)
r2_sdm <- calc_pseudo_r2(model_sdm_tw, y_var)

#===============================================================================
# 7.7 TABEL FINAL (BAB IV READY)
#===============================================================================
model_comp_tbl <- tibble(
  Model = c("SAR FE Two-ways", "SEM FE Two-ways", "SDM FE Two-ways"),
  AIC   = c(aicbic_sar["AIC"], aicbic_sem["AIC"], aicbic_sdm["AIC"]),
  BIC   = c(aicbic_sar["BIC"], aicbic_sem["BIC"], aicbic_sdm["BIC"]),
  Pseudo_R2 = c(r2_sar, r2_sem, r2_sdm)
) |>
  arrange(AIC)

print(model_comp_tbl)

# Simpan ke Excel
dir.create("out_estimasi", showWarnings = FALSE)
writexl::write_xlsx(model_comp_tbl,
                    "out_estimasi/Model_Comparison_AIC_BIC_R2.xlsx")

message("✅ Tabel perbandingan SAR–SEM–SDM selesai & tersimpan.")

