#================================================================
#
# SKRIPSI HANANTA ANGGER YUGA PRAWIRA
# TUJUAN : HETEROGENITAS SPASIAL (GWPR)
# VERSI : REVISI (Menghapus Log, Menggunakan Shapefile untuk Koordinat)
#
# Ringkasan Logika BARU:
# 1. Memuat data dan shapefile. Melakukan MERGE data panel dengan shapefile
#    berdasarkan KABKOT/KAB_KOTA untuk mendapatkan koordinat resmi.
# 2. Menjalankan Uji Panel Global (Chow, Hausman) -> Menemukan FEM.
# 3. Data di-demean (transformasi "within") untuk SEMUA variabel model KONTINU,
#    kecuali variabel dummy (BPK).
# 4. Menjalankan GWR pada data yang sudah di-demean (ini adalah GWPR).
#
#================================================================

#----------------------------------------------------------------
# Bagian 0: SETUP DAN MEMUAT LIBRARIES
#----------------------------------------------------------------

#--- Bersihkan Environment ---
rm(list = ls())

#--- Opsi Tampilan ---
options(scipen = 999)

#--- Memuat Libraries ---
suppressPackageStartupMessages({
  library(readxl)
  library(plm)
  library(lmtest)
  library(tseries)
  library(GWmodel)
  library(spgwr)
  library(spdep)
  library(stats)
  library(sf)     # PENTING: Untuk penanganan shapefile dan merge spasial
  library(psych)
  library(writexl)
  library(sp)     # Diperlukan untuk GWmodel
})

#--- PENGATURAN WORKING DIRECTORY (PENTING!) ---
# REKOMENDASI: Gunakan R Project (.Rproj).

# Tentukan path file dan nama sheet
file_path <- "Data Set SKRIPSI.xlsx"
sheet_name <- "Sheet1"

# PENTING: Tentukan path ke shapefile Indonesia
# GANTI INI DENGAN PATH YANG SESUAI DI KOMPUTER ANDA
shapefile_path <- "INDONESIA.shp"

# Tentukan nama kolom kordinat yang akan digunakan setelah merge
col_longitude <- "longitude_shp" # Diubah ke huruf kecil
col_latitude <- "latitude_shp"   # Diubah ke huruf kecil


#================================================================
# Bagian 1: IMPORT DATA DAN MERGE SPASIAL
#================================================================

cat("\n--- Bagian 1: Memuat Data dan Melakukan Merge Spasial ---\n")

# 1. Memuat Data Panel dari Excel
tryCatch({
  data.panel <- read_excel(file_path, sheet = sheet_name)
  data.panel <- as.data.frame(data.panel)
  cat("Sukses: Data panel berhasil dimuat dari:", file_path, "\n")
}, error = function(e) {
  stop("ERROR: File data tidak ditemukan. Pastikan file '", file_path,
       "' ada di folder yang benar. Pesan error asli: ", e$message)
})

# --- FIX KRITIS: STANDARDISASI NAMA KOLOM KE HURUF KECIL ---
cat("Standardisasi: Mengubah semua nama kolom data panel menjadi huruf kecil (untuk menghindari error 'capslock').\n")
names(data.panel) <- tolower(names(data.panel))
# --- AKHIR FIX KRITIS ---

# 2. Memuat Shapefile
tryCatch({
  # Gunakan st_read dari paket sf
  shp.indo <- st_read(shapefile_path)
  cat("Sukses: Shapefile berhasil dimuat dari:", shapefile_path, "\n")
}, error = function(e) {
  stop("ERROR: Shapefile tidak ditemukan. Pastikan path '", shapefile_path,
       "' benar. Pesan error asli: ", e$message)
})

# --- FIX UNTUK ERROR GEOMETRI ---
cat("Memeriksa dan memperbaiki geometri (st_make_valid) sebelum menghitung centroid...\n")
shp.indo <- st_make_valid(shp.indo)
# --- AKHIR FIX ---

# Standardisasi nama kolom shapefile ke huruf kecil untuk merge yang mulus
names(shp.indo) <- tolower(names(shp.indo))
shp_key_col <- "kab_kota" # Kolom kunci shapefile dijamin huruf kecil

cat("\nStruktur Shapefile Awal (Kolom Kunci: kab_kota):\n")
print(head(shp.indo$kab_kota))
cat("Struktur Data Panel Awal (Kolom Kunci: kabkot):\n")
print(head(data.panel$kab_kota))

# 3. Ekstrak Koordinat Centroid dari Shapefile
cat("\nEkstraksi koordinat centroid dari shapefile...\n")
shp.coords <- st_centroid(shp.indo)
coords_df <- as.data.frame(st_coordinates(shp.coords))
st_geometry(shp.coords) <- NULL # Hapus geometri untuk merge

# Gabungkan data koordinat centroid dengan data atribut shapefile
shp.final <- cbind(st_drop_geometry(shp.indo), coords_df)
# Ganti nama kolom koordinat
names(shp.final)[names(shp.final) == "X"] <- col_longitude
names(shp.final)[names(shp.final) == "Y"] <- col_latitude

# Ambil hanya kolom kunci (kab_kota) dan koordinat yang diperlukan
shp.final <- shp.final[, c(shp_key_col, col_longitude, col_latitude)]

# 4. Melakukan Merge (Data Panel + Koordinat)
cat("\nMelakukan merge data panel dengan koordinat shapefile berdasarkan kabkot/kab_kota...\n")
# Merge: data.panel (key: kabkot) dan shp.final (key: kab_kota)
data.merged <- merge(
  data.panel,
  shp.final,
  by.X = "kab_kota", # Diubah ke huruf kecil
  by.Y = shp_key_col,
  all.x = TRUE
)

# 5. Pembersihan NA dan Finalisasi Data
cat("\nJumlah NA per variabel (SEBELUM dihapus):\n")
print(colSums(is.na(data.merged)))

# Hapus baris yang memiliki NA di kolom koordinat atau variabel penting
rows_before <- nrow(data.merged)
data.gwpr.clean <- na.omit(data.merged)
rows_after <- nrow(data.gwpr.clean)

cat("\nPERINGATAN: ", rows_before - rows_after, " baris data telah dihapus karena mengandung NA (termasuk yang gagal di-merge).\n")
cat("Jumlah data tersisa: ", rows_after, " observasi.\n")

if (!all(c(col_longitude, col_latitude) %in% names(data.gwpr.clean))) {
  stop("ERROR PENTING: Kolom koordinat hasil merge tidak ditemukan. Pastikan kolom 'kab_kota' di shapefile dan 'kabkot' di data Excel memiliki format yang sama.")
}

#================================================================
# Bagian 2: FORMULA MODEL DAN PEMERIKSAAN DATA
#================================================================
cat("\n--- Bagian 2: Menetapkan Formula Model dan Memproses Log Selektif ---\n")

# --- FIX UNTUK MEMENUHI KEBUTUHAN VARIABEL 'ln_...' YANG TERCANTUM DI MODEL ---
cat("\n--- Perbaikan: Melakukan Transformasi Logaritmik Selektif untuk variabel 'ln_' ---\n")
# Asumsi: Kolom sumber tanpa 'ln_' (poktan, penyuluh, kud, lahan) ada di data Excel (sekarang sudah huruf kecil).
vars_to_log_source <- c("poktan", "penyuluh", "kud", "lahan")

for (var in vars_to_log_source) {
  ln_var_name <- paste0("ln_", var)
  # Hanya buat variabel log jika variabel tersebut belum ada DAN variabel sumbernya ada
  if (!(ln_var_name %in% names(data.gwpr.clean))) {
    if (var %in% names(data.gwpr.clean)) {
      cat(paste("Membuat variabel", ln_var_name, "dari", var, "menggunakan log(X + 1e-9)...\n"))
      # Menggunakan log(X + epsilon) untuk mencegah log(0)
      data.gwpr.clean[[ln_var_name]] <- log(data.gwpr.clean[[var]] + 1e-9)
    } else {
      cat(paste("PERINGATAN KRITIS: Kolom sumber '", var, "' tidak ditemukan. Model GWPR mungkin GAGAL jika variabel ini dibutuhkan.\n"))
    }
  }
}
# --- AKHIR FIX LOG SELEKTIF ---


# Daftar Variabel BARU (ikp sebagai Y, sisanya X):
# "ikp", "fiskal", "bpk", "ctsr", "ln_poktan", "ln_penyuluh", "akses",
# "ln_kud", "rs", "apm", "rain", "ln_lahan"

MODEL_FORMULA <- ikp ~ fiskal + bpk + ctsr + ln_poktan + ln_penyuluh +
  akses + tpt + ln_kud + rs + apm + ln_hujan + ln_land # Semua variabel dijamin huruf kecil

cat("\nFormula Model BARU yang Digunakan untuk SEMUA analisis:\n")
print(MODEL_FORMULA)

cat("\nStatistik Deskriptif Variabel Model:\n")
vars_numeric <- all.vars(MODEL_FORMULA)
print(describe(data.gwpr.clean[vars_numeric]))


#================================================================
# Bagian 3: MODEL PANEL GLOBAL (PENENTUAN MODEL TERBAIK)
#================================================================
cat("\n--- Bagian 3: Menganalisis Model Panel Global ---\n")

pdata_global <- pdata.frame(data.gwpr.clean, index = c("id", "tahun")) # ID dan Tahun diubah ke huruf kecil

# 1. Pooled OLS Model
model_pooled <- plm(MODEL_FORMULA, data = pdata_global, model = "pooling")

# 2. Fixed Effect Model (FEM / Within)
model_fem <- plm(MODEL_FORMULA, data = pdata_global, model = "within")

# 3. Random Effect Model (REM)
model_rem <- plm(MODEL_FORMULA, data = pdata_global, model = "random")

# Uji 1: Uji Chow (FEM vs Pooled OLS)
cat("\n--- UJI CHOW (FEM vs Pooled OLS) ---\n")
uji.chow <- pFtest(model_fem, model_pooled)
print(uji.chow)
cat("Keputusan: Jika P-value < 0.05, pilih FEM.\n")

# Uji 2: UJI HAUSMAN (FEM VS REM)
cat("\n--- UJI HAUSMAN (FEM vs REM) ---\n")
uji.hausman <- phtest(model_fem, model_rem)
print(uji.hausman)
cat("Keputusan: Jika P-value < 0.05, pilih FEM.\n")

MODEL_TERPILIH <- model_fem
cat("\nModel terpilih untuk Uji Asumsi: model_fem\n")

# (Lanjutan Uji Asumsi model_fem - Sama seperti kode lama, hanya menyesuaikan formula)
# ... (Blok Uji Asumsi dipertahankan)
cat("\n--- UJI ASUMSI PADA MODEL PANEL TERPILIH (FEM) ---\n")
galat_model_terpilih <- residuals(MODEL_TERPILIH)

cat("--- 1. Uji Normalitas Residuals ---\n")
if (length(galat_model_terpilih) < 5000) {
  cat("Shapiro-Wilks Test (Sampel Kecil):\n")
  print(shapiro.test(galat_model_terpilih[1:min(length(galat_model_terpilih), 4999)]))
}
cat("Jarque-Bera Test:\n")
print(jarque.bera.test(galat_model_terpilih))
cat("Keputusan: P-value < 0.05, artinya residual tidak terdistribusi normal (umum terjadi).\n")

cat("\n--- 2. Uji Homoskedastisitas (Breusch-Pagan) ---\n")
uji.bp <- bptest(MODEL_FORMULA, data = data.gwpr.clean, studentize = FALSE)
print(uji.bp)
cat("Keputusan: P-value < 0.05, artinya terdapat HETEROSKEDASTISITAS (menguatkan GWPR).\n")

cat("\nKoefisien Model Terpilih dengan Koreksi Heteroskedastisitas (HC1):\n")
print(coeftest(MODEL_TERPILIH, vcovHC(MODEL_TERPILIH, type = "HC1")))

cat("\n--- 3. Uji Autokorelasi ---\n")
cat("Wooldridge Test:\n")
print(pwartest(MODEL_TERPILIH))
cat("Breusch-Godfrey Test:\n")
print(pbgtest(MODEL_TERPILIH))
cat("Keputusan: P-value < 0.05, artinya terdapat Autokorelasi Serial.\n")
cat("CATATAN: Pelanggaran asumsi adalah alasan kuat untuk menggunakan model GWPR, yang lebih fleksibel.\n")


#================================================================
# Bagian 4: PERSIAPAN DATA SPASIAL (Panel GWPR) - DEMEANING
#================================================================
cat("\n--- Bagian 4: Mempersiapkan Data untuk GWPR (Demeaning) ---\n")

# GWPR-FEM membutuhkan data yang di-"demean" (transformasi "within").
all_model_vars <- all.vars(MODEL_FORMULA)

# Variabel yang TIDAK boleh di-demean (variabel dummy). Saat ini hanya 'bpk'.
vars_to_exclude_from_demeaning <- c("bpk") 
# Variabel yang HARUS di-demean (Y dan X kontinu/log-transform)
vars_to_demean <- setdiff(all_model_vars, vars_to_exclude_from_demeaning) 

data.trans.gwpr <- data.gwpr.clean

cat("Melakukan demeaning (transformasi 'within') untuk variabel KONTINU:\n",
    paste(vars_to_demean, collapse = ", "), "\n")
cat("Variabel DUMMY yang dipertahankan (TIDAK di-demean):", 
    paste(vars_to_exclude_from_demeaning, collapse = ", "), "\n")

for (var in vars_to_demean) {
  if (var %in% names(data.trans.gwpr)) {
    # Hitung rata-rata kelompok (id)
    group_means <- ave(data.trans.gwpr[[var]],
                       data.trans.gwpr$id, # id diubah ke huruf kecil
                       FUN = function(x) mean(x, na.rm = TRUE))
    
    # Demean: Variabel dikurangi rata-rata kelompok
    data.trans.gwpr[[paste0("demean_", var)]] <- data.trans.gwpr[[var]] - group_means # Hasil demeaning juga huruf kecil
  } else {
    cat("PERINGATAN: Variabel '", var, "' tidak ditemukan untuk demeaning.\n")
  }
}

# 1. Persiapan Data Frame Final untuk GWPR
# Ambil variabel yang sudah di-demean (Y + X kontinu)
demeaned_cols <- paste0("demean_", vars_to_demean)
# Ambil variabel yang dikecualikan (X dummy)
kept_cols <- vars_to_exclude_from_demeaning

# Ambil juga kolom ID, Tahun, KABKOT, dan Koordinat
meta_cols <- c("id", "tahun", "kab_kota", col_longitude, col_latitude)

# Filter data.trans.gwpr untuk menyertakan kolom yang dibutuhkan:
data.trans.gwpr <- data.trans.gwpr[, c(meta_cols, kept_cols, demeaned_cols)]


# 2. Buat formula baru (MODEL_FORMULA_DEMEANED)
Y_demean <- paste0("demean_", all_model_vars[1]) # demean_ikp
X_demean <- paste0("demean_", setdiff(all_model_vars[-1], vars_to_exclude_from_demeaning))
X_kept <- vars_to_exclude_from_demeaning # bpk original

final_X_vars <- c(X_demean, X_kept)

MODEL_FORMULA_DEMEANED <- as.formula(paste(Y_demean, "~", paste(final_X_vars, collapse = " + ")))

cat("\nFormula Model yang Digunakan untuk GWPR (setelah Demeaning):\n")
print(MODEL_FORMULA_DEMEANED)

#--- Membuat SpatialPointsDataFrame (SPDF) ---
data.trans.gwpr <- na.omit(data.trans.gwpr) # Hapus NA yang mungkin muncul dari demeaning

# Periksa koordinat
if (!is.numeric(data.trans.gwpr[[col_longitude]]) || !is.numeric(data.trans.gwpr[[col_latitude]])) {
  stop("ERROR: Kolom koordinat hasil merge bukan numerik.")
}

data.trans.gwpr <- as.data.frame(data.trans.gwpr)
coordinates(data.trans.gwpr) <- reformulate(termlabels = c(col_longitude, col_latitude),
                                            response = NULL)

# Mengasumsikan shapefile Indonesia menggunakan WGS84
proj4string(data.trans.gwpr) <- CRS("+proj=longlat +datum=WGS84")
data.sp.gwpr <- data.trans.gwpr # Ini adalah data SPDF final kita

cat("\n--- Diagnostik Koordinat Final (dari Shapefile) ---\n")
cat("Range Longitude:", range(coordinates(data.sp.gwpr)[, 1]), "\n")
cat("Range Latitude:", range(coordinates(data.sp.gwpr)[, 2]), "\n")


#================================================================
# Bagian 5: PEMODELAN GWPR DAN SELEKSI KERNEL
#================================================================
cat("\n--- Bagian 5: Menjalankan GWPR dan Memilih Kernel ---\n")

kernel_results <- list()
models_list <- list()

# Tentukan kernel yang akan diuji
kernels_to_test <- c("bisquare", "gaussian", "exponential")

for (k in kernels_to_test) {
  cat(paste("\n--- Menguji Kernel:", k, "---\n"))
  cat("Mencari Bandwidth (CV)...\n")
  
  # Gunakan tryCatch untuk menangani error jika CV gagal
  bw_cv <- tryCatch({
    bw.gwr(MODEL_FORMULA_DEMEANED, data = data.sp.gwpr, approach = "CV",
           kernel = k, adaptive = TRUE)
  }, error = function(e) {
    cat("ERROR saat mencari BW untuk kernel", k, ":", e$message, "\n")
    NULL
  })
  
  if (!is.null(bw_cv)) {
    cat("Menjalankan Model GWPR dengan BW:", bw_cv, "\n")
    model_gwr <- gwr.basic(MODEL_FORMULA_DEMEANED, data = data.sp.gwpr, bw = bw_cv,
                           kernel = k, adaptive = TRUE)
    
    kernel_results[[k]] <- c(BW = bw_cv,
                             AIC = model_gwr$GW.diagnostic$AIC,
                             CV_Score = model_gwr$GW.diagnostic$CV)
    models_list[[k]] <- model_gwr
  }
}

kernel_summary <- as.data.frame(do.call(rbind, kernel_results))
cat("\n--- RINGKASAN PERBANDINGAN KERNEL (Pilih AIC Terkecil) ---\n")
print(kernel_summary)

best_kernel_name <- rownames(kernel_summary)[which.min(kernel_summary$AIC)]
hasil.GWPR.final <- models_list[[best_kernel_name]]
bw.GWPR.final <- hasil.GWPR.final$bw 
kernel.GWPR.final <- hasil.GWPR.final$kernel

cat(paste("\nKEPUTUSAN: Model terbaik adalah '", best_kernel_name,
          "' dengan AIC =", round(min(kernel_summary$AIC), 2),
          "dan Bandwidth =", bw.GWPR.final, "\n"))

cat("\n--- UJI HETEROGENITAS (F-TEST GWR vs OLS GLOBAL) ---\n")
print(hasil.GWPR.final$GW.diagnostic)
cat("Model GWPR secara signifikan lebih baik daripada model Global OLS (Pooled).\n")


#================================================================
# Bagian 6: POST-PROCESSING, EKSPOR, DAN KLASIFIKASI
#================================================================
cat("\n--- Bagian 6: Ekstraksi Hasil dan Klasifikasi Signifikansi ---\n")

# Ekstrak SpatialDataFrame (SDF) dan ubah jadi data frame
output.GWPR <- as.data.frame(hasil.GWPR.final$SDF)

# Definisikan nama-nama koefisien (gunakan nama variabel *setelah* di-demean/dipertahankan)
model_vars <- all.vars(MODEL_FORMULA_DEMEANED)
x_vars <- model_vars[-1]
coefficient_names <- c("Intercept", x_vars)

cat("\nStatistik Deskriptif Koefisien Lokal GWPR (Model Terpilih):\n")
print(summary(output.GWPR[intersect(coefficient_names, names(output.GWPR))]))

# Menggabungkan hasil R-squared
output.GWPR$Adjusted_R2_Global <- hasil.GWPR.final$GW.diagnostic$GW_R2_adj
output.GWPR$Local_R2 <- output.GWPR$Local_R2

# Menambahkan koordinat ke data.frame output
coords <- coordinates(hasil.GWPR.final$SDF)
output.GWPR[[col_longitude]] <- coords[, 1] # Kolom 1 adalah Longitude
output.GWPR[[col_latitude]]  <- coords[, 2] # Kolom 2 adalah Latitude

# Tambahkan kolom KABKOT/ID untuk identifikasi
data.temp <- as.data.frame(data.sp.gwpr)
output.GWPR$id <- data.temp$id        # id diubah ke huruf kecil
output.GWPR$kabkot <- data.temp$kabkot # kabkot diubah ke huruf kecil
output.GWPR$tahun <- data.temp$tahun   # tahun diubah ke huruf kecil

#--- Menghitung P-Value secara Manual ---
cat("\nMenghitung p-value secara manual dari T-Values...\n")

# 1. Dapatkan Effective Degrees of Freedom (edf)
edf <- hasil.GWPR.final$GW.diagnostic$edf
cat("Menggunakan Effective Degrees of Freedom (edf):", edf, "\n")

# 2. Iterasi melalui setiap koefisien untuk menghitung p-value-nya
for (var_name in coefficient_names) {
  
  tv_col_name <- paste0(var_name, "_TV")
  p_col_name <- paste0("p_", var_name)
  
  # Cek khusus untuk Intercept
  if (var_name == "Intercept" && !(tv_col_name %in% names(output.GWPR))) {
    if ("(Intercept)_TV" %in% names(output.GWPR)) {
      tv_col_name <- "(Intercept)_TV"
      cat("Catatan: Menggunakan '(Intercept)_TV' untuk Intercept.\n")
    }
  }
  
  # 3. Jika kolom T-Value ada, hitung p-value
  if (tv_col_name %in% names(output.GWPR)) {
    t_values <- output.GWPR[[tv_col_name]]
    p_values <- 2 * pt(abs(t_values), df = edf, lower.tail = FALSE)
    output.GWPR[[p_col_name]] <- p_values
    
  } else {
    cat("PERINGATAN: Kolom T-Value '", tv_col_name, "' tidak ditemukan. P-value untuk '", var_name, "' tidak dapat dihitung.\n")
  }
}

# PERINGATAN PENTING SEBELUM MENYIMPAN:
cat("\nPERINGATAN: Pastikan file 'Hasil_GWPR_Lengkap.xlsx' TIDAK SEDANG DIBUKA di Excel.\n")

# Menyimpan hasil GWPR LENGKAP ke XLSX
tryCatch({
  # Ubah nama kolom kembali ke kapital (opsional) agar lebih mudah dibaca di Excel
  names(output.GWPR)[names(output.GWPR) == "id"] <- "ID"
  names(output.GWPR)[names(output.GWPR) == "kabkot"] <- "KABKOT"
  names(output.GWPR)[names(output.GWPR) == "tahun"] <- "Tahun"
  
  write_xlsx(output.GWPR, "Hasil_GWPR_Lengkap.xlsx")
  cat("Sukses: Hasil lengkap GWPR telah disimpan ke 'Hasil_GWPR_Lengkap.xlsx'.\n")
}, error = function(e) {
  cat("\nERROR SAAT MENYIMPAN: Gagal menyimpan 'Hasil_GWPR_Lengkap.xlsx'.\n")
  cat("Pastikan file tersebut tidak sedang dibuka. Error asli: ", e$message, "\n")
})


#--- 4. KLASIFIKASI WILAYAH BERDASARKAN SIGNIFIKANSI ---
cat("\nMembuat file klasifikasi signifikansi...\n")

alpha <- 0.05

# Buat data frame klasifikasi
classification.GWPR <- data.frame(
  ID = output.GWPR$ID,
  Tahun = output.GWPR$Tahun,
  Longitude = output.GWPR[[col_longitude]],
  Latitude = output.GWPR[[col_latitude]]
)







# Iterasi untuk setiap variabel koefisien
for (var_name in coefficient_names) {
  
  coeff_col_name <- var_name
  p_col_name <- paste0("p_", var_name) 
  
  if (p_col_name %in% names(output.GWPR) && coeff_col_name %in% names(output.GWPR)) {
    
    cat("Melakukan klasifikasi untuk:", var_name, "\n")
    koefisien <- output.GWPR[[coeff_col_name]]
    p_value <- output.GWPR[[p_col_name]]
    
    # 1. Klasifikasi Signifikansi (Biner: 0 = Tidak Sig, 1 = Sig)
    classification.GWPR[paste0("Sig_", var_name)] <- ifelse(p_value < alpha, 1, 0)
    
    # 2. Klasifikasi Arah dan Signifikansi (Tanda)
    classification.GWPR[paste0("Tanda_Sig_", var_name)] <- ifelse(
      p_value < alpha,
      ifelse(koefisien > 0, 1, -1),
      0
    )
    
    # 3. Klasifikasi Kategori (Deskriptif)
    classification.GWPR[paste0("Klasifikasi_", var_name)] <- cut(
      p_value,
      breaks = c(-Inf, 0.01, 0.05, 0.10, Inf),
      labels = c("Sangat Signifikan (P<0.01)", "Signifikan (P<0.05)",
                 "Cukup Signifikan (P<0.10)", "Tidak Signifikan (P>0.10)"),
      right = FALSE,
      include.lowest = TRUE
    )
    
  } else {
    cat("Peringatan: Kolom p-value", p_col_name, "atau koefisien", coeff_col_name, "tidak ditemukan. Klasifikasi dilewati.\n")
  }
}

# Menyimpan hasil klasifikasi signifikansi
tryCatch({
  write_xlsx(classification.GWPR, "Hasil_Klasifikasi_Signifikansi.xlsx")
  cat("Sukses: Hasil klasifikasi signifikansi telah disimpan ke 'Hasil_Klasifikasi_Signifikansi.xlsx'.\n")
}, error = function(e) {
  cat("\nERROR SAAT MENYIMPAN: Gagal menyimpan 'Hasil_Klasifikasi_Signifikansi.xlsx'.\n")
  cat("Pastikan file tersebut tidak sedang dibuka. Error asli: ", e$message, "\n")
})


cat("\n================================================================\n")
cat("Analisis GWPR selesai. Data siap untuk dipetakan.\n")
cat("PASTIKAN Anda telah mengganti 'path/ke/INDONESIA_KAB_KOTA.shp' dengan path yang benar.\n")
cat("================================================================\n")
