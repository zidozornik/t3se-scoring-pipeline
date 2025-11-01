#!/usr/bin/env Rscript

# run_scoring_combined.R
# Scoring pipeline (Swiss-model, AlphaFold, HHpred, genomic context, interactions)
# Usage:
#   Rscript run_scoring_combined.R input.xlsx output.xlsx
# If no args provided it uses: data/sample_input.xlsx -> results/scoring_combined.xlsx
#
# Requirements: readxl, dplyr, stringr, purrr, janitor, openxlsx
#
# Authors:
# Date:

# ---- parse CLI args / defaults ----
args <- commandArgs(trailingOnly = TRUE)
if (length(args) >= 1) {
  in_file <- args[1]
} else {
  in_file <- "data/sample_input.xlsx"
}
if (length(args) >= 2) {
  out_file <- args[2]
} else {
  out_file <- "results/scoring_combined.xlsx"
}
sheet_to_read <- 1

cat("Input file:", in_file, "\n")
cat("Output file:", out_file, "\n")

# ---- packages: install if missing, then load ----
pkgs <- c("readxl","dplyr","stringr","purrr","janitor","openxlsx")
for(p in pkgs) if(!requireNamespace(p, quietly=TRUE)) install.packages(p, repos = "https://cloud.r-project.org")
suppressPackageStartupMessages({
  library(readxl); library(dplyr); library(stringr); library(purrr); library(janitor); library(openxlsx)
})

# ---- override cols (optional) ----
override_cols <- list(
  size = NULL,
  distance = NULL,
  swiss_identity = NULL,
  swiss_gmqe = NULL,
  swiss_template = NULL,
  alpha_identity = NULL,
  alpha_plddt = NULL,
  alpha_template = NULL,
  hh_prob = NULL,
  hh_cov = NULL,
  hh_template = NULL
)

# ---- thresholds / keywords ----
size_min_aa <- 150; size_max_aa <- 900
eel_threshold_bp <- 7000; pai_threshold_bp <- 100000
interaction_threshold <- 0.5

# AlphaFold thresholds (pLDDT A = 70)
af_plddt_A <- 70
af_plddt_B <- 80
af_identity_rescue <- 95

# HHpred thresholds
hh_prob_primary <- 90; hh_cov_primary <- 0.4; hh_alnlen_min <- 40
hh_prob_unchar <- 95; hh_cov_unchar <- 0.6

# Keyword lists
unchar_keywords <- c("uncharacter","unknown","hypothetical","uniprot","putative","probable",
                     "predicted protein","conserved hypothetical","conserved protein",
                     "unnamed product","unnamed protein","unannotated","family protein","\\borf\\b","yp_","locus_tag","gene of unknown function")
toxin_keywords <- c("yenb","abc toxin","toxin complex","tc holotoxin","insecticidal toxin","tc toxin","yenb;")
effector_keywords <- c("effector","hop","avr","t3se","t3ss","type iii","type3","ascx","t3s")
enzyme_keywords <- c("phospholipase","pld","phospholipase d","phosphodiesterase","phosphatase","lipase","esterase",
                     "protease","peptidase","metalloprotease","acetyltransferase","adp-ribosyltransferase","deamidase",
                     "glycosyltransferase","glycosidase","glyoxylase","hydrolase","nuclease","rnase","ribonuclease")

# ---- parsing helpers ----
parse_num_vec <- function(x) {
  s <- as.character(x); s[is.na(s)] <- NA_character_
  s <- str_trim(s); s[s == ""] <- NA_character_
  s <- ifelse(!is.na(s), str_replace_all(s, ",", "."), NA_character_)
  tok <- str_extract(s, "-?\\d+\\.?\\d*")
  as.numeric(tok)
}
parse_int_vec <- function(x) { v <- parse_num_vec(x); ifelse(is.na(v), NA_integer_, as.integer(round(v))) }
parse_pct100_vec <- function(x) { v <- parse_num_vec(x); idx <- !is.na(v) & v <= 1; v[idx] <- v[idx] * 100; v }
contains_any_vec <- function(text_vec, kw_list) {
  s <- tolower(as.character(text_vec)); s[is.na(s)] <- ""
  map_lgl(s, function(tt) any(vapply(kw_list, function(kw) str_detect(tt, regex(kw, ignore_case = TRUE)), logical(1))))
}
norm_match_vec <- function(text_vec, keywords) {
  s <- tolower(as.character(text_vec)); s[is.na(s)] <- ""
  s_norm <- gsub("[^a-z0-9 ]", " ", s)
  s_norm <- gsub("\\s+", " ", s_norm)
  vapply(s_norm, function(tt) any(vapply(keywords, function(k) grepl(tolower(k), tt, fixed = TRUE), logical(1))), logical(1))
}

# ---- read input and auto-detect headers (single read) ----
if (!file.exists(in_file) && interactive()) {
  message("Default input not found. Please choose an input file to copy to 'data/sample_input.xlsx' (interactive).")
  src <- file.choose()                       # choose file interactively
  dir.create("data", showWarnings = FALSE)   # ensure folder exists
  ok <- file.copy(from = src, to = "data/sample_input.xlsx", overwrite = TRUE)
  if (isTRUE(ok)) {
    message("Copied file to data/sample_input.xlsx")
  } else {
    stop("Failed to copy file to data/sample_input.xlsx")
  }
}

# --- ADD THIS READ --- (paste these lines directly after the copy block)
if (!file.exists(in_file)) {
  stop("Input file still not found after interactive selection: ", in_file)
}
# read the Excel into orig_raw (must exist before any use of orig_raw)
orig_raw <- readxl::read_excel(in_file, sheet = sheet_to_read, col_names = TRUE)
if (!exists("orig_raw") || nrow(orig_raw) == 0) stop("Failed to read input into orig_raw or file is empty: ", in_file)
message("Read ", nrow(orig_raw), " rows from input: ", in_file)
# ----------------------
clean_names_vec <- names(orig_raw) %>% tolower() %>% str_replace_all("[\\s\\n\\r]+", " ") %>% str_replace_all("[^a-z0-9 ]", " ") %>% str_squish()

find_col <- function(tokens, override_val = NULL) {
  if(!is.null(override_val) && override_val %in% names(orig_raw)) return(override_val)
  for(i in seq_along(clean_names_vec)) {
    nm <- clean_names_vec[i]
    for(t in tokens) if(str_detect(nm, fixed(t))) return(names(orig_raw)[i])
  }
  return(NULL)
}

size_col       <- find_col(c("size","size (aa)","aa","length"), override_cols$size)
distance_col   <- find_col(c("distance","distance (bp)","bp","dist"), override_cols$distance)
swiss_id_col   <- find_col(c("identity","sequence identity","% id","%identity","highest sequence identity"), override_cols$swiss_identity)
swiss_gmqe_col <- find_col(c("gmqe"), override_cols$swiss_gmqe)
swiss_temp_col <- find_col(c("swiss model","swiss_model","swiss","template","homology","homolog"), override_cols$swiss_template)
alpha_id_col   <- find_col(c("identity (alphafold)","identity alphafold","identity (alpha)","identity"), override_cols$alpha_identity)
alpha_plddt_col<- find_col(c("plddt","average plddt","avg plddt","average  plddt"), override_cols$alpha_plddt)
alpha_temp_col <- find_col(c("template (alphafold)","alpha_template","alphafold","alpha template","template (alpha)","template"), override_cols$alpha_template)
hh_prob_col    <- find_col(c("probability (hhpred)","probability","prob (hhpred)","probability (hhpred domains)"), override_cols$hh_prob)
hh_cov_col     <- find_col(c("coverage (hhpred)","coverage","cov (hhpred)","coverage (hhpred domains)"), override_cols$hh_cov)
hh_temp_col    <- find_col(c("name (hhpred)","hhpred","hh template","hhpred domains","name (hhpred domains)","domain"), override_cols$hh_template)

cat("Detected columns:\n")
print(list(size_col=size_col, distance_col=distance_col, swiss_identity=swiss_id_col, swiss_gmqe=swiss_gmqe_col, swiss_template=swiss_temp_col, alpha_identity=alpha_id_col, alpha_plddt=alpha_plddt_col, alpha_template=alpha_temp_col, hh_prob=hh_prob_col, hh_cov=hh_cov_col, hh_template=hh_temp_col))
cat("If anything is NULL or wrong, set override_cols[...] near top and re-run.\n\n")

# ---- prepare df (keep orig for output) ----
df <- orig_raw %>% mutate(.row = row_number())
n <- nrow(df)

# size & distance
df$size_aa <- if(!is.null(size_col)) parse_int_vec(df[[size_col]]) else NA_integer_
df$distance_bp <- if(!is.null(distance_col)) parse_int_vec(df[[distance_col]]) else NA_integer_

df <- df %>%
  mutate(Size_bin = ifelse(!is.na(size_aa) & size_aa >= size_min_aa & size_aa <= size_max_aa, 1L, 0L),
         EEL = ifelse(!is.na(distance_bp) & distance_bp <= eel_threshold_bp, 1L, 0L),
         PAI = ifelse(!is.na(distance_bp) & distance_bp <= pai_threshold_bp, 1L, 0L))

# Swiss-model parsing & scoring
df$swiss_identity_pct <- if(!is.null(swiss_id_col)) parse_pct100_vec(df[[swiss_id_col]]) else NA_real_
df$swiss_gmqe <- if(!is.null(swiss_gmqe_col)) parse_num_vec(df[[swiss_gmqe_col]]) else NA_real_
df$swiss_template <- if(!is.null(swiss_temp_col)) as.character(df[[swiss_temp_col]]) else NA_character_
df$swiss_unchar <- contains_any_vec(df$swiss_template, unchar_keywords)

df$Swiss_bin <- purrr::pmap_int(list(df$swiss_identity_pct, df$swiss_gmqe, df$swiss_unchar), function(idv,g,unch) {
  if(!is.na(idv) && idv > 95) return(0L)
  sc <- 1L
  if(!is.na(idv) && idv >= 80 && !is.na(g) && g >= 0.70 && !isTRUE(unch)) sc <- 0L
  if(!is.na(idv) && idv >= 80 && !is.na(g) && g >= 0.80 && isTRUE(unch)) sc <- 0L
  sc
})

# AlphaFold parsing
df$af_identity_pct <- if(!is.null(alpha_id_col)) parse_pct100_vec(df[[alpha_id_col]]) else df$swiss_identity_pct
df$af_plddt <- if(!is.null(alpha_plddt_col)) parse_num_vec(df[[alpha_plddt_col]]) else NA_real_
df$af_template <- if(!is.null(alpha_temp_col)) as.character(df[[alpha_temp_col]]) else NA_character_
df$af_unchar <- contains_any_vec(df$af_template, unchar_keywords)

# HHpred parsing
df$hh_prob <- if(!is.null(hh_prob_col)) parse_num_vec(df[[hh_prob_col]]) else NA_real_
if(!is.null(hh_cov_col)) {
  tmpcov <- parse_num_vec(df[[hh_cov_col]])
  tmpcov <- ifelse(!is.na(tmpcov) & tmpcov > 1, tmpcov/100, tmpcov)
  df$hh_cov <- tmpcov
} else df$hh_cov <- NA_real_
df$hh_template <- if(!is.null(hh_temp_col)) as.character(df[[hh_temp_col]]) else NA_character_
df$hh_alnlen <- ifelse(!is.na(df$hh_cov) & !is.na(df$size_aa), round(df$hh_cov * df$size_aa), NA_integer_)

# Build combined template (AF | Swiss | HHpred) for keyword checks
df$combined_template <- paste0(
  ifelse(is.na(df$af_template),"", df$af_template), " ",
  ifelse(is.na(df$swiss_template),"", df$swiss_template), " ",
  ifelse(is.na(df$hh_template),"", df$hh_template)
)

# detect keywords in combined template
has_enzyme_any <- norm_match_vec(df$combined_template, enzyme_keywords)
has_effector_any <- norm_match_vec(df$combined_template, effector_keywords)
has_toxin_any <- norm_match_vec(df$combined_template, toxin_keywords)

# AlphaFold scoring using combined-template & pLDDT >= 70
df$AlphaFold_bin <- mapply(function(idv, plddt, enz, eff, tox, unch) {
  if(!is.na(tox) && tox) return(0L)
  if(!is.na(plddt) && plddt >= af_plddt_A && (isTRUE(enz) || isTRUE(eff))) return(1L)
  if(!is.na(plddt) && plddt >= af_plddt_B && !is.na(idv) && idv >= 80 && isTRUE(unch)) return(1L)
  if(!is.na(idv) && idv > af_identity_rescue && !is.na(plddt) && plddt >= af_plddt_A) return(1L)
  return(0L)
}, idv = df$af_identity_pct, plddt = df$af_plddt, enz = has_enzyme_any, eff = has_effector_any, tox = has_toxin_any, unch = df$af_unchar, SIMPLIFY = TRUE)

# HHpred scoring
df$hh_is_toxin <- norm_match_vec(df$hh_template, toxin_keywords)
df$hh_is_unchar <- contains_any_vec(df$hh_template, unchar_keywords)

df$HHpred_bin <- purrr::pmap_int(list(df$hh_prob, df$hh_cov, df$hh_alnlen, df$hh_is_toxin, df$hh_is_unchar),
                                 function(prob,cov,aln,is_tox,is_unchar) {
                                   if(!is.na(is_tox) && is_tox) return(0L)
                                   if(!is.na(prob) && prob >= hh_prob_primary && !is.na(cov) && cov >= hh_cov_primary && !is.na(aln) && aln >= hh_alnlen_min) return(1L)
                                   if(!is.na(is_unchar) && is_unchar && !is.na(prob) && prob >= hh_prob_unchar && !is.na(cov) && cov >= hh_cov_unchar) return(1L)
                                   return(0L)
                                 })
# enforce aligned length for characterized hits
df$HHpred_bin <- ifelse(df$HHpred_bin == 1 & is.na(df$hh_alnlen) & df$hh_is_unchar == FALSE, 0L, df$HHpred_bin)

# Interaction detection: try known names then heuristic
orig_names_norm <- tolower(names(orig_raw)) %>% str_replace_all("[\\s\\n\\r]+"," ") %>% str_replace_all("[^a-z0-9 ]"," ") %>% str_squish()
expected_targets <- c("exo70b1","mkk2","mkk5","rin4","bak1")
interaction_cols <- c()
for(tok in expected_targets) {
  idx <- which(str_detect(orig_names_norm, fixed(tok)))
  if(length(idx) > 0) interaction_cols <- c(interaction_cols, names(orig_raw)[idx[1]])
}
used_cols <- c(size_col, distance_col, swiss_id_col, swiss_gmqe_col, swiss_temp_col, alpha_id_col, alpha_plddt_col, alpha_temp_col, hh_prob_col, hh_cov_col, hh_temp_col)
used_cols <- used_cols[!sapply(used_cols, is.null)]
if(length(interaction_cols) == 0) {
  candidates <- setdiff(names(orig_raw), used_cols)
  for(c in candidates) {
    nums <- parse_num_vec(orig_raw[[c]])
    if(mean(!is.na(nums)) >= 0.4) {
      nn <- nums[!is.na(nums)]
      if(length(nn)>0) {
        frac01 <- mean(nn >= 0 & nn <= 1)
        frac0100 <- mean(nn >= 0 & nn <= 100)
        if(frac01 >= 0.4 || frac0100 >= 0.6) interaction_cols <- c(interaction_cols, c)
      }
    }
  }
}
interaction_cols <- unique(interaction_cols)
cat("Interaction columns detected:\n"); print(interaction_cols)

# compute binary interactions and sum
for(col in interaction_cols) {
  numcol <- paste0(col, "_num"); binc <- paste0(col, "_bin")
  df[[numcol]] <- parse_num_vec(orig_raw[[col]])
  df[[numcol]] <- ifelse(!is.na(df[[numcol]]) & df[[numcol]] > 1, df[[numcol]] / 100, df[[numcol]])
  df[[binc]] <- ifelse(!is.na(df[[numcol]]) & df[[numcol]] > interaction_threshold, 1L, 0L)
}
bin_cols <- paste0(interaction_cols, "_bin")
df$Interaction_score <- if(length(bin_cols) > 0) rowSums(df[, bin_cols, drop = FALSE], na.rm = TRUE) else 0L

# final total score
score_cols <- c("Swiss_bin","AlphaFold_bin","HHpred_bin","Size_bin","EEL","PAI","Interaction_score")
df$Total_score <- rowSums(df[, score_cols], na.rm = TRUE)

# prepare output: original + score cols
out_df <- bind_cols(orig_raw, df %>% select(all_of(score_cols), Total_score))

# write Excel
dir.create(dirname(out_file), recursive = TRUE, showWarnings = FALSE)

wb <- createWorkbook()
addWorksheet(wb, "scored")
writeData(wb, "scored", out_df)
saveWorkbook(wb, out_file, overwrite = TRUE)

cat("Wrote output to:", out_file, "\n\n")

wb <- createWorkbook(); addWorksheet(wb, "scored"); writeData(wb, "scored", out_df); saveWorkbook(wb, out_file, overwrite = TRUE)

# summary
cat("Summary positives per feature:\n"); print(sapply(score_cols, function(cn) sum(out_df[[cn]] == 1, na.rm = TRUE)))
cat("Total rows:", nrow(out_df), "\n")