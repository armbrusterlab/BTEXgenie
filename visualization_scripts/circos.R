#!/usr/bin/env Rscript

# Circos visualization of pathway hits and GC skew.

suppressPackageStartupMessages({
  library(optparse)
  library(readr)
  library(dplyr)
  library(stringr)
  library(tidyr)
  library(purrr)
  library(circlize)
  library(GenomicRanges)
})

ensure_dir <- function(p) {
  dir.create(p, recursive = TRUE, showWarnings = FALSE)
}

hex_to_rgb_str <- function(hexcode) {
  h <- gsub("^#", "", trimws(hexcode))
  if (nchar(h) != 6) stop(sprintf("Bad hex color: %s", hexcode), call. = FALSE)
  vals <- c(
    strtoi(substr(h, 1, 2), 16L),
    strtoi(substr(h, 3, 4), 16L),
    strtoi(substr(h, 5, 6), 16L)
  )
  paste(vals, collapse = ",")
}

sanitize_color_name <- function(token) {
  token |>
    tolower() |>
    gsub("[^a-z0-9]+", "_", x = _)
}

normalize_color_token <- function(token) {
  if (length(token) != 1 || is.na(token)) {
    stop("Color token is missing.", call. = FALSE)
  }
  t <- trimws(token)
  if (!nzchar(t)) {
    stop("Color token is blank.", call. = FALSE)
  }
  if (grepl(",", t, fixed = TRUE)) {
    parts <- trimws(strsplit(t, ",", fixed = TRUE)[[1]])
    if (length(parts) != 3) stop(sprintf("Unsupported color token: %s", token), call. = FALSE)
    vals <- suppressWarnings(as.integer(parts))
    if (any(is.na(vals)) || any(vals < 0 | vals > 255)) {
      stop(sprintf("Unsupported color token: %s", token), call. = FALSE)
    }
    return(paste(vals, collapse = ","))
  }
  if (startsWith(t, "#")) return(hex_to_rgb_str(t))
  stop(sprintf("Unsupported color token: %s", token), call. = FALSE)
}

rgb_str_to_hex <- function(rgb_str) {
  parts <- trimws(strsplit(rgb_str, ",", fixed = TRUE)[[1]])
  vals <- as.integer(parts)
  if (length(vals) != 3 || any(is.na(vals))) stop(sprintf("Bad rgb string: %s", rgb_str), call. = FALSE)
  grDevices::rgb(vals[1], vals[2], vals[3], maxColorValue = 255)
}

mix_with_white <- function(hex_color, weight = 0.5) {
  weight <- max(0, min(1, weight))
  rgb_mat <- grDevices::col2rgb(hex_color)
  vals <- as.numeric(rgb_mat[, 1])
  mixed <- round(255 * (1 - weight) + vals * weight)
  grDevices::rgb(mixed[1], mixed[2], mixed[3], maxColorValue = 255)
}

make_density_shades <- function(base_hex, n_levels = 4) {
  if (n_levels < 1) stop("n_levels must be >= 1", call. = FALSE)
  if (n_levels == 1) return(base_hex)
  weights <- seq(0.30, 1.00, length.out = n_levels)
  vapply(weights, function(w) mix_with_white(base_hex, w), character(1))
}

infer_default_pathway_map <- function() {
  candidate1 <- file.path("/home/juneq/BTEXgenie/btexhmm/data/pathway_map.tsv")
  if (file.exists(candidate1)) return(candidate1)
  stop(
    sprintf("[ERROR] Default pathway map not found. Tried:\n  %s\nUse --pathway-map to provide it explicitly.", candidate1),
    call. = FALSE
  )
}

infer_default_ko_category_map <- function() {
  candidate1 <- file.path("/home/juneq/BTEX-HMMs/btexhmm/data/A_ko_category_map.tsv")
  if (file.exists(candidate1)) return(candidate1)
  stop(
    sprintf("[ERROR] Default KO category map not found. Tried:\n  %s\nUse --ko-category-map to provide it explicitly.", candidate1),
    call. = FALSE
  )
}

read_fasta_lengths <- function(dna_path) {
  lines <- readLines(dna_path, warn = FALSE)
  headers <- which(startsWith(lines, ">"))
  if (length(headers) == 0) stop(sprintf("[ERROR] No FASTA headers found in %s", dna_path), call. = FALSE)

  contigs <- vector("list", length(headers))
  for (i in seq_along(headers)) {
    start <- headers[i]
    end <- if (i < length(headers)) headers[i + 1] - 1 else length(lines)
    header <- sub("^>", "", lines[start])
    contig_name <- strsplit(trimws(header), "\\s+")[[1]][1]
    seq_lines <- if (end > start) lines[(start + 1):end] else character(0)
    seq_len <- sum(nchar(gsub("\\s+", "", seq_lines)))
    contigs[[i]] <- data.frame(contig = contig_name, length = seq_len, stringsAsFactors = FALSE)
  }
  bind_rows(contigs)
}

read_fasta_sequences <- function(dna_path) {
  lines <- readLines(dna_path, warn = FALSE)
  headers <- which(startsWith(lines, ">"))
  if (length(headers) == 0) stop(sprintf("[ERROR] No FASTA headers found in %s", dna_path), call. = FALSE)

  seqs <- vector("list", length(headers))
  for (i in seq_along(headers)) {
    start <- headers[i]
    end <- if (i < length(headers)) headers[i + 1] - 1 else length(lines)
    header <- sub("^>", "", lines[start])
    contig_name <- strsplit(trimws(header), "\\s+")[[1]][1]
    seq_lines <- if (end > start) lines[(start + 1):end] else character(0)
    seq <- paste0(gsub("\\s+", "", seq_lines), collapse = "")
    seqs[[i]] <- tibble(contig = contig_name, sequence = seq)
  }
  bind_rows(seqs)
}

generate_contig_lengths_table_from_dna <- function(dna_path, genome, out_path) {
  if (!file.exists(dna_path)) stop(sprintf("[ERROR] DNA FASTA not found: %s", dna_path), call. = FALSE)
  entries <- read_fasta_lengths(dna_path)
  if (nrow(entries) == 0) stop(sprintf("[ERROR] No contig lengths gathered from %s", dna_path), call. = FALSE)
  out_df <- entries |>
    mutate(sample = genome) |>
    select(sample, contig, length)
  ensure_dir(dirname(out_path))
  write_tsv(out_df, out_path)
  message(sprintf("[log] wrote %d contig entries to %s", nrow(out_df), out_path))
  out_path
}

generate_contig_lengths_for_sample <- function(opt, output_dir) {
  default_out <- file.path(output_dir, "contig_length.tsv")
  ensure_dir(output_dir)
  if (is.null(opt$genome) || !nzchar(opt$genome)) {
    stop("[ERROR] Provide -g/--genome.", call. = FALSE)
  }
  message(sprintf("[log] Auto-generating contig lengths -> %s", default_out))
  generate_contig_lengths_table_from_dna(opt$genome, opt$sample, default_out)
}

read_contig_lengths <- function(path, genome) {
  read_tsv(path, show_col_types = FALSE) |>
    filter(sample == genome) |>
    mutate(length = as.integer(length)) |>
    select(contig, length)
}

normalize_contig_id <- function(token) {
  if (is.null(token) || !nzchar(token)) return("")
  name <- token |>
    trimws() |>
    sub("^>", "", x = _) |>
    trimws()
  if (grepl(" #", name, fixed = TRUE)) name <- strsplit(name, " #", fixed = TRUE)[[1]][1]
  if (grepl("|", name, fixed = TRUE)) {
    parts <- strsplit(name, "|", fixed = TRUE)[[1]]
    name <- parts[length(parts)]
  }
  trimws(name)
}

parse_hit_headers <- function(sample, hmm, headers) {
  if (is.na(headers) || !nzchar(headers)) return(tibble())
  raw_parts <- strsplit(headers, "|", fixed = TRUE)[[1]]
  out <- vector("list", length(raw_parts))
  keep <- 0

  for (raw in raw_parts) {
    if (!nzchar(trimws(raw))) next
    hash_fields <- trimws(strsplit(raw, "#", fixed = TRUE)[[1]])
    left <- hash_fields[1]
    contig <- normalize_contig_id(left)
    if (grepl(" #", raw, fixed = TRUE)) {
      contig <- sub("_\\d+$", "", contig)
    }
    if (!nzchar(contig)) next

    if (length(hash_fields) >= 4) {
      nums <- suppressWarnings(as.integer(hash_fields[2:3]))
      strand_num <- suppressWarnings(as.integer(hash_fields[4]))
    } else {
      nums <- as.integer(stringr::str_extract_all(raw, "\\b(\\d+)\\b")[[1]])
      strand_num <- if (length(nums) >= 2 && nums[1] > nums[2]) -1L else 1L
    }
    if (length(nums) < 2 || any(is.na(nums[1:2]))) next
    if (is.na(strand_num) || !strand_num %in% c(-1L, 1L)) {
      strand_num <- if (nums[1] > nums[2]) -1L else 1L
    }

    a <- nums[1]
    b <- nums[2]
    keep <- keep + 1
    out[[keep]] <- tibble(
      sample = sample,
      hmm = hmm,
      contig = contig,
      start = min(a, b),
      end = max(a, b),
      strand = strand_num
    )
  }

  if (keep == 0) return(tibble())
  bind_rows(out[seq_len(keep)])
}

read_hmmscan_features <- function(path, genome) {
  hmmscan_df <- read_csv(path, show_col_types = FALSE)

  # Current BTEX-HMM output: one row per best protein hit.
  if (all(c("sample", "hmm", "hit_header") %in% colnames(hmmscan_df))) {
    selected_rows <- hmmscan_df |>
      filter(sample == genome, !is.na(hit_header), nzchar(hit_header))

    if (nrow(selected_rows) == 0) return(tibble())

    features_list <- pmap(
      list(selected_rows$sample, selected_rows$hmm, selected_rows$hit_header),
      parse_hit_headers
    )
    return(bind_rows(features_list))
  }

  # Legacy BTEX-HMM output: one row per sample/HMM with aggregated hit headers.
  required_cols <- c("sample", "hmm", "hits", "hit_headers")
  missing <- setdiff(required_cols, colnames(hmmscan_df))
  if (length(missing) > 0) {
    stop(sprintf("[ERROR] hmmscan CSV missing required columns. Expected either sample,hmm,hit_header or sample,hmm,hits,hit_headers. Missing: %s", paste(sort(missing), collapse = ", ")), call. = FALSE)
  }

  selected_rows <- hmmscan_df |>
    filter(sample == genome, as.integer(hits) >= 1)

  if (nrow(selected_rows) == 0) return(tibble())

  features_list <- pmap(
    list(selected_rows$sample, selected_rows$hmm, selected_rows$hit_headers),
    parse_hit_headers
  )
  bind_rows(features_list)
}

load_pathway_map <- function(pathway_map_path) {
  pm <- read_tsv(pathway_map_path, show_col_types = FALSE)
  required <- c("hmm", "pathway", "color")
  missing <- setdiff(required, colnames(pm))
  if (length(missing) > 0) {
    stop(sprintf("[ERROR] Pathway map missing required columns: %s", paste(sort(missing), collapse = ", ")), call. = FALSE)
  }

  pm <- pm |>
    mutate(
      hmm = trimws(hmm),
      pathway = trimws(pathway),
      color = trimws(color)
    )

  bad <- pm |>
    filter(!nzchar(hmm) | !nzchar(pathway) | !nzchar(color))
  if (nrow(bad) > 0) {
    stop("[ERROR] Pathway map row must include hmm, pathway, and color values.", call. = FALSE)
  }
  pm
}

build_color_map <- function(hmms_with_hits, pathway_df) {
  pathway_map <- pathway_df |>
    select(hmm, pathway) |>
    distinct()
  pathway_colors <- pathway_df |>
    select(pathway, color) |>
    distinct() |>
    mutate(rgb = map_chr(color, normalize_color_token))

  final <- lapply(hmms_with_hits, function(hmm) {
    pathway_row <- pathway_map |> filter(tolower(hmm) == tolower(!!hmm))
    if (nrow(pathway_row) == 0) {
      stop(sprintf("[ERROR] No pathway mapping found for HMM: %s", hmm), call. = FALSE)
    }
    pathway <- pathway_row$pathway[1]
    color_row <- pathway_colors |> filter(tolower(pathway) == tolower(!!pathway))
    if (nrow(color_row) == 0) {
      stop(sprintf("[ERROR] No color defined in pathway_map.tsv for pathway '%s' (HMM: %s)", pathway, hmm), call. = FALSE)
    }
    tibble(
      hmm = hmm,
      color_name = sanitize_color_name(pathway),
      rgb = color_row$rgb[1],
      color_hex = rgb_str_to_hex(color_row$rgb[1]),
      pathway = pathway
    )
  })
  bind_rows(final)
}

compute_gc_skew_windows <- function(dna_path, contigs_to_keep = NULL, window_size = 5000L, step_size = 5000L) {
  seq_df <- read_fasta_sequences(dna_path)
  if (!is.null(contigs_to_keep)) {
    seq_df <- seq_df |>
      filter(contig %in% contigs_to_keep)
  }

  out <- list()
  idx <- 1L

  for (i in seq_len(nrow(seq_df))) {
    contig <- seq_df$contig[i]
    seq <- toupper(seq_df$sequence[i])
    n <- nchar(seq)
    if (n == 0) next

    starts <- seq.int(1L, n, by = step_size)
    for (s in starts) {
      e <- min(s + window_size - 1L, n)
      subseq <- substr(seq, s, e)
      chars <- strsplit(subseq, "", fixed = TRUE)[[1]]
      g_count <- sum(chars == "G")
      c_count <- sum(chars == "C")
      denom <- g_count + c_count
      skew <- if (denom == 0L) 0 else (g_count - c_count) / denom

      out[[idx]] <- tibble(
        contig = contig,
        start = s,
        end = e,
        gc_skew = skew
      )
      idx <- idx + 1L
    }
  }

  if (length(out) == 0) {
    return(tibble(contig = character(), start = integer(), end = integer(), gc_skew = numeric()))
  }
  bind_rows(out)
}

extract_contig_from_definition <- function(line) {
  m_seqhdr <- stringr::str_match(line, 'seqhdr="([^\"]+)"')
  if (!is.na(m_seqhdr[1,2])) {
    seqhdr <- trimws(m_seqhdr[1,2])
    first_token <- strsplit(seqhdr, "\\s+")[[1]][1]
    if (nzchar(first_token)) return(first_token)
  }

  m_contig <- stringr::str_match(line, 'contig:\\s*([^,;\\"]+)')
  if (!is.na(m_contig[1,2])) return(trimws(m_contig[1,2]))

  NA_character_
}

parse_location_to_range <- function(loc_text) {
  loc <- trimws(loc_text)
  strand <- 1L
  if (startsWith(loc, "complement(")) {
    strand <- -1L
    loc <- sub('^complement\\((.*)\\)$', '\\1', loc)
  }
  loc <- gsub('[<>]', '', loc)
  loc <- sub('^join\\((.*)\\)$', '\\1', loc)
  parts <- trimws(strsplit(loc, ',', fixed = TRUE)[[1]])
  starts <- integer(0)
  ends <- integer(0)
  for (part in parts) {
    mm <- stringr::str_match(part, '([0-9]+)\\.\\.([0-9]+)')
    if (!is.na(mm[1,2]) && !is.na(mm[1,3])) {
      starts <- c(starts, as.integer(mm[1,2]))
      ends <- c(ends, as.integer(mm[1,3]))
    } else {
      ms <- stringr::str_match(part, '^([0-9]+)$')
      if (!is.na(ms[1,2])) {
        starts <- c(starts, as.integer(ms[1,2]))
        ends <- c(ends, as.integer(ms[1,2]))
      }
    }
  }
  if (length(starts) == 0) return(NULL)
  list(start = min(starts), end = max(ends), strand = strand)
}

extract_qualifier_value <- function(x, key) {
  pattern <- sprintf('/%s="([^\"]+)"', key)
  m <- stringr::str_match(x, pattern)
  val <- m[, 2]
  val[is.na(val)] <- ""
  val
}

extract_prodigal_note_id <- function(x) {
  m <- stringr::str_match(x, 'ID=([^;"]+)')
  val <- m[, 2]
  val[is.na(val)] <- ""
  val
}

read_prodigal_gbk_cds <- function(gbk_path) {
  if (!file.exists(gbk_path)) stop(sprintf("[ERROR] Prodigal GBK not found: %s", gbk_path), call. = FALSE)
  lines <- readLines(gbk_path, warn = FALSE)
  if (length(lines) == 0) {
    return(tibble(contig = character(), start = integer(), end = integer(), strand = integer(), cds_id = character(), locus_tag = character(), protein_id = character(), prodigal_id = character()))
  }

  current_contig <- NA_character_
  records <- list()
  idx <- 1L
  per_contig_idx <- integer(0)
  names(per_contig_idx) <- character(0)

  i <- 1L
  while (i <= length(lines)) {
    line <- lines[i]
    if (startsWith(line, "DEFINITION")) {
      current_contig <- extract_contig_from_definition(line)
      i <- i + 1L
      next
    }

    if (grepl('^\\s{5}CDS\\s+', line)) {
      if (is.na(current_contig) || !nzchar(current_contig)) {
        i <- i + 1L
        next
      }

      loc_text <- trimws(sub('^\\s{5}CDS\\s+', '', line))
      loc_parsed <- parse_location_to_range(loc_text)
      if (is.null(loc_parsed)) {
        i <- i + 1L
        next
      }

      qualifiers <- character(0)
      j <- i + 1L
      while (j <= length(lines) && grepl('^\\s{21}/', lines[j])) {
        qualifiers <- c(qualifiers, trimws(lines[j]))
        j <- j + 1L
      }
      qual_text <- paste(qualifiers, collapse = " ")
      locus_tag <- extract_qualifier_value(qual_text, "locus_tag")
      protein_id <- extract_qualifier_value(qual_text, "protein_id")
      prodigal_id <- extract_prodigal_note_id(qual_text)

      if (!current_contig %in% names(per_contig_idx)) per_contig_idx[[current_contig]] <- 0L
      per_contig_idx[[current_contig]] <- per_contig_idx[[current_contig]] + 1L
      cds_id <- paste0(current_contig, "_", per_contig_idx[[current_contig]])

      records[[idx]] <- tibble(
        contig = current_contig,
        start = loc_parsed$start,
        end = loc_parsed$end,
        strand = loc_parsed$strand,
        cds_id = cds_id,
        locus_tag = locus_tag,
        protein_id = protein_id,
        prodigal_id = prodigal_id
      )
      idx <- idx + 1L
      i <- j
      next
    }

    i <- i + 1L
  }

  if (length(records) == 0) {
    return(tibble(contig = character(), start = integer(), end = integer(), strand = integer(), cds_id = character(), locus_tag = character(), protein_id = character(), prodigal_id = character()))
  }
  bind_rows(records)
}

read_delim_guess <- function(path) {
  first_non_comment <- readLines(path, warn = FALSE)
  first_non_comment <- first_non_comment[!startsWith(trimws(first_non_comment), "#")]
  first_non_comment <- first_non_comment[nzchar(trimws(first_non_comment))]
  if (length(first_non_comment) == 0) {
    stop(sprintf("[ERROR] No data lines found in %s", path), call. = FALSE)
  }
  line <- first_non_comment[1]
  delim <- if (grepl("\\t", line)) "\t" else if (grepl(",", line)) "," else "\t"
  suppressMessages(readr::read_delim(path, delim = delim, col_types = cols(.default = col_character()), show_col_types = FALSE))
}

normalize_colnames <- function(x) {
  x |>
    tolower() |>
    gsub("[^a-z0-9]+", "_", x = _)
}

read_kofam_assignments <- function(path) {
  if (!file.exists(path)) stop(sprintf("[ERROR] KOfam file not found: %s", path), call. = FALSE)

  lines <- readLines(path, warn = FALSE)
  data_lines <- lines[!startsWith(trimws(lines), "#")]
  data_lines <- data_lines[nzchar(trimws(data_lines))]
  if (length(data_lines) == 0) {
    return(tibble(query_id = character(), ko = character(), score = numeric(), threshold = numeric(), definition = character()))
  }

  first <- data_lines[1]

  if (startsWith(first, "*") || startsWith(first, "+") || startsWith(first, "-")) {
    split_rows <- strsplit(data_lines, "\\t")
    parsed <- lapply(split_rows, function(parts) {
      parts <- trimws(parts)
      if (length(parts) < 5) return(NULL)
      if (!grepl("^K\\d+$", parts[2])) return(NULL)
      tibble(
        query_id = parts[5],
        ko = parts[2],
        score = suppressWarnings(as.numeric(parts[3])),
        threshold = suppressWarnings(as.numeric(parts[4])),
        definition = if (length(parts) >= 6) paste(parts[6:length(parts)], collapse = " ") else ""
      )
    })
    out <- bind_rows(parsed)
    return(out |> filter(nzchar(query_id), nzchar(ko)))
  }

  df <- read_delim_guess(path)
  if (nrow(df) == 0) {
    return(tibble(query_id = character(), ko = character(), score = numeric(), threshold = numeric(), definition = character()))
  }
  names(df) <- normalize_colnames(names(df))

  query_candidates <- c("query_id", "query", "protein_id", "gene_id", "gene_name", "gene", "target", "sequence_name")
  ko_candidates <- c("ko", "k_number", "knum", "k_number_id")
  score_candidates <- c("score", "bit_score")
  threshold_candidates <- c("threshold", "score_threshold")
  def_candidates <- c("definition", "description", "function")

  query_col <- query_candidates[query_candidates %in% names(df)][1]
  ko_col <- ko_candidates[ko_candidates %in% names(df)][1]
  score_col <- score_candidates[score_candidates %in% names(df)][1]
  threshold_col <- threshold_candidates[threshold_candidates %in% names(df)][1]
  def_col <- def_candidates[def_candidates %in% names(df)][1]

  if (is.na(query_col) || is.na(ko_col)) {
    stop("[ERROR] Could not identify query ID and KO columns in the KOfam file. Use a KOfamScan detail output or a table with query_id and ko columns, or a table with gene_name and KO columns.", call. = FALSE)
  }

  tibble(
    query_id = df[[query_col]],
    ko = df[[ko_col]],
    score = if (!is.na(score_col)) suppressWarnings(as.numeric(df[[score_col]])) else NA_real_,
    threshold = if (!is.na(threshold_col)) suppressWarnings(as.numeric(df[[threshold_col]])) else NA_real_,
    definition = if (!is.na(def_col)) df[[def_col]] else ""
  ) |>
    filter(nzchar(query_id), nzchar(ko))
}

read_ko_category_map <- function(path) {
  if (!file.exists(path)) stop(sprintf("[ERROR] KO category map not found: %s", path), call. = FALSE)
  df <- read_delim_guess(path)
  names(df) <- normalize_colnames(names(df))

  ko_col <- c("ko", "k_number", "knum")[c("ko", "k_number", "knum") %in% names(df)][1]
  cat_col <- c("category", "group", "pathway", "class")[c("category", "group", "pathway", "class") %in% names(df)][1]
  color_col <- c("color", "hex", "hex_color")[c("color", "hex", "hex_color") %in% names(df)][1]

  if (is.na(ko_col) || is.na(cat_col) || is.na(color_col)) {
    stop("[ERROR] KO category map must contain KO, category, and color columns.", call. = FALSE)
  }

  out <- df |>
    transmute(
      ko = ifelse(is.na(.data[[ko_col]]), "", trimws(.data[[ko_col]])),
      category = ifelse(is.na(.data[[cat_col]]), "", trimws(.data[[cat_col]])),
      color = ifelse(is.na(.data[[color_col]]), "", trimws(.data[[color_col]]))
    ) |>
    filter(nzchar(ko), nzchar(category), nzchar(color)) |>
    distinct()

  if (nrow(out) == 0) {
    stop("[ERROR] KO category map contains no usable rows after removing blank KO/category/color entries.", call. = FALSE)
  }

  bad_colors <- which(!grepl("^#[0-9A-Fa-f]{6}$", out$color) & !grepl("^\\s*\\d+\\s*,\\s*\\d+\\s*,\\s*\\d+\\s*$", out$color))
  if (length(bad_colors) > 0) {
    idx <- bad_colors[1]
    stop(sprintf("[ERROR] Invalid color token in KO category map for KO '%s' and category '%s': %s", out$ko[idx], out$category[idx], out$color[idx]), call. = FALSE)
  }

  out |>
    mutate(
      rgb = map_chr(color, normalize_color_token),
      color_hex = map_chr(rgb, rgb_str_to_hex)
    )
}

build_cds_id_lookup <- function(cds_df) {
  bind_rows(
    cds_df |>
      transmute(query_id = cds_id, contig, start, end, strand, cds_id, locus_tag, protein_id, prodigal_id),
    cds_df |>
      filter(nzchar(locus_tag)) |>
      transmute(query_id = locus_tag, contig, start, end, strand, cds_id, locus_tag, protein_id, prodigal_id),
    cds_df |>
      filter(nzchar(protein_id)) |>
      transmute(query_id = protein_id, contig, start, end, strand, cds_id, locus_tag, protein_id, prodigal_id),
    cds_df |>
      filter(nzchar(prodigal_id)) |>
      transmute(query_id = prodigal_id, contig, start, end, strand, cds_id, locus_tag, protein_id, prodigal_id)
  ) |>
    distinct(query_id, .keep_all = TRUE)
}

join_kofam_to_cds <- function(kofam_df, cds_df, ko_cat_df) {
  if (nrow(kofam_df) == 0) {
    return(tibble(contig = character(), start = integer(), end = integer(), strand = integer(), query_id = character(), cds_id = character(), locus_tag = character(), protein_id = character(), prodigal_id = character(), ko = character(), category = character(), color_hex = character(), score = numeric(), threshold = numeric(), definition = character()))
  }

  lookup <- build_cds_id_lookup(cds_df)

  kofam_df |>
    left_join(lookup, by = "query_id") |>
    filter(!is.na(contig), !is.na(start), !is.na(end)) |>
    left_join(ko_cat_df |> select(ko, category, color_hex), by = "ko") |>
    filter(!is.na(category), !is.na(color_hex)) |>
    distinct(query_id, ko, category, contig, start, end, .keep_all = TRUE)
}

make_fixed_windows <- function(contig_lengths, window_size = 5000L) {
  out <- vector("list", nrow(contig_lengths))
  for (i in seq_len(nrow(contig_lengths))) {
    ctg <- contig_lengths$contig[i]
    len <- contig_lengths$length[i]
    starts <- seq.int(1L, len, by = window_size)
    ends <- pmin(starts + window_size - 1L, len)
    out[[i]] <- tibble(contig = ctg, start = starts, end = ends)
  }
  bind_rows(out) |>
    mutate(window_id = row_number())
}

compute_category_density_windows <- function(interval_df, contig_lengths, window_size = 5000L) {
  windows <- make_fixed_windows(contig_lengths, window_size = window_size)
  if (nrow(interval_df) == 0) {
    empty_summary <- windows |>
      mutate(category = NA_character_, category_color = NA_character_, window_value = 0, density_level = 0L, fill = "#FFFFFF")
    return(list(summary = empty_summary, matrix = tibble()))
  }

  color_col <- c("color_hex", "category_color", "color")[c("color_hex", "category_color", "color") %in% names(interval_df)][1]
  if (is.na(color_col)) {
    stop("[ERROR] KOfam interval data must contain a color column (expected one of: color_hex, category_color, color).", call. = FALSE)
  }

  interval_df <- interval_df |>
    mutate(category_color = .data[[color_col]])

  win_gr <- GenomicRanges::GRanges(
    seqnames = windows$contig,
    ranges = IRanges::IRanges(start = windows$start, end = windows$end),
    window_id = windows$window_id
  )

  cats <- sort(unique(interval_df$category))
  per_cat <- vector("list", length(cats))
  names(per_cat) <- cats

  for (cat in cats) {
    sub_df <- interval_df |>
      filter(category == cat) |>
      distinct(contig, start, end, query_id, ko, category, category_color)
    if (nrow(sub_df) == 0) next

    reg_gr <- GenomicRanges::GRanges(
      seqnames = sub_df$contig,
      ranges = IRanges::IRanges(start = sub_df$start, end = sub_df$end)
    )
    hits <- GenomicRanges::findOverlaps(win_gr, reg_gr, ignore.strand = TRUE)
    counts <- tabulate(S4Vectors::queryHits(hits), nbins = length(win_gr))
    per_cat[[cat]] <- windows |>
      mutate(category = cat, value = counts)
  }

  matrix_df <- bind_rows(per_cat) |>
    left_join(interval_df |> distinct(category, category_color), by = "category") |>
    relocate(category_color, .after = category)

  if (nrow(matrix_df) == 0) {
    empty_summary <- windows |>
      mutate(category = NA_character_, category_color = NA_character_, window_value = 0, density_level = 0L, fill = "#FFFFFF")
    return(list(summary = empty_summary, matrix = tibble()))
  }

  max_value <- max(matrix_df$value, na.rm = TRUE)
  n_levels <- 4L

  summary_df <- matrix_df |>
    arrange(contig, start, end, desc(value), category) |>
    group_by(window_id, contig, start, end) |>
    dplyr::slice(1) |>
    ungroup() |>
    mutate(
      category = if_else(value > 0, category, NA_character_),
      category_color = if_else(value > 0, category_color, NA_character_),
      window_value = value,
      density_level = if_else(
        value > 0,
        pmax(1L, ceiling(value / max(1, max_value) * n_levels)),
        0L
      )
    ) |>
    rowwise() |>
    mutate(
      fill = if (is.na(category_color) || density_level <= 0L) {
        "#FFFFFF"
      } else {
        make_density_shades(category_color, n_levels = n_levels)[density_level]
      }
    ) |>
    ungroup() |>
    select(contig, start, end, category, category_color, window_value, density_level, fill)

  matrix_df <- matrix_df |>
    mutate(
      density_level = if_else(
        value > 0,
        pmax(1L, ceiling(value / max(1, max_value) * n_levels)),
        0L
      )
    )

  list(summary = summary_df, matrix = matrix_df)
}

write_hmm_colors_tsv <- function(path, cmap) {
  if (nrow(cmap) == 0) {
    write_tsv(tibble(hmm = character(), color_name = character(), rgb = character()), path)
    return(invisible(NULL))
  }
  write_tsv(cmap |> select(hmm, color_name, rgb), path)
}

write_karyotype_tsv <- function(path, contig_lengths) {
  write_tsv(contig_lengths, path)
}

write_gene_hits_tsv <- function(path, features, cmap) {
  out <- features |>
    left_join(cmap |> select(hmm, color_name, rgb, pathway), by = "hmm")
  write_tsv(out, path)
}

write_gene_labels_tsv <- function(path, features) {
  out <- features |>
    distinct(contig, start, end, hmm, strand)
  write_tsv(out, path)
}

write_genbank_like_hits_tsv <- function(path, features, contig_lengths) {
  out <- features |>
    left_join(contig_lengths, by = c("contig" = "contig")) |>
    arrange(contig, start, end, hmm)
  write_tsv(out, path)
}

format_genbank_location <- function(start, end, strand) {
  loc <- sprintf("%d..%d", as.integer(start), as.integer(end))
  if (!is.na(strand) && as.integer(strand) == -1L) {
    return(sprintf("complement(%s)", loc))
  }
  loc
}

wrap_genbank_text <- function(prefix, text, width = 80) {
  if (is.na(text) || !nzchar(text)) return(character(0))
  available <- max(1, width - nchar(prefix))
  wrapped <- strwrap(text, width = available, exdent = 0)
  c(
    paste0(prefix, wrapped[1]),
    if (length(wrapped) > 1) paste0(strrep(" ", nchar(prefix)), wrapped[-1]) else character(0)
  )
}

format_genbank_qualifier <- function(key, value) {
  if (is.na(value) || !nzchar(value)) return(character(0))
  escaped <- gsub('"', "'", as.character(value), fixed = TRUE)
  wrap_genbank_text(sprintf('                     /%s="', key), paste0(escaped, '"'))
}

format_origin_lines <- function(seq) {
  seq <- tolower(seq)
  if (!nzchar(seq)) return("        1")

  seq_len <- nchar(seq)
  starts <- seq.int(1L, seq_len, by = 60L)
  ends <- pmin(starts + 59L, seq_len)
  chunks <- substring(seq, starts, ends)

  vapply(seq_along(chunks), function(i) {
    chunk <- chunks[i]
    chunk_len <- nchar(chunk)
    grp_starts <- seq.int(1L, chunk_len, by = 10L)
    grp_ends <- pmin(grp_starts + 9L, chunk_len)
    groups <- substring(chunk, grp_starts, grp_ends)
    sprintf("%9d %s", starts[i], paste(groups, collapse = " "))
  }, character(1))
}

write_genbank_hits_gbk <- function(path, genome, features, cmap, dna_path) {
  seq_df <- read_fasta_sequences(dna_path)
  if (nrow(seq_df) == 0) {
    writeLines(character(0), path)
    return(invisible(NULL))
  }

  feature_df <- features |>
    left_join(cmap |> select(hmm, pathway), by = "hmm") |>
    arrange(contig, start, end, hmm)

  message(sprintf("[info] Writing BTEXgenie GenBank export -> %s", path))
  con <- file(path, open = "w")
  on.exit(close(con), add = TRUE)

  for (i in seq_len(nrow(seq_df))) {
    contig_name <- seq_df$contig[i]
    sequence <- seq_df$sequence[i]
    contig_features <- feature_df |>
      filter(contig == contig_name)

    locus_line <- sprintf(
      "LOCUS       %-16s %7d bp    DNA     linear   UNK 01-JAN-1980",
      substr(contig_name, 1, 16),
      nchar(sequence)
    )
    header_lines <- c(
      locus_line,
      sprintf("DEFINITION  BTEXgenie hits for %s contig %s.", genome, contig_name),
      sprintf("ACCESSION   %s", contig_name),
      sprintf("VERSION     %s", contig_name),
      "KEYWORDS    .",
      sprintf("SOURCE      %s", genome),
      sprintf("  ORGANISM  %s", genome),
      "            .",
      "FEATURES             Location/Qualifiers"
    )

    feature_lines <- c()
    if (nrow(contig_features) > 0) {
      for (j in seq_len(nrow(contig_features))) {
        feature_lines <- c(
          feature_lines,
          sprintf("     misc_feature    %s", format_genbank_location(contig_features$start[j], contig_features$end[j], contig_features$strand[j])),
          format_genbank_qualifier("gene", contig_features$hmm[j]),
          format_genbank_qualifier("note", contig_features$pathway[j]),
          format_genbank_qualifier("inference", "profile:BTEXgenie")
        )
      }
    }

    record_lines <- c(
      header_lines,
      feature_lines,
      "ORIGIN",
      format_origin_lines(sequence),
      "//"
    )
    writeLines(record_lines, con = con, sep = "\n")
  }

  message(sprintf("[info] Finished BTEXgenie GenBank export -> %s", path))
  invisible(path)
}

write_gc_skew_tsv <- function(path, gc_skew_df) {
  write_tsv(gc_skew_df, path)
}

write_kofam_hits_tsv <- function(path, kofam_hit_df) {
  write_tsv(kofam_hit_df, path)
}

write_kofam_density_summary_tsv <- function(path, density_summary_df) {
  write_tsv(density_summary_df, path)
}

write_kofam_density_matrix_tsv <- function(path, density_matrix_df) {
  write_tsv(density_matrix_df, path)
}

prepare_sector_data <- function(contig_lengths) {
  contig_lengths |>
    mutate(contig = factor(contig, levels = contig))
}

plot_ideogram_track <- function(sector_df, tick_major = 5e5, tick_minor = 1e5) {
  circos.trackPlotRegion(
    ylim = c(0, 1),
    bg.border = NA,
    track.height = 0.08,
    panel.fun = function(x, y) {
      sector.index <- get.cell.meta.data("sector.index")
      xlim <- get.cell.meta.data("xlim")

      circos.rect(xlim[1], 0, xlim[2], 1, col = "lightblue", border = "black", lwd = 1)

      max_x <- xlim[2]
      major_ticks <- sort(unique(c(seq(0, max_x, by = tick_major), max_x)))
      minor_ticks <- seq(0, max_x, by = tick_minor)

      if (length(minor_ticks) > 0) {
        circos.segments(minor_ticks, 1, minor_ticks, 1.08, col = "black", lwd = 0.5)
      }
      if (length(major_ticks) > 0) {
        circos.segments(major_ticks, 1, major_ticks, 1.15, col = "black", lwd = 1)
        circos.text(
          x = major_ticks,
          y = rep(1.22, length(major_ticks)),
          labels = paste0(round(major_ticks / 1e3, 1), " kb"),
          cex = 0.45,
          facing = "clockwise",
          niceFacing = TRUE,
          adj = c(0, 0.5)
        )
      }

      contig_name <- sector.index
      contig_center <- mean(xlim)

      circos.text(
        x = contig_center,
        y = CELL_META$ylim[2] + mm_y(2),
        labels = contig_name,
        sector.index = CELL_META$sector.index,
        track.index = CELL_META$track.index,
        facing = "clockwise",
        niceFacing = TRUE,
        adj = c(0, 0.5),
        cex = 0.5
      )
    }
  )
}

plot_gc_skew_track <- function(gc_skew_df, track_height = 0.10) {
  if (nrow(gc_skew_df) == 0) return(invisible(NULL))

  max_abs_skew <- max(abs(gc_skew_df$gc_skew), na.rm = TRUE)
  plot_limit <- max(0.05, max_abs_skew)

  circos.genomicTrackPlotRegion(
    gc_skew_df[, c("contig", "start", "end", "gc_skew")],
    ylim = c(-plot_limit, plot_limit),
    bg.border = "black",
    track.height = track_height,
    panel.fun = function(region, value, ...) {
      if (nrow(region) == 0) return()
      skew_vals <- value[[1]]
      colors <- ifelse(skew_vals >= 0, "#E67E22", "#2E86C1")
      circos.genomicRect(
        region,
        value,
        ybottom = pmin(skew_vals, 0),
        ytop = pmax(skew_vals, 0),
        col = colors,
        border = NA,
        ...
      )
      xlim <- get.cell.meta.data("xlim")
      circos.segments(xlim[1], 0, xlim[2], 0, lty = 2, col = "grey35")
    }
  )
}

plot_kofam_density_track <- function(kofam_density_summary_df, track_height = 0.07) {
  if (is.null(kofam_density_summary_df) || nrow(kofam_density_summary_df) == 0) return(invisible(NULL))

  plot_df <- kofam_density_summary_df |>
    filter(window_value > 0) |>
    transmute(contig, start, end, fill)

  if (nrow(plot_df) == 0) {
    circos.trackPlotRegion(ylim = c(0, 1), bg.border = "black", track.height = track_height, panel.fun = function(x, y) {})
    return(invisible(NULL))
  }

  circos.genomicTrackPlotRegion(
    plot_df[, c("contig", "start", "end", "fill")],
    ylim = c(0, 1),
    bg.border = "black",
    track.height = track_height,
    panel.fun = function(region, value, ...) {
      if (nrow(region) == 0) return()
      circos.genomicRect(region, value, ybottom = 0, ytop = 1, col = value[[1]], border = NA, ...)
    }
  )
}

plot_feature_track <- function(track_features, cmap, track_height = 0.08, label_cex = 0.45, show_links = FALSE) {
  if (nrow(track_features) == 0) {
    circos.trackPlotRegion(ylim = c(0, 1), bg.border = NA, track.height = track_height, panel.fun = function(x, y) {})
    return(invisible(NULL))
  }

  feat <- track_features |>
    left_join(cmap |> select(hmm, color_hex), by = "hmm")

  circos.genomicTrackPlotRegion(
    feat[, c("contig", "start", "end", "hmm", "color_hex")],
    ylim = c(0, 1),
    bg.border = NA,
    track.height = track_height,
    panel.fun = function(region, value, ...) {
      if (nrow(region) == 0) return()
      circos.genomicRect(region, value, ytop = 1, ybottom = 0, col = value[["color_hex"]], border = NA, ...)
    }
  )

  lab_df <- feat |>
    distinct(contig, start, end, hmm) %>%
    mutate(chr = contig) %>%
    select(chr, start, end, hmm)

  if (nrow(lab_df) > 0) {
    circos.genomicLabels(
      lab_df,
      labels.column = 4,
      labels.side = "inside",
      cex = label_cex,
      col = "black",
      line_col = "black",
      connection_height = mm_h(5)
    )
  }
}

make_circlize_plot <- function(sample_dir, genome, contig_lengths, features, cmap, gc_skew_df = NULL, kofam_density_summary_df = NULL, pdf_width = 10, pdf_height = 10, show_kofam_density = TRUE) {
  pdf_path <- file.path(sample_dir, "circos_plot.pdf")
  png_path <- file.path(sample_dir, "circos_plot.png")
  gc_skew_pos_color <- "#E67E22"
  gc_skew_neg_color <- "#2E86C1"

  if (!is.numeric(pdf_width) || length(pdf_width) != 1 || is.na(pdf_width) || !is.finite(pdf_width) || pdf_width <= 0) {
    stop(sprintf("[ERROR] --pdf-width must be a positive finite number; got: %s", as.character(pdf_width)), call. = FALSE)
  }
  if (!is.numeric(pdf_height) || length(pdf_height) != 1 || is.na(pdf_height) || !is.finite(pdf_height) || pdf_height <= 0) {
    stop(sprintf("[ERROR] --pdf-height must be a positive finite number; got: %s", as.character(pdf_height)), call. = FALSE)
  }

  sector_df <- prepare_sector_data(contig_lengths)
  xlim_mat <- cbind(rep(0, nrow(sector_df)), sector_df$length)
  rownames(xlim_mat) <- sector_df$contig

  draw_once <- function(device_fun) {
    device_fun()
    on.exit(dev.off(), add = TRUE)

    par(mar = c(5, 5, 5, 5))

    circos.clear()
    circos.par(
      start.degree = 90,
      gap.degree = rep(1, nrow(sector_df)),
      cell.padding = c(0, 0, 0, 0),
      track.margin = c(0.01, 0.01),
      points.overflow.warning = FALSE,
      canvas.xlim = c(-1.35, 1.35),
      canvas.ylim = c(-1.35, 1.35)
    )

    circos.initialize(factors = as.character(sector_df$contig), xlim = xlim_mat)
    plot_ideogram_track(sector_df)

    if (!is.null(gc_skew_df) && nrow(gc_skew_df) > 0) {
      plot_gc_skew_track(gc_skew_df, track_height = 0.10)
    }
    if (isTRUE(show_kofam_density) && !is.null(kofam_density_summary_df) && nrow(kofam_density_summary_df) > 0) {
      plot_kofam_density_track(kofam_density_summary_df, track_height = 0.07)
    }
    plot_feature_track(features, cmap, track_height = 0.09, label_cex = 0.45, show_links = FALSE)

    text(
      x = 0,
      y = 0,
      labels = genome,
      cex = 1.2,
      font = 2
    )
    rect(xleft = -0.17, ybottom = -0.1125, xright = -0.1625, ytop = -0.0675, col = gc_skew_pos_color, border = NA)
    text(x = -0.145, y = -0.09, labels = "Positive GC skew", adj = c(0, 0.5), cex = 0.55)
    rect(xleft = -0.17, ybottom = -0.2025, xright = -0.1625, ytop = -0.1575, col = gc_skew_neg_color, border = NA)
    text(x = -0.145, y = -0.18, labels = "Negative GC skew", adj = c(0, 0.5), cex = 0.55)
    circos.clear()
  }

  draw_once(function() grDevices::cairo_pdf(pdf_path, width = pdf_width, height = pdf_height))
  draw_once(function() {
    png_width_px <- max(1L, as.integer(round(pdf_width * 300)))
    png_height_px <- max(1L, as.integer(round(pdf_height * 300)))
    png(png_path, width = png_width_px, height = png_height_px, res = 300, type = "cairo")
  })

  invisible(list(pdf = pdf_path, png = png_path))
}

DEFAULT_SHOW_KOFAM_DENSITY <- TRUE
DEFAULT_PDF_WIDTH <- 20
DEFAULT_PDF_HEIGHT <- 20

run <- function(opt) {
  genome <- opt$sample
  if (is.null(genome) || !nzchar(genome)) stop("[ERROR] --sample is required", call. = FALSE)
  if (is.null(opt$hmmscan) || !file.exists(opt$hmmscan)) stop(sprintf("[ERROR] hmmscan file not found: %s", opt$hmmscan), call. = FALSE)
  if (is.null(opt$outdir) || !nzchar(opt$outdir)) stop("[ERROR] --outdir is required", call. = FALSE)

  pathway_map_path <- if (!is.null(opt$pathway_map) && nzchar(opt$pathway_map)) opt$pathway_map else infer_default_pathway_map()
  if (!file.exists(pathway_map_path)) stop(sprintf("[ERROR] Pathway map file not found: %s", pathway_map_path), call. = FALSE)
  message(sprintf("[log] Using pathway map: %s", pathway_map_path))

  sample_dir <- file.path(opt$outdir, paste0(genome, "_circos_plot"))
  ensure_dir(sample_dir)

  pathway_df <- load_pathway_map(pathway_map_path)

  contig_lengths_path <- generate_contig_lengths_for_sample(opt, sample_dir)
  contig_lengths <- read_contig_lengths(contig_lengths_path, genome)
  if (nrow(contig_lengths) == 0) stop(sprintf("[ERROR] No contig lengths found for genome: %s", genome), call. = FALSE)

  features <- read_hmmscan_features(opt$hmmscan, genome)

  hmms_with_hits <- sort(unique(features$hmm))

  cmap <- if (length(hmms_with_hits) > 0) build_color_map(hmms_with_hits, pathway_df) else tibble(
    hmm = character(), color_name = character(), rgb = character(), color_hex = character(), pathway = character()
  )

  gc_skew_df <- compute_gc_skew_windows(
    dna_path = opt$genome,
    contigs_to_keep = contig_lengths$contig,
    window_size = opt$window_size,
    step_size = opt$window_size
  )

  kofam_hit_df <- tibble()
  kofam_density_summary_df <- tibble()
  kofam_density_matrix_df <- tibble()
  if (!is.null(opt$kofam) && nzchar(opt$kofam)) {
    if (is.null(opt$prodigal_gbk) || !nzchar(opt$prodigal_gbk)) {
      stop("[ERROR] --prodigal-gbk is required when using --kofam and --ko-category-map.", call. = FALSE)
    }
    ko_category_map_path <- if (!is.null(opt$ko_category_map) && nzchar(opt$ko_category_map)) opt$ko_category_map else infer_default_ko_category_map()
    message(sprintf("[info] Using KOfam assignments: %s", opt$kofam))
    message(sprintf("[info] Using KO category map: %s", ko_category_map_path))
    cds_df <- read_prodigal_gbk_cds(opt$prodigal_gbk)
    kofam_df <- read_kofam_assignments(opt$kofam)
    ko_cat_df <- read_ko_category_map(ko_category_map_path)
    kofam_hit_df <- join_kofam_to_cds(kofam_df, cds_df, ko_cat_df)
    density_out <- compute_category_density_windows(kofam_hit_df, contig_lengths, window_size = opt$window_size)
    kofam_density_summary_df <- density_out$summary
    kofam_density_matrix_df <- density_out$matrix
  }

  filter_to_plot_contigs <- NULL

  if (isTRUE(opt$plot_btex_hit_contigs_only)) {
    message("[info] Restricting plot to contigs with BTEX-HMM hits only")
    filter_to_plot_contigs <- unique(features$contig)
    filter_to_plot_contigs <- filter_to_plot_contigs[!is.na(filter_to_plot_contigs) & nzchar(filter_to_plot_contigs)]

    if (length(filter_to_plot_contigs) == 0) {
      stop("[ERROR] --plot-btex-hit-contigs-only was set, but no contigs with BTEX-HMM hits were found to plot.", call. = FALSE)
    }
  } else if (nrow(contig_lengths) > 10) {
    message("[WARNING] More than 10 contigs detected; to enhance readability, only the contigs with BTEXgenie hits are shown")
    filter_to_plot_contigs <- unique(features$contig)
    filter_to_plot_contigs <- filter_to_plot_contigs[!is.na(filter_to_plot_contigs) & nzchar(filter_to_plot_contigs)]

    if (length(filter_to_plot_contigs) == 0) {
      stop("[ERROR] More than 10 contigs detected, but no contigs with BTEX-HMM hits were found to plot.", call. = FALSE)
    }
  }

  if (!is.null(filter_to_plot_contigs)) {
    contig_lengths <- contig_lengths |>
      filter(contig %in% filter_to_plot_contigs)

    if (nrow(contig_lengths) == 0) {
      stop("[ERROR] No contigs remained after filtering the contig set for plotting.", call. = FALSE)
    }

    features <- features |>
      filter(contig %in% contig_lengths$contig)

    gc_skew_df <- gc_skew_df |>
      filter(contig %in% contig_lengths$contig)

    kofam_hit_df <- kofam_hit_df |>
      filter(contig %in% contig_lengths$contig)

    kofam_density_summary_df <- kofam_density_summary_df |>
      filter(contig %in% contig_lengths$contig)

    kofam_density_matrix_df <- kofam_density_matrix_df |>
      filter(contig %in% contig_lengths$contig)
  }

  write_hmm_colors_tsv(file.path(sample_dir, "hmm_colors.tsv"), cmap)
  write_karyotype_tsv(file.path(sample_dir, "karyotype.tsv"), contig_lengths)
  write_gene_labels_tsv(file.path(sample_dir, "gene_labels.tsv"), features)
  write_genbank_hits_gbk(file.path(sample_dir, "btex_genie_hits.gbk"), genome, features, cmap, opt$genome)
  write_gc_skew_tsv(file.path(sample_dir, "gc_skew_windows.tsv"), gc_skew_df)

  if (nrow(kofam_hit_df) > 0) {
    write_kofam_hits_tsv(file.path(sample_dir, "kofam_category_hits.tsv"), kofam_hit_df)
  }
  plot_paths <- make_circlize_plot(
    sample_dir = sample_dir,
    genome = genome,
    contig_lengths = contig_lengths,
    features = features,
    cmap = cmap,
    gc_skew_df = gc_skew_df,
    kofam_density_summary_df = kofam_density_summary_df,
    pdf_width = DEFAULT_PDF_WIDTH,
    pdf_height = DEFAULT_PDF_HEIGHT,
    show_kofam_density = DEFAULT_SHOW_KOFAM_DENSITY
  )

  if (nrow(kofam_density_summary_df) > 0) {
    message("[OK] Added KOfam category density track using dominant category per window and density-shaded fills")
  }
  message(sprintf("[OK] circlize project for %s -> %s", genome, sample_dir))
  message(sprintf("[OK] Plot PDF -> %s", plot_paths$pdf))
  message(sprintf("[OK] Plot PNG -> %s", plot_paths$png))

  invisible(sample_dir)
}

option_list <- list(
  make_option("--hmmscan", type = "character", help = "Path to btex_hmm_summary.csv or similar hmmscan summary CSV"),
  make_option(c("-o", "--outdir"), type = "character", help = "Output directory"),
  make_option(c("-g", "--genome"), type = "character", default = NULL, help = "Genome FASTA file for the selected sample"),
  make_option("--prodigal-gbk", dest = "prodigal_gbk", type = "character", default = NULL, help = "Prodigal GenBank file used to map KOfam hits to genomic coordinates"),
  make_option(c("-s", "--sample"), type = "character", help = "Sample/genome name to plot"),
  make_option("--pathway-map", dest = "pathway_map", type = "character", default = NULL, help = "Optional override for pathway_map.tsv"),
  make_option("--kofam", type = "character", default = NULL, help = "Optional KOfamScan detail output or a table with query_id and ko columns, including gene_name and KO style tables"),
  make_option("--ko-category-map", dest = "ko_category_map", type = "character", default = NULL, help = "Optional table with columns KO, category, color for KOfam density plotting"),
  make_option("--plot-btex-hit-contigs-only", dest = "plot_btex_hit_contigs_only", action = "store_true", default = FALSE, help = "Only plot contigs with BTEX-HMM hits, excluding contigs retained solely by KOfam hits"),
  make_option("--window-size", dest = "window_size", type = "integer", default = 5000L, help = "Window size for GC skew and KOfam category density [default %default]")
)

parser <- OptionParser(option_list = option_list)
opt <- parse_args(parser)

if (is.null(opt$hmmscan) || is.null(opt$outdir) || is.null(opt$sample) || is.null(opt$genome)) {
  print_help(parser)
  quit(status = 1)
}

run(opt)
