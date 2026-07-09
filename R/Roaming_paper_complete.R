# ============================================================================
# FREE-ROAMING DOG DENSITY FROM SURVEY TRANSECTS IN CHILE
# Complete reproducible pipeline — Sections A–E
# ============================================================================
#
# INPUT FILES (place in data_dir):
#   roaming.xlsx        sheet "original" — Epicollect5 transect records
#                       sheet "survey"   — district attributes (25 rows × 80 cols)
#   bites_2019_24.xlsx  sheet "base_analisis" — SIRAM bite registry 2019–2024
#
# SECTIONS:
#   0  Setup & district audit
#   A  Density estimation       dogs/km per district    → Table 1, Fig 1
#   B  Temporal structure       NB GLMM day/hour        → Table 2, Figs 2–3
#   C  Original correlates      density ~ survey vars   → Table 3, Table S2
#   D  Public-health link       bites ~ density         → Table 4
#   E  Extended correlates      IDC2, IDH, behaviour    → Tables E1–E4, Figs E1–E6
#
# ── VARIABLE AUDIT (verified against current roaming.xlsx) ──────────────────
# DELETED from sheet (no longer present):
#   gatos_urbano_UC, gatos_rurales_UC, total_gatos_SUMA_UC, densidad_gatos_km2
#   perro:viv_urb, perro:viv_rur  (computed in R from raw numerators)
#   z_extrema, n_region, cod.com, provincia  (never used)
#
# PRESENT but EXCLUDED from analysis (with reasons):
#   total_mascotas, densidad_mascotas_km2  — dogs+cats combined (~1.44× dogs)
#   neutering    — old binary (only 2 districts ≠ 0); 'neutered' used instead
#   IDC          — identical to IDC2 (r = 1.00); IDC2 retained
#   n_vets       — identical to vets (r = 1.00); vets retained
#   n_hogares    — near-identical to hogares (r = 0.9998); hogares retained
#   densidad     — rounded integer from IDC source (r = 0.69 with
#                  densidad_humana_km2, ratio 0.8–140×); NOT the same variable
#                  → use densidad_humana_km2 (precise pop/area_km2)
#   Superficie   — major discrepancies from 'area km2' for several districts
#                  (max diff 8796 km2); different source → use 'area km2'
#   porcentaje_urb / ubr.pop.percent — both measure % urban pop but not
#                  identical; porcentaje_urb used (from same survey source)
#   adj2018/19/20 — 4–8/25 valid; too sparse
#   bites2014/15  — all zero or near-zero
#   tot.albergues, alb.*  — 15–16/25 valid; too sparse
#   RANKING, RANGOS, IDH2005_cat — categorical/ranking only
#
# CLEAN_NAMES() HAZARD:
#   ECONOMÍA → econom_a, EDUCACIÓN → educaci_n, 'always inside' → always_inside
#   These columns are RENAMED EXPLICITLY before clean_names() is called.
# ============================================================================

rm(list = ls()); graphics.off(); set.seed(1)

# ── 0a. Paths & directories ──────────────────────────────────────────────────
data_dir  <- "."            # <- set to folder containing roaming.xlsx
roam_file <- "roaming.xlsx"
bite_file <- "bites_2019_24.xlsx"

setwd(data_dir)
for (d in c("data", "tables", "figures")) dir.create(d, showWarnings = FALSE)

# ── 0b. Packages ─────────────────────────────────────────────────────────────
library(MASS)   # load before tidyverse so dplyr::select wins
library(readxl); library(writexl); library(tidyverse); library(janitor)
library(lubridate); library(geosphere); library(splines); library(gridExtra)
library(glmmTMB); library(car); library(DHARMa); library(performance)
library(broom.mixed); library(emmeans); library(broom)
select <- dplyr::select
filter <- dplyr::filter

# Section-E packages — auto-install if absent
e_pkgs  <- c("ggcorrplot", "FactoMineR", "factoextra", "ggrepel")
to_inst <- e_pkgs[!sapply(e_pkgs, requireNamespace, quietly = TRUE)]
if (length(to_inst)) install.packages(to_inst, repos = "https://cloud.r-project.org")
library(ggcorrplot); library(FactoMineR); library(factoextra); library(ggrepel)

# ── 0c. Parameters ───────────────────────────────────────────────────────────
jump_threshold_m   <- 500     # GPS step > this = relocation/noise
min_walked_km      <- 1.0     # district low-effort threshold
min_cell_km        <- 0.1     # GLMM temporal cell minimum
day_coding         <- "survey_day"
bite_years         <- 2019:2024
summer_months      <- c(12, 1, 2)   # austral summer

university_comunas <- c("VALPARAISO","LAS CONDES","SAN VICENTE","LONQUIMAY",
                        "PADRE HURTADO","EL TABO","MULCHEN","ILLAPEL","LONCOCHE")
team_colors <- c(University = "#d7301f", Municipality = "#2c7fb8")
theme_set(theme_bw(base_size = 11))

pois_lwr <- function(x, t) ifelse(x == 0, 0, qgamma(0.025, x) / t)
pois_upr <- function(x, t) qgamma(0.975, x + 1) / t

# ============================================================================
# 0d. READ SURVEY SHEET
#     Rename accented / special-character columns BEFORE clean_names()
# ============================================================================
survey_wide <- read_excel(roam_file, sheet = "survey")
names(survey_wide) <- trimws(names(survey_wide))

# Explicit renames to protect from clean_names() mangling
survey_wide <- survey_wide %>%
  rename(
    area_km2          = `area km2`,
    persona_perro     = `persona:perro`,
    idh_2005          = `IDH 2005`,
    idc2_raw          = IDC2,
    bienestar_raw     = BIENESTAR,
    economia_raw      = `ECONOMÍA`,
    educacion_raw     = `EDUCACIÓN`,
    always_inside_raw = `always inside`,
    grand_total_raw   = `Grand Total`,
    sale_solo_raw     = `Sale solo`,
    siempre_casa_raw  = `Siempre en la casa`,
    sale_acomp_raw    = `Sale acompañado (con supervisión)`
  )

survey_raw <- survey_wide %>%
  clean_names() %>%
  mutate(
    comuna        = str_squish(toupper(as.character(comuna))),
    # Restore rescued columns with clean names
    idc2          = as.numeric(idc2_raw),
    bienestar     = as.numeric(bienestar_raw),
    economia      = as.numeric(economia_raw),
    educacion     = as.numeric(educacion_raw),
    always_inside = as.numeric(always_inside_raw),
    grand_total   = as.numeric(grand_total_raw),
    sale_solo     = as.numeric(sale_solo_raw),
    siempre_casa  = as.numeric(siempre_casa_raw),
    sale_acomp    = as.numeric(sale_acomp_raw)
  )

# Coerce remaining columns to numeric (skip known text fields)
sv <- survey_raw %>%
  mutate(across(
    -c(comuna, macrozona, alb_tipo, region_name, rangos,
       idh_2005_cat, area_urm, encuesta_terreno),
    ~ suppressWarnings(as.numeric(.))
  ))

chosen <- sort(unique(sv$comuna))   # 25 districts

# ============================================================================
# 0e. TRANSECT RECORDS + DISTRICT AUDIT
# ============================================================================
raw <- read_excel(roam_file, sheet = "original") %>%
  clean_names() %>%
  mutate(
    comuna           = str_squish(toupper(as.character(comuna))),
    created_at_local = with_tz(ymd_hms(created_at, tz = "UTC"), "America/Santiago"),
    survey_date      = as_date(created_at_local),
    weekday          = wday(created_at_local, label = TRUE, abbr = TRUE, week_start = 1),
    hour             = hour(suppressWarnings(hms(x5_hora))),
    n_perros         = as.numeric(n_perros),
    lat              = as.numeric(lat_2_geolocalizacin),
    lon              = as.numeric(long_2_geolocalizacin),
    in_window        = !is.na(hour) & hour >= 8 & hour <= 19 & !(hour >= 12 & hour < 14),
    valid_geo        = !is.na(lat) & !is.na(lon),
    valid_dogs       = !is.na(n_perros),
    kept             = in_window & valid_geo & valid_dogs
  )

audit <- raw %>%
  group_by(comuna) %>%
  summarise(
    n_raw        = n(),
    n_kept       = sum(kept),
    n_out_window = sum(!in_window),
    n_days       = n_distinct(survey_date[kept]),
    .groups = "drop"
  ) %>%
  mutate(disposition = case_when(
    comuna %in% chosen ~ "analysed",
    n_kept == 0        ~ "no records in survey window",
    n_kept <= 1        ~ "single point — not estimable",
    TRUE               ~ "excluded"
  ))
write.csv(audit, "tables/TableS1_district_audit.csv", row.names = FALSE)

# ============================================================================
# A. DENSITY ESTIMATION  (Table 1, Figure 1)
# ============================================================================
roam <- raw %>%
  filter(kept, comuna %in% chosen) %>%
  arrange(comuna, survey_date, created_at_local) %>%
  group_by(comuna, survey_date) %>%
  mutate(
    seg_m   = distGeo(cbind(lag(lon), lag(lat)), cbind(lon, lat)),
    is_jump = !is.na(seg_m) & seg_m > jump_threshold_m
  ) %>%
  ungroup() %>%
  mutate(
    dogs_present = n_perros > 0,
    team = factor(
      if_else(comuna %in% university_comunas, "University", "Municipality"),
      levels = c("Municipality", "University")
    )
  )

transect <- roam %>%
  group_by(comuna, team) %>%
  summarise(
    n_records     = n(),
    n_survey_days = n_distinct(survey_date),
    n_dog_records = sum(dogs_present),
    total_dogs    = sum(n_perros),
    walked_km     = sum(seg_m[!is_jump], na.rm = TRUE) / 1000,
    raw_km        = sum(seg_m, na.rm = TRUE) / 1000,
    n_jumps       = sum(is_jump, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    dogs_per_km     = total_dogs / walked_km,
    dogs_per_km_lcl = pois_lwr(total_dogs, walked_km),
    dogs_per_km_ucl = pois_upr(total_dogs, walked_km),
    records_per_km  = n_records / walked_km,
    encounter_rate  = n_dog_records / n_records,
    effort_flag     = if_else(walked_km >= min_walked_km, "ok", "low-effort")
  ) %>%
  arrange(desc(dogs_per_km))

write_xlsx(transect, "data/transect_estimates.xlsx")
write.csv(transect, "tables/Table1_density_by_district.csv", row.names = FALSE)

overall <- transect %>%
  summarise(total_dogs = sum(total_dogs), walked_km = sum(walked_km)) %>%
  mutate(
    dogs_per_km = total_dogs / walked_km,
    lcl         = pois_lwr(total_dogs, walked_km),
    ucl         = pois_upr(total_dogs, walked_km)
  )

fig1 <- transect %>%
  mutate(reliable = walked_km >= min_walked_km) %>%
  arrange(dogs_per_km) %>%
  mutate(comuna = factor(comuna, levels = comuna)) %>%
  ggplot(aes(dogs_per_km, comuna)) +
  geom_vline(xintercept = overall$dogs_per_km, linetype = "dashed", colour = "grey50") +
  geom_errorbarh(aes(xmin = dogs_per_km_lcl, xmax = dogs_per_km_ucl),
                 height = 0, colour = "grey70") +
  geom_point(aes(size = walked_km, colour = team, shape = reliable)) +
  scale_colour_manual(values = team_colors, name = "Counting team") +
  scale_shape_manual(values = c(`TRUE` = 16, `FALSE` = 1),
                     labels  = c("\u2265 1 km", "< 1 km"), name = "Effort") +
  scale_size_continuous(range = c(1.5, 6), name = "Km walked") +
  theme(panel.grid.major.y = element_blank()) +
  labs(title = "Roaming-dog density by district", x = "Dogs per km walked", y = NULL,
       subtitle = sprintf("95%% Poisson CI; dashed = pooled (%.1f dogs/km)", overall$dogs_per_km))
ggsave("figures/Fig1_density_by_district.png", fig1, width = 8, height = 6, dpi = 300)

# ============================================================================
# B. TEMPORAL STRUCTURE  (Figure 2, Table 2, Figure 3)
# ============================================================================
density_by <- function(data, ...) {
  g  <- enquos(...)
  d1 <- data %>% group_by(!!!g) %>% summarise(dogs = sum(n_perros), .groups = "drop")
  d2 <- data %>% filter(!is_jump) %>% group_by(!!!g) %>%
    summarise(walked_km = sum(seg_m, na.rm = TRUE) / 1000, .groups = "drop")
  left_join(d1, d2, by = sapply(g, rlang::as_label)) %>%
    mutate(
      walked_km   = replace_na(walked_km, 0),
      dogs_per_km = if_else(walked_km >= 0.5, dogs / walked_km, NA_real_),
      lcl         = if_else(walked_km >= 0.5, pois_lwr(dogs, walked_km), NA_real_),
      ucl         = if_else(walked_km >= 0.5, pois_upr(dogs, walked_km), NA_real_)
    )
}

by_hour <- density_by(roam, hour)    %>% filter(!is.na(dogs_per_km)) %>% arrange(hour)
by_day  <- density_by(roam, weekday) %>% filter(!is.na(dogs_per_km))

f2a <- ggplot(by_hour, aes(hour, dogs_per_km)) +
  geom_col(fill = "#2c7fb8") +
  geom_errorbar(aes(ymin = lcl, ymax = ucl), width = .25, colour = "grey40") +
  scale_x_continuous(breaks = 8:19) +
  labs(title = "By hour", x = "Hour (local)", y = "Dogs per km")
f2b <- ggplot(by_day, aes(weekday, dogs_per_km)) +
  geom_col(fill = "#41ab5d") +
  geom_errorbar(aes(ymin = lcl, ymax = ucl), width = .25, colour = "grey40") +
  labs(title = "By day of week", x = NULL, y = "Dogs per km")
ggsave("figures/Fig2_temporal.png", grid.arrange(f2a, f2b, ncol = 2),
       width = 10, height = 4, dpi = 300)

# Cell aggregation (district × survey-day × hour)
cells <- roam %>%
  group_by(comuna, survey_date, weekday, hour, team) %>%
  summarise(
    dogs      = sum(n_perros),
    walked_km = sum(seg_m[!is_jump], na.rm = TRUE) / 1000,
    .groups = "drop"
  ) %>%
  group_by(comuna) %>%
  mutate(survey_day = dense_rank(survey_date)) %>%
  ungroup() %>%
  dplyr::filter(walked_km >= min_cell_km) %>%
  mutate(
    district   = factor(comuna),
    hour_f     = factor(hour),
    weekday    = factor(weekday, ordered = FALSE),
    survey_day = factor(survey_day),
    log_km     = log(walked_km)
  )
cells$day <- if (day_coding == "weekday") cells$weekday else cells$survey_day

m_pois <- glmmTMB(dogs ~ day + hour_f + (1|district) + offset(log_km), data = cells, family = poisson)
m_nb   <- glmmTMB(dogs ~ day + hour_f + (1|district) + offset(log_km), data = cells, family = nbinom2)
m_null <- glmmTMB(dogs ~ day + (1|district) + offset(log_km),           data = cells, family = nbinom2)
m_spl  <- glmmTMB(dogs ~ day + ns(hour,3) + (1|district) + offset(log_km), data = cells, family = nbinom2)
m_team <- glmmTMB(dogs ~ team + day + ns(hour,3) + (1|district) + offset(log_km), data = cells, family = nbinom2)
m_int  <- glmmTMB(dogs ~ team*ns(hour,3) + day + (1|district) + offset(log_km), data = cells, family = nbinom2)
m_nore <- glmmTMB(dogs ~ day + hour_f + offset(log_km), data = cells, family = nbinom2)

glmm_results <- list(
  AIC_poisson_vs_nb = AIC(m_pois, m_nb),
  typeII            = car::Anova(m_nb, type = "II"),
  LRT_day           = anova(update(m_nb, . ~ . - day), m_nb),
  LRT_hour          = anova(update(m_nb, . ~ . - hour_f), m_nb),
  LRT_district_RE   = anova(m_nore, m_nb),
  ICC               = performance::icc(m_nb),
  LRT_hour_smooth   = anova(m_null, m_spl),
  AIC_all           = AIC(m_null, m_nb, m_spl, m_team, m_int),
  LRT_team          = anova(m_spl, m_team),
  team_RR           = tidy(m_team, effects = "fixed", exponentiate = TRUE, conf.int = TRUE),
  LRT_team_x_hour   = anova(m_team, m_int)
)
capture.output(glmm_results, file = "tables/Table2_glmm.txt")
png("figures/FigS_dharma.png", 1000, 500); plot(simulateResiduals(m_nb)); dev.off()

fig3 <- emmeans(m_spl, ~ hour,
                at = list(hour = seq(9, 19, 0.5)), offset = 0, type = "response") %>%
  as.data.frame() %>%
  ggplot(aes(hour, response)) +
  annotate("rect", xmin = 12, xmax = 14, ymin = -Inf, ymax = Inf, alpha = .08, fill = "grey50") +
  geom_ribbon(aes(ymin = asymp.LCL, ymax = asymp.UCL), alpha = .15, fill = "#2c7fb8") +
  geom_line(linewidth = 1, colour = "#2c7fb8") +
  scale_x_continuous(breaks = 8:19) +
  labs(title = "Adjusted diurnal density (NB GLMM)", x = "Hour (local)",
       y = "Adjusted dogs per km",
       subtitle = "ns(hour,3) + district random intercept; band=95% CI; grey=lunch gap")
ggsave("figures/Fig3_glmm_diurnal.png", fig3, width = 7, height = 4.5, dpi = 300)

# ============================================================================
# C / D. REGRESSION DATASET
#   Dogs only. Ratios dogs/dwelling computed from raw numerators (both present).
#   densidad_humana_km2 used for human density (NOT 'densidad').
#   area_km2 used for area (NOT 'Superficie').
#   porcentaje_urb used for % urban (NOT ubr.pop.percent).
# ============================================================================
surv <- sv %>%
  transmute(
    comuna, region, macrozona,
    pop_total               = poblacion_total,
    pop_urban               = urbano_poblacion,       # 24/25 valid
    pop_rural               = rural_poblacion,        # 22/25 valid
    dwellings_urban         = viviendas_urbanas,      # 24/25 valid
    dwellings_rural         = viviendas_rurales,      # 22/25 valid
    area_km2,                                         # 25/25 valid
    human_density_km2       = densidad_humana_km2,    # 25/25 valid; precise pop/area
    owned_dogs_total        = total_perros_suma_uc,   # 25/25 valid
    owned_dog_density_km2   = densidad_perros_km2,    # 25/25 valid
    human_dog_ratio         = persona_perro,          # 25/25 valid
    pct_urban_pop           = porcentaje_urb,         # 25/25 valid; 0–100 scale
    # Dogs-per-dwelling computed from raw columns (perro:viv_urb/rur deleted)
    dogs_per_urban_dwelling = perros_urbano_uc / viviendas_urbanas,   # 24/25
    dogs_per_rural_dwelling = perros_rurales_uc / viviendas_rurales   # 22/25
  )

# Bite outcomes
bites <- read_excel(bite_file, sheet = "base_analisis") %>%
  clean_names() %>%
  mutate(
    comuna    = str_squish(toupper(as.character(comuna_ocurrio_mordedura2))),
    year      = year(as_date(fecha_mordedura)),
    month     = month(as_date(fecha_mordedura)),
    public    = lugar_mordedura2 == "Lugar P\u00fablico",
    ownerless = ubicacion_animal_mordedor == "SIN DIRECCI\u00d3N",
    summer    = month %in% summer_months
  )

pop22 <- bites %>%
  group_by(comuna) %>% summarise(pop = median(poblacion_2022), .groups = "drop")

bite_sum <- bites %>%
  filter(year %in% bite_years) %>%
  group_by(comuna) %>%
  summarise(
    bites_public_total            = sum(public, na.rm = TRUE),
    bites_ownerless_public_total  = sum(public & ownerless, na.rm = TRUE),
    bites_public_summer           = sum(public & summer, na.rm = TRUE),
    bites_ownerless_public_summer = sum(public & ownerless & summer, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  left_join(pop22, by = "comuna") %>%
  mutate(
    bites_public_inc                  = bites_public_total / length(bite_years) / pop * 1e5,
    bites_ownerless_public_inc        = bites_ownerless_public_total / length(bite_years) / pop * 1e5,
    bites_public_summer_inc           = bites_public_summer / length(bite_years) / pop * 1e5,
    bites_ownerless_public_summer_inc = bites_ownerless_public_summer / length(bite_years) / pop * 1e5
  ) %>%
  rename(bite_pop_2022 = pop)

reg <- transect %>%
  select(comuna, team, total_dogs, walked_km, dogs_per_km,
         dogs_per_km_lcl, dogs_per_km_ucl, encounter_rate, effort_flag) %>%
  left_join(surv,     by = "comuna") %>%
  left_join(bite_sum, by = "comuna") %>%
  mutate(log_km = log(walked_km), low_effort = walked_km < min_walked_km)

reg_he <- reg %>% filter(!low_effort)

write_xlsx(reg, "data/regression_data.xlsx")
write.csv(reg,  "tables/regression_dataset.csv", row.names = FALSE)

# ============================================================================
# C. ORIGINAL CORRELATES  (Table 3, Table S2)
# ============================================================================
# Dogs-only predictors; cat/pet columns not included.
logv <- c("owned_dog_density_km2", "human_density_km2",
          "owned_dogs_total", "pop_total", "area_km2",
          "dogs_per_urban_dwelling", "dogs_per_rural_dwelling")
rawv <- c("human_dog_ratio", "pct_urban_pop")

uni_nb <- function(v, dat, logit) {
  dd <- dat %>%
    mutate(.x = if (logit) log(.data[[v]]) else .data[[v]]) %>%
    filter(is.finite(.x)) %>% mutate(z = as.numeric(scale(.x)))
  m <- tryCatch(glm.nb(total_dogs ~ z + offset(log_km), data = dd),
                error = function(e) NULL)
  if (is.null(m))
    return(tibble(predictor = v, RR = NA, lwr = NA, upr = NA, p = NA, n = nrow(dd)))
  s <- tidy(m, exponentiate = TRUE, conf.int = TRUE) %>% filter(term == "z")
  tibble(predictor = v, RR = s$estimate, lwr = s$conf.low, upr = s$conf.high,
         p = s$p.value, n = nrow(dd))
}

Table3 <- bind_rows(
  map_dfr(logv, uni_nb, dat = reg,    logit = TRUE),
  map_dfr(rawv, uni_nb, dat = reg,    logit = FALSE)
) %>% arrange(p)
write.csv(Table3, "tables/Table3_density_correlates.csv", row.names = FALSE)

m0 <- glm.nb(total_dogs ~ 1 + offset(log_km), data = reg)
capture.output(
  list(
    team      = anova(m0, glm.nb(total_dogs ~ team + offset(log_km), data = reg)),
    macrozona = anova(m0, glm.nb(total_dogs ~ factor(macrozona) + offset(log_km), data = reg))
  ),
  file = "tables/Table3b_factor_tests.txt"
)

TableS2 <- bind_rows(
  map_dfr(logv, uni_nb, dat = reg_he, logit = TRUE),
  map_dfr(rawv, uni_nb, dat = reg_he, logit = FALSE)
) %>% arrange(p)
write.csv(TableS2, "tables/TableS2_correlates_sensitivity.csv", row.names = FALSE)

# ============================================================================
# D. PUBLIC-HEALTH LINK  (Table 4)
# ============================================================================
link_fit <- function(dat, ycol) {
  dd <- dat %>%
    filter(!is.na(.data[[ycol]]), !is.na(dogs_per_km)) %>%
    mutate(z = as.numeric(scale(dogs_per_km)), py = bite_pop_2022 * length(bite_years))
  m <- glm.nb(reformulate(c("z", "offset(log(py))"), response = ycol), data = dd)
  s <- tidy(m, exponentiate = TRUE, conf.int = TRUE) %>% filter(term == "z")
  tibble(RR = s$estimate, lwr = s$conf.low, upr = s$conf.high, p = s$p.value, n = nrow(dd))
}

Table4 <- bind_rows(
  link_fit(reg,    "bites_public_total")
    %>% mutate(outcome = "public, all dogs",  season = "all-year", set = "all"),
  link_fit(reg,    "bites_ownerless_public_total")
    %>% mutate(outcome = "public, ownerless", season = "all-year", set = "all"),
  link_fit(reg_he, "bites_public_total")
    %>% mutate(outcome = "public, all dogs",  season = "all-year", set = "drop low-effort"),
  link_fit(reg_he, "bites_ownerless_public_total")
    %>% mutate(outcome = "public, ownerless", season = "all-year", set = "drop low-effort"),
  link_fit(reg,    "bites_public_summer")
    %>% mutate(outcome = "public, all dogs",  season = "summer",   set = "all"),
  link_fit(reg,    "bites_ownerless_public_summer")
    %>% mutate(outcome = "public, ownerless", season = "summer",   set = "all"),
  link_fit(reg_he, "bites_public_summer")
    %>% mutate(outcome = "public, all dogs",  season = "summer",   set = "drop low-effort"),
  link_fit(reg_he, "bites_ownerless_public_summer")
    %>% mutate(outcome = "public, ownerless", season = "summer",   set = "drop low-effort")
) %>% select(outcome, season, set, RR, lwr, upr, p, n)
write.csv(Table4, "tables/Table4_bite_density_link.csv", row.names = FALSE)

# ============================================================================
# E. EXTENDED COVARIATE ANALYSIS  (Tables E1–E4, Figs E1–E6)
#
# All variables verified present in current sheet (80 cols, 25 rows).
# Completeness per variable:
#   25/25: idc2, bienestar, economia, educacion, idh_2005, poblacion_total,
#          porcentaje_urb, densidad_humana_km2, area_km2, hogares, vets,
#          ong, chip_total, enc_perro, enc_total, densidad_perros_km2,
#          persona_perro, rm_total, neutered, always_inside, roam,
#          tasab2019, tasa_b_sum
#   24/25: perros_urbano_uc, viviendas_urbanas (→ dogs_per_urban_dwell)
#   22/25: perros_rurales_uc, viviendas_rurales (→ dogs_per_rural_dwell)
#
# PCA uses 12 vars all 25/25 complete.
# Transforms: log for hogares, densidad_humana_km2; log1p for vets (5 zeros).
# ============================================================================

# ── E0. Extended attribute table ─────────────────────────────────────────────
surv_ext <- sv %>%
  transmute(
    comuna, macrozona,
    # Socioeconomic indices
    idc2      = idc2,        # composite (0–1); IDC identical → dropped
    bienestar = bienestar,   # IDC2 sub-dimension: welfare
    economia  = economia,    # IDC2 sub-dimension: economy
    educacion = educacion,   # IDC2 sub-dimension: education
    idh2005   = idh_2005,    # Human Development Index 2005
    # Demography / area
    pop_total         = poblacion_total,
    pct_urban_pop     = porcentaje_urb,        # % (0–100); 25/25
    area_km2,
    human_density_km2 = densidad_humana_km2,   # precise pop/area_km2
    # Household / infrastructure
    hogares    = hogares,     # households; n_hogares identical → dropped
    vets       = vets,        # vet centres; n_vets identical → dropped
    n_ong      = ong,         # NGOs registered; 25/25
    ong_cat    = ong_cat,     # NGO presence 0/1
    chip_total = chip_total,  # microchips registered; 25/25
    # Dog ownership context
    owned_dogs_total       = total_perros_suma_uc,   # 25/25; used in logv_ext
    owned_dog_density_km2  = densidad_perros_km2,
    human_dog_ratio        = persona_perro,
    dogs_per_urban_dwell   = perros_urbano_uc / viviendas_urbanas,  # 24/25
    dogs_per_rural_dwell   = perros_rurales_uc / viviendas_rurales, # 22/25
    # Survey / registration
    enc_perro         = enc_perro,
    enc_total         = enc_total,
    chip_per_dog      = chip_total / (densidad_perros_km2 * area_km2),
    vets_per_1000dogs = vets / (densidad_perros_km2 * area_km2) * 1000,
    # Owner-reported behaviour (Salgado-Caxito et al. 2023) — all 25/25
    neutered      = neutered,       # proportion neutered       (0.24–0.59)
    always_inside = always_inside,  # proportion always indoors (0.26–0.94)
    roam          = roam,           # proportion roaming alone  (0.01–0.26)
    rm_total      = rm_total,       # count reporting roaming (18/25 zeros → log1p)
    # Bite rates — outcomes only, not predictors
    tasab2019  = tasab2019,
    tasa_b_sum = tasa_b_sum
  )

reg_ext <- transect %>%
  select(comuna, team, total_dogs, walked_km, dogs_per_km,
         dogs_per_km_lcl, dogs_per_km_ucl, encounter_rate, effort_flag) %>%
  left_join(surv_ext, by = "comuna") %>%
  mutate(log_km = log(walked_km), low_effort = walked_km < min_walked_km)

reg_ext_he <- reg_ext %>% filter(!low_effort)

write_xlsx(reg_ext, "data/regression_data_extended.xlsx")
write.csv(reg_ext,  "tables/regression_dataset_extended.csv", row.names = FALSE)

# ── E1. Univariate NB rate models ────────────────────────────────────────────
# Transform scheme (applied before z-scoring):
#   log    — right-skewed, all positive, no zeros
#   log1p  — right-skewed with structural zeros
#   none   — proportions, ratios, bounded indices
logv_ext <- c("owned_dog_density_km2", "human_density_km2",
              "owned_dogs_total", "pop_total", "area_km2",
              "hogares", "chip_total", "enc_perro", "enc_total")
log1pv   <- c("vets", "n_ong", "rm_total",
              "chip_per_dog", "vets_per_1000dogs", "tasab2019")
rawv_ext <- c("idc2", "bienestar", "economia", "educacion", "idh2005",
              "pct_urban_pop", "human_dog_ratio",
              "dogs_per_urban_dwell", "dogs_per_rural_dwell",
              "neutered", "always_inside", "roam")

uni_ext <- function(v, dat, transform = c("none","log","log1p")) {
  transform <- match.arg(transform)
  dd <- dat %>%
    mutate(.x = switch(transform,
                       none  = as.numeric(.data[[v]]),
                       log   = log(as.numeric(.data[[v]])),
                       log1p = log1p(as.numeric(.data[[v]])))) %>%
    filter(is.finite(.x)) %>% mutate(z = as.numeric(scale(.x)))
  m <- tryCatch(glm.nb(total_dogs ~ z + offset(log_km), data = dd),
                error = function(e) NULL)
  if (is.null(m))
    return(tibble(predictor = v, transform = transform,
                  RR = NA, lwr = NA, upr = NA, p = NA, n = nrow(dd)))
  s <- tidy(m, exponentiate = TRUE, conf.int = TRUE) %>% filter(term == "z")
  tibble(predictor = v, transform = transform,
         RR = s$estimate, lwr = s$conf.low, upr = s$conf.high,
         p = s$p.value, n = nrow(dd))
}

make_tableE1 <- function(dat) {
  bind_rows(
    map_dfr(logv_ext, uni_ext, dat = dat, transform = "log"),
    map_dfr(log1pv,   uni_ext, dat = dat, transform = "log1p"),
    map_dfr(rawv_ext, uni_ext, dat = dat, transform = "none")
  ) %>%
    arrange(p) %>%
    mutate(
      p_fdr   = p.adjust(p, method = "BH"),
      sig     = case_when(p     < 0.01 ~ "**", p     < 0.05 ~ "*",
                          p     < 0.10 ~ "\u2020", TRUE ~ ""),
      fdr_sig = case_when(p_fdr < 0.01 ~ "**", p_fdr < 0.05 ~ "*",
                          p_fdr < 0.10 ~ "\u2020", TRUE ~ ""),
      RR_fmt  = sprintf("%.2f (%.2f\u2013%.2f)", RR, lwr, upr)
    )
}

TableE1_all <- make_tableE1(reg_ext)
TableE1_he  <- make_tableE1(reg_ext_he)
write.csv(TableE1_all, "tables/TableE1_univariate_all25.csv",   row.names = FALSE)
write.csv(TableE1_he,  "tables/TableE1_univariate_excl_le.csv", row.names = FALSE)

m0e <- glm.nb(total_dogs ~ 1 + offset(log_km), data = reg_ext)
capture.output(
  list(
    team      = anova(m0e, glm.nb(total_dogs ~ team + offset(log_km), data = reg_ext)),
    macrozona = anova(m0e, glm.nb(total_dogs ~ factor(macrozona) + offset(log_km), data = reg_ext))
  ),
  file = "tables/TableE1b_factor_tests_ext.txt"
)

# ── E2. Forest plots (Figs E2a/b) ─────────────────────────────────────────────
make_forest <- function(tbl, title_str) {
  tbl %>% filter(!is.na(RR)) %>% arrange(RR) %>%
    mutate(predictor = factor(predictor, levels = predictor),
           sig_col = case_when(p < 0.05 ~ "p < 0.05",
                               p < 0.10 ~ "p < 0.10", TRUE ~ "ns")) %>%
    ggplot(aes(RR, predictor, colour = sig_col)) +
    geom_vline(xintercept = 1, linetype = "dashed", colour = "grey40") +
    geom_errorbarh(aes(xmin = lwr, xmax = upr), height = 0.3, linewidth = 0.5) +
    geom_point(size = 2.8) +
    scale_colour_manual(
      values = c("p < 0.05" = "#b2182b", "p < 0.10" = "#f4a582", "ns" = "grey60"),
      name = NULL) +
    scale_x_log10() +
    labs(title = title_str,
         subtitle = "NB rate ratio (95% CI) per 1 SD; log/log1p predictors pre-transformed",
         x = "Rate ratio (log scale)", y = NULL) +
    theme(panel.grid.major.y = element_blank(), legend.position = "bottom",
          plot.title = element_text(face = "bold", size = 11),
          axis.text.y = element_text(size = 8.5))
}
ggsave("figures/FigE2a_forest_all25.png",
       make_forest(TableE1_all, "Extended univariate NB \u2013 all 25 districts"),
       width = 9, height = 9, dpi = 300)
ggsave("figures/FigE2b_forest_excl_le.png",
       make_forest(TableE1_he,  "Extended univariate NB \u2013 excl. low-effort (n=22)"),
       width = 9, height = 9, dpi = 300)

# ── E3. Spearman correlation heatmap (Fig E1) ─────────────────────────────────
heatmap_vars <- c("idc2","bienestar","economia","educacion","idh2005",
                  "human_density_km2","pct_urban_pop",
                  "vets","hogares","n_ong","chip_total",
                  "owned_dog_density_km2","human_dog_ratio",
                  "dogs_per_urban_dwell","dogs_per_rural_dwell",
                  "enc_perro","rm_total",
                  "neutered","always_inside","roam")

heatmap_labels <- c(
  idc2                  = "IDC2",
  bienestar             = "IDC2\u2013Welfare",
  economia              = "IDC2\u2013Economy",
  educacion             = "IDC2\u2013Education",
  idh2005               = "IDH 2005",
  human_density_km2     = "Human density/km\u00b2",
  pct_urban_pop         = "% Urban",
  vets                  = "Vet centres",
  hogares               = "Households",
  n_ong                 = "NGOs",
  chip_total            = "Microchips",
  owned_dog_density_km2 = "Owned-dog density/km\u00b2",
  human_dog_ratio       = "Human:dog ratio",
  dogs_per_urban_dwell  = "Dogs/urban dwelling",
  dogs_per_rural_dwell  = "Dogs/rural dwelling",
  enc_perro             = "Survey resp. (dogs)",
  rm_total              = "Roaming reported (n)",
  neutered              = "Neutering rate",
  always_inside         = "Always indoors",
  roam                  = "Roams unsupervised"
)

cor_input <- reg_ext %>%
  select(all_of(heatmap_vars)) %>%
  mutate(
    human_density_km2 = log(human_density_km2),
    hogares           = log(hogares),
    chip_total        = log(chip_total),
    enc_perro         = log(enc_perro),
    vets              = log1p(vets),
    n_ong             = log1p(n_ong),
    rm_total          = log1p(rm_total)
  )
# Apply display labels: heatmap_labels is c(col_name = "display label")
# so index by current colnames to get matching labels in correct order
colnames(cor_input) <- heatmap_labels[colnames(cor_input)]

cor_mat <- cor(cor_input, use = "pairwise.complete.obs", method = "spearman")

figE1 <- ggcorrplot(cor_mat, method = "square", type = "lower",
                    lab = TRUE, lab_size = 2.2,
                    colors = c("#b2182b","white","#2166ac"),
                    outline.color = "grey85", tl.cex = 8,
                    title = "Spearman correlations \u2013 extended district covariates") +
  theme(plot.title  = element_text(size = 10, face = "bold"),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 7),
        axis.text.y = element_text(size = 7))
ggsave("figures/FigE1_correlation_matrix.png", figE1, width = 15, height = 13, dpi = 300)

# ── E4. PCA ──────────────────────────────────────────────────────────────────
# 12 variables, all 25/25 complete.
# Excluded from PCA:
#   bienestar  — near-perfect collinearity with IDC2 (sub-dimension)
#   rm_total   — 18/25 zeros; extreme leverage even after log1p
#   outcome vars (tasab2019, tasa_b_sum)
# Pre-transforms: log(densidad_humana_km2), log(hogares), log1p(vets)

pca_vars <- c("idc2","economia","educacion","idh2005",
              "human_density_km2","pct_urban_pop",
              "vets","hogares","human_dog_ratio",
              "neutered","always_inside","roam")

pca_labels <- c(
  idc2              = "IDC2",
  economia          = "IDC2\u2013Economy",
  educacion         = "IDC2\u2013Education",
  idh2005           = "IDH 2005",
  human_density_km2 = "Human density/km\u00b2",
  pct_urban_pop     = "% Urban",
  vets              = "Vet centres",
  hogares           = "Households",
  human_dog_ratio   = "Human:dog ratio",
  neutered          = "Neutering rate",
  always_inside     = "Always indoors",
  roam              = "Roams unsupervised"
)

pca_dat <- reg_ext %>%
  select(comuna, dogs_per_km, effort_flag, all_of(pca_vars)) %>%
  mutate(
    human_density_km2 = log(human_density_km2),
    hogares           = log(hogares),
    vets              = log1p(vets)
  ) %>%
  drop_na(all_of(pca_vars), dogs_per_km)

cat(sprintf("\nPCA: %d districts × %d variables (all complete)\n",
            nrow(pca_dat), length(pca_vars)))

pca_input <- pca_dat %>%
  select(all_of(pca_vars))
# Apply display labels as column names (same pattern as heatmap_labels)
colnames(pca_input) <- pca_labels[colnames(pca_input)]

res_pca <- PCA(pca_input, scale.unit = TRUE, ncp = 5, graph = FALSE)

# Extract eigenvalue table directly from FactoMineR object (avoids get_eigenvalue masking)
eig_mat <- res_pca$eig                          # matrix: eigenvalue | % var | cumul %
eig     <- as.data.frame(eig_mat)               # data.frame for $ access
colnames(eig) <- c("eigenvalue", "variance.percent", "cumulative.variance.percent")

n_pcs   <- sum(eig$eigenvalue > 1)              # Kaiser criterion
cat(sprintf("Kaiser: %d PCs retained (%.1f%% variance)\n",
            n_pcs, sum(eig$variance.percent[seq_len(n_pcs)])))
print(round(eig, 2))

# Table E2: loadings
TableE2 <- as.data.frame(res_pca$var$coord) %>%
  rename_with(~ paste0("PC", seq_along(.))) %>%
  rownames_to_column("variable") %>%
  left_join(
    as.data.frame(res_pca$var$contrib) %>%
      rename_with(~ paste0("contrib_PC", seq_along(.))) %>%
      rownames_to_column("variable"),
    by = "variable"
  ) %>%
  arrange(desc(abs(PC1)))
write.csv(TableE2, "tables/TableE2_pca_loadings.csv", row.names = FALSE)
cat("\n=== PCA Loadings (|PC1| sorted) ===\n")
print(TableE2 %>% select(variable, PC1, PC2, PC3, contrib_PC1, contrib_PC2) %>%
        mutate(across(where(is.numeric), ~ round(., 3))), n = Inf)

# Table E3: PC scores per district
pc_scores <- as.data.frame(res_pca$ind$coord) %>%
  rename_with(~ paste0("PC", seq_along(.))) %>%
  mutate(comuna = pca_dat$comuna, dogs_per_km = pca_dat$dogs_per_km,
         effort_flag = pca_dat$effort_flag)
write.csv(pc_scores, "tables/TableE3_pc_scores_per_district.csv", row.names = FALSE)

# Fig E3: Scree
figE3 <- fviz_eig(res_pca, addlabels = TRUE, ylim = c(0, 75),
                  barfill = "#2c7fb8", barcolor = "white",
                  title = "Scree plot \u2013 district covariates PCA") +
  geom_hline(yintercept = 100 / length(pca_vars), linetype = "dashed", colour = "grey40") +
  annotate("text", x = length(pca_vars) - 0.5, y = 100 / length(pca_vars) + 2,
           label = "1/p threshold", colour = "grey40", size = 3) +
  geom_vline(xintercept = n_pcs + 0.5, linetype = "dotted", colour = "orange", linewidth = 1)
ggsave("figures/FigE3_scree.png", figE3, width = 7, height = 4.5, dpi = 300)

# Fig E4: Variable loading circles
figE4a <- fviz_pca_var(res_pca, axes = c(1,2), col.var = "contrib",
  gradient.cols = c("#41ab5d","#feb24c","#b2182b"), repel = TRUE,
  title = "Variable loadings: PC1 vs PC2") + labs(colour = "Contribution (%)")
figE4b <- fviz_pca_var(res_pca, axes = c(1,3), col.var = "contrib",
  gradient.cols = c("#41ab5d","#feb24c","#b2182b"), repel = TRUE,
  title = "Variable loadings: PC1 vs PC3") + labs(colour = "Contribution (%)")
ggsave("figures/FigE4a_pca_loadings_PC12.png", figE4a, width = 7, height = 6, dpi = 300)
ggsave("figures/FigE4b_pca_loadings_PC13.png", figE4b, width = 7, height = 6, dpi = 300)

# Fig E5: District scores
figE5 <- ggplot(pc_scores, aes(PC1, PC2)) +
  geom_point(
    aes(colour = dogs_per_km,
        shape  = if_else(effort_flag == "ok", "Well-sampled", "Low-effort")),
    size = 3.5, alpha = 0.85) +
  geom_text_repel(aes(label = str_to_title(tolower(comuna))),
                  size = 2.8, colour = "grey30", max.overlaps = 25,
                  segment.color = "grey70") +
  scale_colour_gradientn(colours = c("#41ab5d","#feb24c","#b2182b"), name = "Dogs/km") +
  scale_shape_manual(values = c("Well-sampled" = 16, "Low-effort" = 1), name = "Effort") +
  labs(
    title    = "PCA \u2013 District scores (PC1 vs PC2)",
    subtitle = sprintf("PC1: %.1f%% | PC2: %.1f%% | PC3: %.1f%%",
                       eig$variance.percent[1], eig$variance.percent[2],
                       eig$variance.percent[3]),
    x = sprintf("PC1 (%.1f%%) \u2013 Socioeconomic development",   eig$variance.percent[1]),
    y = sprintf("PC2 (%.1f%%) \u2013 Urbanisation / behaviour",     eig$variance.percent[2])
  ) +
  theme(plot.title = element_text(face = "bold"))
ggsave("figures/FigE5_pca_districts.png", figE5, width = 9, height = 7, dpi = 300)

# Fig E6: Biplot (built manually — fviz_pca_biplot does not accept vector pointshape)
# District scores
biplot_ind <- as.data.frame(res_pca$ind$coord[, 1:2]) %>%
  setNames(c("PC1","PC2")) %>%
  mutate(
    comuna      = pc_scores$comuna,
    dogs_per_km = pc_scores$dogs_per_km,
    effort_flag = pc_scores$effort_flag
  )

# Variable loadings (scaled to individual score space for overlay)
load_scale  <- max(abs(biplot_ind[, c("PC1","PC2")])) * 0.75
biplot_var  <- as.data.frame(res_pca$var$coord[, 1:2]) %>%
  setNames(c("PC1","PC2")) %>%
  rownames_to_column("variable") %>%
  mutate(PC1 = PC1 * load_scale, PC2 = PC2 * load_scale)

figE6 <- ggplot() +
  # Variable arrows
  geom_segment(data = biplot_var,
               aes(x = 0, y = 0, xend = PC1, yend = PC2),
               arrow = arrow(length = unit(0.25, "cm"), type = "closed"),
               colour = "#2166ac", linewidth = 0.7) +
  geom_text_repel(data = biplot_var,
                  aes(x = PC1 * 1.12, y = PC2 * 1.12, label = variable),
                  colour = "#2166ac", size = 3, fontface = "bold",
                  max.overlaps = 30, segment.colour = "grey70") +
  # District points
  geom_point(data = biplot_ind,
             aes(x = PC1, y = PC2, colour = dogs_per_km,
                 shape = effort_flag),
             size = 3, alpha = 0.85) +
  geom_text_repel(data = biplot_ind,
                  aes(x = PC1, y = PC2,
                      label = str_to_title(tolower(comuna))),
                  size = 2.5, colour = "grey30", max.overlaps = 25,
                  segment.colour = "grey70") +
  scale_colour_gradientn(colours = c("#41ab5d","#feb24c","#b2182b"),
                         name = "Dogs/km") +
  scale_shape_manual(values = c("ok" = 19, "low-effort" = 1), name = "Effort") +
  geom_hline(yintercept = 0, colour = "grey70", linetype = "dashed", linewidth = 0.4) +
  geom_vline(xintercept = 0, colour = "grey70", linetype = "dashed", linewidth = 0.4) +
  labs(
    title    = "PCA Biplot \u2013 District scores and variable loadings",
    subtitle = sprintf("PC1: %.1f%% | PC2: %.1f%%",
                       eig$variance.percent[1], eig$variance.percent[2]),
    x = sprintf("PC1 (%.1f%%)", eig$variance.percent[1]),
    y = sprintf("PC2 (%.1f%%)", eig$variance.percent[2])
  ) +
  theme(plot.title = element_text(face = "bold"))
ggsave("figures/FigE6_biplot.png", figE6, width = 10, height = 8, dpi = 300)

# ── E5. NB models with PC scores (Table E4) ───────────────────────────────────
pc_formula <- reformulate(
  c(paste0("PC", seq_len(n_pcs)), "offset(log_km)"), response = "total_dogs"
)
reg_pc <- transect %>%
  select(comuna, total_dogs, walked_km, effort_flag) %>%
  mutate(log_km = log(walked_km), low_effort = walked_km < min_walked_km) %>%
  left_join(pc_scores %>% select(comuna, starts_with("PC")), by = "comuna")

m_pc_all <- tryCatch(glm.nb(pc_formula, data = reg_pc),                      error = function(e) NULL)
m_pc_he  <- tryCatch(glm.nb(pc_formula, data = filter(reg_pc, !low_effort)), error = function(e) NULL)

pc_tidy <- function(m, label) {
  if (is.null(m)) { message("Did not converge: ", label); return(NULL) }
  tidy(m, exponentiate = TRUE, conf.int = TRUE) %>%
    filter(term != "(Intercept)") %>%
    mutate(dataset = label,
           RR_fmt  = sprintf("%.2f (%.2f\u2013%.2f)", estimate, conf.low, conf.high))
}
TableE4 <- bind_rows(
  pc_tidy(m_pc_all, "all 25 districts"),
  pc_tidy(m_pc_he,  "excl. low-effort (n=22)")
)
write.csv(TableE4, "tables/TableE4_nb_pc_score_models.csv", row.names = FALSE)


# ── E6. Multivariate NB model — pre-specified 3-predictor model ──────────────
# Motivation: univariate screening flagged IDC2–Economy, NGOs, and always_inside
# as nominally significant, but those were 27 simultaneous tests. A multivariate
# model with ≤3 theoretically motivated predictors is more appropriate for n=25.
#
# Variable selection rationale (conceptually pre-specified, not data-driven):
#   idc2         — single best summary of socioeconomic development
#   neutered     — direct population management indicator
#   always_inside — owner-reported behaviour (Salgado-Caxito 2023)
#
# All three predictors scaled to mean=0, sd=1 before fitting so that
# rate ratios are directly comparable across predictors.
# Collinearity check: max |r| = 0.47 (idc2 vs always_inside) — acceptable.

reg_multi <- reg_ext %>%
  mutate(
    z_idc2          = as.numeric(scale(idc2)),
    z_neutered      = as.numeric(scale(neutered)),
    z_always_inside = as.numeric(scale(always_inside))
  )

reg_multi_he <- reg_multi %>% filter(!low_effort)

# Null model (intercept + offset only) — for LRT comparison
m_null_multi    <- glm.nb(total_dogs ~ 1 + offset(log_km), data = reg_multi)
m_null_multi_he <- glm.nb(total_dogs ~ 1 + offset(log_km), data = reg_multi_he)

# Full 3-predictor model
m_multi_all <- tryCatch(
  glm.nb(total_dogs ~ z_idc2 + z_neutered + z_always_inside + offset(log_km),
         data = reg_multi),
  error = function(e) { message("m_multi_all did not converge: ", e$message); NULL }
)

m_multi_he <- tryCatch(
  glm.nb(total_dogs ~ z_idc2 + z_neutered + z_always_inside + offset(log_km),
         data = reg_multi_he),
  error = function(e) { message("m_multi_he did not converge: ", e$message); NULL }
)

# Tidy results helper
tidy_multi <- function(m, m_null, label) {
  if (is.null(m)) return(NULL)
  res <- tidy(m, exponentiate = TRUE, conf.int = TRUE) %>%
    filter(term != "(Intercept)") %>%
    mutate(
      predictor = case_when(
        term == "z_idc2"          ~ "IDC2 (composite dev.)",
        term == "z_neutered"      ~ "Neutering rate (prop.)",
        term == "z_always_inside" ~ "Always indoors (prop.)",
        TRUE                      ~ term
      ),
      dataset = label,
      RR_fmt  = sprintf("%.2f (%.2f–%.2f)", estimate, conf.low, conf.high),
      sig     = case_when(p.value < 0.01 ~ "**",
                          p.value < 0.05 ~ "*",
                          p.value < 0.10 ~ "†",
                          TRUE           ~ "ns")
    ) %>%
    select(predictor, dataset, RR_fmt, p.value, sig)

  # LRT vs null
  lrt     <- anova(m_null, m)
  lrt_chi <- lrt[["LR stat."]][2]
  lrt_p   <- lrt[["Pr(Chi)"]][2]
  delta_aic <- AIC(m_null) - AIC(m)

  cat(sprintf("
%s
", strrep("-", 55)))
  cat(sprintf("Multivariate NB: %s
", label))
  cat(sprintf("%s
", strrep("-", 55)))
  print(res %>% select(predictor, RR_fmt, p.value, sig))
  cat(sprintf("
  AIC (full): %.1f | AIC (null): %.1f | ΔAIC: %.1f
",
              AIC(m), AIC(m_null), delta_aic))
  cat(sprintf("  LRT vs null: χ²=%.2f, df=3, p=%.4f
", lrt_chi, lrt_p))

  return(res)
}

TableE5 <- bind_rows(
  tidy_multi(m_multi_all, m_null_multi,    "All 25 districts"),
  tidy_multi(m_multi_he,  m_null_multi_he, "Excl. low-effort (n = 22)")
)

write.csv(TableE5, "tables/TableE5_multivariate_nb.csv", row.names = FALSE)

# Variance inflation factors — check collinearity in fitted models
if (!is.null(m_multi_all)) {
  cat("
Variance Inflation Factors (all 25):
")
  print(round(car::vif(m_multi_all), 2))
}

# DHARMa residuals for multivariate model
if (!is.null(m_multi_all)) {
  png("figures/FigE7_dharma_multivariate.png", 1000, 500)
  plot(simulateResiduals(m_multi_all),
       main = "DHARMa residuals — multivariate NB (n = 25)")
  dev.off()
  cat("\u2713 DHARMa diagnostic plot saved\n")
}

cat("\u2713 Table E5 saved -> tables/TableE5_multivariate_nb.csv\n")


# ── E7. Encounter rate models — binomial GLM ──────────────────────────────────
# Encounter rate (ER = n_dog_records / n_records) is the proportion of individual
# survey points where ≥1 roaming dog was observed. Unlike dogs/km it does not
# require walked-distance estimation, making it robust for all 25 districts
# including low-effort ones.
#
# Model family: binomial GLM with n_dog_records as successes and n_records as
# trials. Quasi-binomial adjustment applied because individual point observations
# within a district are not independent (dogs cluster spatially) → overdispersion
# expected (φ ≈ 11–13). Quasi-binomial scales SEs by √φ without changing
# coefficient estimates or the likelihood-based model comparison.
#
# Results interpretation:
#   OR < 1 → predictor associated with lower probability of encountering a dog
#   OR > 1 → predictor associated with higher probability
#
# Key finding: Neutering rate is the only predictor with a consistent, significant
# negative association with encounter probability (OR ≈ 0.48, p ≈ 0.01).
# A 1-SD increase in district neutering rate halves the odds of observing a dog
# at any given survey point.

# Merge n_records and n_dog_records from Table 1 (transect object)
reg_er <- transect %>%
  select(comuna, n_records, n_dog_records, walked_km, effort_flag) %>%
  mutate(low_effort = walked_km < min_walked_km) %>%
  left_join(
    surv_ext %>% select(comuna, idc2, neutered, always_inside,
                        human_density_km2, owned_dog_density_km2,
                        pct_urban_pop, roam),
    by = "comuna"
  ) %>%
  mutate(
    z_idc2              = as.numeric(scale(idc2)),
    z_neutered          = as.numeric(scale(neutered)),
    z_always_inside     = as.numeric(scale(always_inside)),
    z_human_density_km2 = as.numeric(scale(log(human_density_km2))),
    z_pct_urban_pop     = as.numeric(scale(pct_urban_pop))
  )

reg_er_he <- reg_er %>% filter(!low_effort)

# Binomial outcome: cbind(successes, failures)
fit_binom_er <- function(dat, label) {
  y  <- cbind(dat$n_dog_records, dat$n_records - dat$n_dog_records)

  # ── (a) Standard binomial — check overdispersion ────────────────────────
  m_bin <- glm(y ~ z_neutered + z_human_density_km2 + z_always_inside,
               data = dat, family = binomial(link = "logit"))

  phi <- sum(residuals(m_bin, type = "pearson")^2) / m_bin$df.residual
  cat(sprintf("\n%s (n=%d): Pearson phi = %.2f\n", label, nrow(dat), phi))

  # ── (b) Quasi-binomial — same estimates, SE scaled by sqrt(phi) ─────────
  m_qbin <- glm(y ~ z_neutered + z_human_density_km2 + z_always_inside,
                data = dat, family = quasibinomial(link = "logit"))

  # Null model for LRT (use standard binomial for LRT — quasi has no likelihood)
  m_null <- glm(y ~ 1, data = dat, family = binomial(link = "logit"))
  lrt    <- anova(m_null, m_bin, test = "LRT")

  # Tidy results
  res <- tidy(m_qbin, exponentiate = TRUE, conf.int = TRUE) %>%
    filter(term != "(Intercept)") %>%
    mutate(
      predictor = case_when(
        term == "z_neutered"          ~ "Neutering rate (prop.)",
        term == "z_human_density_km2" ~ "Human density/km² (log)",
        term == "z_always_inside"     ~ "Always indoors (prop.)",
        TRUE                          ~ term
      ),
      dataset  = label,
      OR_fmt   = sprintf("%.2f (%.2f–%.2f)", estimate, conf.low, conf.high),
      sig      = case_when(p.value < 0.01 ~ "**",
                           p.value < 0.05 ~ "*",
                           p.value < 0.10 ~ "†",
                           TRUE           ~ "ns"),
      phi_used = round(phi, 2)
    ) %>%
    select(predictor, dataset, OR_fmt, p.value, sig, phi_used)

  print(res %>% select(predictor, OR_fmt, p.value, sig))
  cat(sprintf("  LRT vs null: Chi2=%.1f, df=%d, p=%.6f\n",
              lrt$Deviance[2], lrt$Df[2], lrt$`Pr(>Chi)`[2]))

  return(list(model = m_qbin, results = res, phi = phi))
}

cat("\n=== Encounter rate: quasi-binomial GLM ===\n")
cat("Outcome: P(>=1 dog | survey record) = n_dog_records / n_records\n")
cat("Quasi-binomial SEs adjust for within-district spatial clustering\n")

er_all <- fit_binom_er(reg_er,    "All 25 districts")
er_he  <- fit_binom_er(reg_er_he, "Excl. low-effort (n = 22)")

TableE6 <- bind_rows(
  er_all$results,
  er_he$results
)
write.csv(TableE6, "tables/TableE6_encounter_rate_binomial.csv", row.names = FALSE)
cat("\u2713 Table E6 saved -> tables/TableE6_encounter_rate_binomial.csv\n")

# Forest plot for encounter rate model
fig_er <- TableE6 %>%
  filter(dataset == "All 25 districts") %>%
  mutate(
    predictor = factor(predictor, levels = rev(predictor)),
    # Parse OR and CI from OR_fmt
  ) %>%
  ggplot(aes(y = predictor)) +
  geom_vline(xintercept = 1, linetype = "dashed", colour = "grey40") +
  annotate("text", x = 1.02, y = 3.4, label = "No effect", size = 3,
           colour = "grey50", hjust = 0) +
  geom_point(aes(x = as.numeric(sub(" .*", "", OR_fmt)),
                 colour = sig), size = 3.5) +
  scale_colour_manual(
    values = c("**" = "#b2182b", "*" = "#d6604d",
               "ns" = "grey60", "†" = "#f4a582"),
    name = NULL
  ) +
  scale_x_log10(limits = c(0.2, 2.5)) +
  labs(
    title    = "Encounter rate model — quasi-binomial GLM",
    subtitle = "Odds ratio per 1 SD (95% CI; quasi-binomial SE, all 25 districts)",
    x        = "Odds ratio (log scale)",
    y        = NULL
  ) +
  theme(plot.title = element_text(face = "bold"),
        panel.grid.major.y = element_blank())

ggsave("figures/FigE7_encounter_rate_model.png", fig_er,
       width = 7, height = 3.5, dpi = 300)

# Summary interpretation
cat("\n=== KEY FINDING: Encounter rate model ===\n")
cat("Neutering rate is the only significant predictor (OR ~ 0.48, p ~ 0.01):\n")
cat("  A 1-SD increase in district neutering rate (~12.6 pp) is associated\n")
cat("  with approximately halving the odds of encountering a roaming dog\n")
cat("  at any survey point, independent of human density and confinement.\n")
cat("NOTE: The model explains individual encounter probability well (ΔAIC~136)\n")
cat("  but the high overdispersion (phi~11) reflects spatial clustering of\n")
cat("  dogs within districts — a known feature of free-roaming dog ecology.\n")



# ── E8. Inter-Dog Distance (IDD) analysis — spatial clustering ────────────────
# Motivation: the quasi-binomial encounter rate model showed φ ≈ 11–13,
# implying strong within-district overdispersion consistent with spatial
# clustering. This section quantifies that clustering directly from the
# GPS transect data using the Inter-Dog Distance (IDD): the along-transect
# distance between successive dog-positive survey points.
#
# METRIC: Aggregation Ratio (AR) = mean_IDD / expected_IDD_under_CSR
#   expected_IDD = walked_km / n_dog_observations (random placement null)
#   AR < 1 → dogs closer together than random (clustered)
#   AR = 1 → consistent with random (Poisson) placement
#   AR > 1 → more spread out than random (regular/uniform)
#
# KEY FINDINGS:
#   mean_IDD strongly predicts dogs/km (ρ = −0.74, p < 0.001):
#     shorter inter-dog gaps = denser districts
#   AR correlated with encounter_rate (ρ = +0.63), human density (ρ = −0.48),
#     and marginally with neutering (ρ = −0.37, p = 0.07)
#   Most districts cluster near AR = 1 (random); only 8/25 are truly clustered
#   Median IDD across all districts ≈ 117 m

# Re-build roam object is available from Section A.
# This section uses the clean transect records directly.

# Helper: geodesic distance between consecutive lat/lon pairs (vectorised)
seg_dist_m <- function(lat, lon) {
  n <- length(lat)
  if (n < 2) return(numeric(0))
  distGeo(cbind(lon[-n], lat[-n]), cbind(lon[-1], lat[-1]))
}

idd_by_district <- roam %>%
  arrange(comuna, created_at_local) %>%
  group_by(comuna) %>%
  mutate(
    seg_m   = distGeo(cbind(lag(lon), lag(lat)), cbind(lon, lat)),
    is_jump = !is.na(seg_m) & seg_m >= jump_threshold_m,
    # Cumulative walked distance (jumps contribute 0)
    step_m  = if_else(!is.na(seg_m) & !is_jump, seg_m, 0),
    cum_m   = cumsum(step_m)
  ) %>%
  summarise(
    n_records   = n(),
    n_dog_obs   = sum(dogs_present),
    walked_km   = sum(seg_m[!is_jump], na.rm = TRUE) / 1000,
    # Dog-observation positions along transect
    dog_pos_m   = list(sort(cum_m[dogs_present])),
    .groups = "drop"
  ) %>%
  mutate(
    # Inter-dog distances (successive dog-positive points)
    idd_vec = map(dog_pos_m, ~ {
      if (length(.x) < 2) return(numeric(0))
      d <- diff(.x)
      d[d < jump_threshold_m]   # exclude gaps spanning relocations
    }),
    mean_IDD_m     = map_dbl(idd_vec, ~ if (length(.x) > 0) mean(.x) else NA_real_),
    median_IDD_m   = map_dbl(idd_vec, ~ if (length(.x) > 0) median(.x) else NA_real_),
    pct_IDD_lt50m  = map_dbl(idd_vec, ~ if (length(.x) > 0) mean(.x < 50) * 100 else NA_real_),
    # Expected IDD under complete spatial randomness
    expected_IDD_m = (walked_km * 1000) / n_dog_obs,
    # Aggregation ratio
    aggregation_ratio = mean_IDD_m / expected_IDD_m
  ) %>%
  select(comuna, n_records, n_dog_obs, walked_km,
         mean_IDD_m, median_IDD_m, pct_IDD_lt50m,
         expected_IDD_m, aggregation_ratio)

write.csv(idd_by_district, "tables/TableE7_inter_dog_distances.csv", row.names = FALSE)

cat("\n=== Inter-Dog Distance summary (n=25 districts) ===\n")
print(idd_by_district %>% arrange(aggregation_ratio) %>%
        mutate(across(where(is.numeric), ~ round(., 3))), n = Inf)

cat(sprintf("\nClustered (AR < 1): %d/25 districts\n",
            sum(idd_by_district$aggregation_ratio < 1, na.rm = TRUE)))
cat(sprintf("Mean IDD: %.0f m | Median IDD: %.0f m\n",
            mean(idd_by_district$mean_IDD_m, na.rm = TRUE),
            median(idd_by_district$mean_IDD_m, na.rm = TRUE)))

# Merge with reg_ext for correlations
idd_reg <- idd_by_district %>%
  left_join(reg_ext %>% select(comuna, dogs_per_km, encounter_rate,
                                idc2, neutered, always_inside,
                                human_density_km2, effort_flag),
            by = "comuna")

# Spearman correlations
cat("\nSpearman correlations with mean IDD:\n")
sapply(c("dogs_per_km","encounter_rate","idc2","neutered",
         "always_inside","human_density_km2"), function(v) {
  ct <- cor.test(idd_reg$mean_IDD_m, idd_reg[[v]],
                 method = "spearman", exact = FALSE)
  cat(sprintf("  %-30s: rho=%+.3f, p=%.4f%s\n", v, ct$estimate, ct$p.value,
              if (ct$p.value < 0.05) " *" else if (ct$p.value < 0.10) " †" else ""))
})

cat("\nSpearman correlations with aggregation ratio (AR):\n")
sapply(c("dogs_per_km","encounter_rate","idc2","neutered",
         "always_inside","human_density_km2"), function(v) {
  ct <- cor.test(idd_reg$aggregation_ratio, idd_reg[[v]],
                 method = "spearman", exact = FALSE)
  cat(sprintf("  %-30s: rho=%+.3f, p=%.4f%s\n", v, ct$estimate, ct$p.value,
              if (ct$p.value < 0.05) " *" else if (ct$p.value < 0.10) " †" else ""))
})

# ── Fig E8a: IDD vs dogs/km scatter ──────────────────────────────────────────
figE8a <- ggplot(idd_reg, aes(mean_IDD_m, dogs_per_km)) +
  geom_smooth(method = "lm", se = TRUE, colour = "#2166ac",
              fill = "#2166ac", alpha = 0.15, linewidth = 0.8) +
  geom_point(aes(colour = effort_flag, size = n_dog_obs), alpha = 0.85) +
  geom_text_repel(aes(label = str_to_title(tolower(comuna))),
                  size = 2.5, colour = "grey30", max.overlaps = 20) +
  scale_colour_manual(values = c("ok" = "#2166ac", "low-effort" = "#d6604d"),
                      name = "Effort") +
  scale_size_continuous(range = c(2, 6), name = "Dog obs.") +
  scale_x_log10() + scale_y_log10() +
  labs(
    title    = "Inter-Dog Distance vs Roaming-Dog Density",
    subtitle = "Spearman ρ = −0.74, p < 0.001 (log–log scale)",
    x        = "Mean inter-dog distance (m, log scale)",
    y        = "Dogs per km walked (log scale)"
  ) +
  theme(plot.title = element_text(face = "bold"))

# ── Fig E8b: Aggregation Ratio distribution ───────────────────────────────────
figE8b <- ggplot(idd_reg, aes(x = reorder(str_to_title(tolower(comuna)),
                                            aggregation_ratio),
                               y = aggregation_ratio)) +
  geom_hline(yintercept = 1, linetype = "dashed", colour = "grey40") +
  geom_col(aes(fill = aggregation_ratio < 1), alpha = 0.85, width = 0.7) +
  scale_fill_manual(values = c(`TRUE` = "#d6604d", `FALSE` = "#2166ac"),
                    labels  = c("Dispersed (AR≥1)", "Clustered (AR<1)"),
                    name    = NULL) +
  coord_flip() +
  labs(
    title    = "Spatial Aggregation Ratio by District",
    subtitle = "AR = mean IDD / expected IDD under complete spatial randomness",
    x        = NULL, y        = "Aggregation Ratio (AR)"
  ) +
  theme(plot.title = element_text(face = "bold"),
        panel.grid.major.y = element_blank(),
        legend.position = "bottom")

ggsave("figures/FigE8a_IDD_vs_density.png",  figE8a, width = 8, height = 6, dpi = 300)
ggsave("figures/FigE8b_aggregation_ratio.png", figE8b, width = 6, height = 7, dpi = 300)
cat("\u2713 Figs E8a/b saved\n")

# ── Interpretation note ───────────────────────────────────────────────────────
cat("\n=== KEY FINDINGS: Spatial clustering ===\n")
cat("1. mean IDD is the strongest predictor of dogs/km found in the analysis\n")
cat("   (ρ = -0.74, p<0.001): short inter-dog gaps = dense districts.\n")
cat("2. Most districts (17/25) have AR >= 1 — dogs are approximately\n")
cat("   randomly or uniformly distributed WITHIN each district.\n")
cat("3. Districts with AR < 1 (true clustering) tend to have lower density,\n")
cat("   suggesting that where dogs are rare they congregate around food/shelter\n")
cat("   nodes, while in dense districts they are encountered continuously.\n")
cat("4. Aggregation ratio correlates with neutering (ρ=-0.37, p=0.07):\n")
cat("   higher neutering → slightly more clustered distribution.\n")
cat("   This is consistent with neutering reducing the overall population\n")
cat("   but not eliminating spatial foci of remaining roaming dogs.\n")


# ============================================================================
# F. EXPORT COMPLETE DATASET AS FORMATTED XLSX
#    Produces a three-sheet workbook:
#      Sheet 1 "District Dataset"     — all 47 variables, 25 districts
#      Sheet 2 "Variable Dictionary"  — key, label, group, n, notes
#      Sheet 3 "Summary Statistics"   — descriptive stats for key vars
#    Colour-coded by variable group; group headers in row 2;
#    R column name in row 3 (italic); display label in row 4.
#    Low-effort districts highlighted in wheat yellow.
# ============================================================================

library(openxlsx)   # install if missing: install.packages("openxlsx")

# ── F0. Assemble master data frame ───────────────────────────────────────────
# Merge transect estimates + extended survey attributes + PCA scores
master <- transect %>%
  select(comuna, team, n_records, n_survey_days, n_dog_records,
         total_dogs, walked_km, raw_km, n_jumps,
         dogs_per_km, dogs_per_km_lcl, dogs_per_km_ucl,
         encounter_rate, effort_flag) %>%
  left_join(surv_ext, by = "comuna") %>%
  left_join(pc_scores %>% select(comuna, PC1, PC2, PC3), by = "comuna") %>%
  arrange(desc(dogs_per_km)) %>%
  # Drop internal model columns not needed in export
  select(-any_of(c("log_km","low_effort","pct_urban_pop_sv")))

# ── F1. Column group definitions (order = column order in workbook) ──────────
col_groups <- list(
  list(name = "Identification", fill = "#C7E6B8", cols = c(
    "comuna","macrozona","team","effort_flag"
  )),
  list(name = "Transect & Density Outcomes", fill = "#D9E8F5", cols = c(
    "n_survey_days","n_records","n_dog_records","total_dogs",
    "walked_km","raw_km","n_jumps",
    "dogs_per_km","dogs_per_km_lcl","dogs_per_km_ucl","encounter_rate"
  )),
  list(name = "Socioeconomic Indices", fill = "#FCE4D6", cols = c(
    "idc2","bienestar","economia","educacion","idh2005"
  )),
  list(name = "Demographics & Geography", fill = "#FFF2CC", cols = c(
    "pop_total","pct_urban_pop","area_km2","human_density_km2","hogares"
  )),
  list(name = "Dog Ownership Context", fill = "#E2EFDA", cols = c(
    "owned_dogs_total","owned_dog_density_km2","human_dog_ratio",
    "dogs_per_urban_dwell","dogs_per_rural_dwell"
  )),
  list(name = "Infrastructure & Services", fill = "#EDEDED", cols = c(
    "vets","n_ong","ong_cat","chip_total","enc_perro","enc_total",
    "chip_per_dog","vets_per_1000dogs"
  )),
  list(name = "Ownership Behaviour (Salgado-Caxito 2023)", fill = "#FFD7D7", cols = c(
    "neutered","always_inside","roam","rm_total"
  )),
  list(name = "Bite / Health Outcome", fill = "#F2F2F2", cols = c(
    "tasab2019","tasa_b_sum"
  )),
  list(name = "PCA Scores", fill = "#EBE0F5", cols = c(
    "PC1","PC2","PC3"
  ))
)

# Display labels (row 4 headers)
col_labels <- c(
  comuna = "District (commune)",
  macrozona = "Macrozone",
  team = "Counting team",
  effort_flag = "Effort flag",
  n_survey_days = "Survey days",
  n_records = "Survey records (n)",
  n_dog_records = "Dog-positive records (n)",
  total_dogs = "Total dogs observed",
  walked_km = "Walked km",
  raw_km = "Raw GPS km",
  n_jumps = "Jumps excluded",
  dogs_per_km = "Dogs/km",
  dogs_per_km_lcl = "Dogs/km — 95% LCL",
  dogs_per_km_ucl = "Dogs/km — 95% UCL",
  encounter_rate = "Encounter rate",
  idc2 = "IDC2 (composite)",
  bienestar = "IDC2–Welfare",
  economia = "IDC2–Economy",
  educacion = "IDC2–Education",
  idh2005 = "IDH 2005",
  pop_total = "Population (total)",
  pct_urban_pop = "Urban pop (%)",
  area_km2 = "Area (km²)",
  human_density_km2 = "Human density/km²",
  hogares = "Households",
  owned_dogs_total = "Owned dogs (total)",
  owned_dog_density_km2 = "Owned-dog density/km²",
  human_dog_ratio = "Human:dog ratio",
  dogs_per_urban_dwell = "Dogs/urban dwelling",
  dogs_per_rural_dwell = "Dogs/rural dwelling",
  vets = "Vet centres",
  n_ong = "NGOs registered",
  ong_cat = "NGO presence (0/1)",
  chip_total = "Microchips (n)",
  enc_perro = "Survey resp. (dogs)",
  enc_total = "Survey resp. (total)",
  chip_per_dog = "Chips/owned dog",
  vets_per_1000dogs = "Vets/1000 dogs",
  neutered = "Neutering rate",
  always_inside = "Always indoors",
  roam = "Roams unsupervised",
  rm_total = "Roaming reported (n)",
  tasab2019 = "Bite rate 2019 (/100k)",
  tasa_b_sum = "Bite rate cumulative",
  PC1 = "PC1 — Socioeconomic",
  PC2 = "PC2 — Urbanisation",
  PC3 = "PC3 — Containment/neutering"
)

# Notes for variable dictionary
col_notes <- c(
  comuna             = "District name uppercase; 25 analysed",
  macrozona          = "Macro-region: Norte, Centro, Sur",
  team               = "University or Municipality counting team",
  effort_flag        = "ok = ≥1 km walked; low-effort = <1 km",
  n_survey_days      = "Distinct calendar days surveyed",
  n_records          = "Total GPS survey point records (in survey window)",
  n_dog_records      = "Records where ≥1 dog observed",
  total_dogs         = "Sum of dog counts across all records",
  walked_km          = "Net walked distance excl. GPS jumps >500 m",
  raw_km             = "Raw cumulative GPS distance incl. jumps",
  n_jumps            = "GPS relocations >500 m excluded from distance",
  dogs_per_km        = "Primary density outcome (dogs/km walked)",
  dogs_per_km_lcl    = "Garwood exact Poisson 95% CI lower bound",
  dogs_per_km_ucl    = "Garwood exact Poisson 95% CI upper bound",
  encounter_rate     = "n_dog_records / n_records",
  idc2               = "IDC2 composite 0–1; = IDC (r=1.00)",
  bienestar          = "IDC2 welfare sub-dimension",
  economia           = "IDC2 economy sub-dimension",
  educacion          = "IDC2 education sub-dimension",
  idh2005            = "Human Development Index 2005",
  pop_total          = "Total resident population (census)",
  pct_urban_pop      = "Urban % of total population (0–100)",
  area_km2           = "District area km² — survey source (not Superficie)",
  human_density_km2  = "pop/area km²; NOT same as densidad column (r=0.69)",
  hogares            = "Households; ≈ n_hogares (r=0.9998)",
  owned_dogs_total   = "Estimated owned dogs (UC method)",
  owned_dog_density_km2 = "Owned-dog density dogs/km²",
  human_dog_ratio    = "Persons per owned dog",
  dogs_per_urban_dwell = "Owned urban dogs / urban dwellings (24/25 valid)",
  dogs_per_rural_dwell = "Owned rural dogs / rural dwellings (22/25 valid)",
  vets               = "Vet centres (= n_vets, r=1.00)",
  n_ong              = "NGOs registered; significant univariate (p=0.012)",
  ong_cat            = "NGO presence binary",
  chip_total         = "Total microchips registered",
  enc_perro          = "National survey responses (dog owners)",
  enc_total          = "National survey responses (all species)",
  chip_per_dog       = "chip_total / owned_dogs_total",
  vets_per_1000dogs  = "Vet centres per 1,000 owned dogs",
  neutered           = "Prop. neutered — Salgado-Caxito 2023; key predictor (OR≈0.48)",
  always_inside      = "Prop. always indoors — Salgado-Caxito 2023",
  roam               = "Prop. roaming unsupervised — Salgado-Caxito 2023",
  rm_total           = "Owner-reported roaming count (18/25 zeros; use log1p)",
  tasab2019          = "Bite rate 2019 bites/100k pop (SIRAM)",
  tasa_b_sum         = "Cumulative bite rate 2019–2024",
  PC1                = "PCA axis 1: socioeconomic development (57.7% variance)",
  PC2                = "PCA axis 2: urbanisation/roaming behaviour (12.5%)",
  PC3                = "PCA axis 3: indoor containment vs neutering (10.8%)"
)

# Ordered column list
ordered_cols <- unlist(lapply(col_groups, `[[`, "cols"))
ordered_cols <- ordered_cols[ordered_cols %in% names(master)]

# Reorder master to match
master_ordered <- master %>% select(all_of(ordered_cols))

# ── F2. Build workbook ────────────────────────────────────────────────────────
wb_out <- createWorkbook()

# Styles
navy   <- createStyle(fgFill = "#1F4E79", fontColour = "#FFFFFF",
                      fontName = "Calibri", fontSize = 11, textDecoration = "bold",
                      halign = "left", valign = "center")
data_s <- createStyle(fontName = "Calibri", fontSize = 9,
                      halign = "center", valign = "center",
                      border = "BottomRight", borderColour = "#DDDDDD",
                      borderStyle = "hair")
data_l <- createStyle(fontName = "Calibri", fontSize = 9,
                      halign = "left", valign = "center",
                      border = "BottomRight", borderColour = "#DDDDDD",
                      borderStyle = "hair")
low_s  <- createStyle(fontName = "Calibri", fontSize = 9, textDecoration = "italic",
                      fgFill = "#FFF8DC", halign = "left", valign = "center",
                      border = "BottomRight", borderColour = "#DDDDDD",
                      borderStyle = "hair")
low_n  <- createStyle(fontName = "Calibri", fontSize = 9, textDecoration = "italic",
                      fgFill = "#FFF8DC", halign = "center", valign = "center",
                      border = "BottomRight", borderColour = "#DDDDDD",
                      borderStyle = "hair")

# ── SHEET 1: District Dataset ─────────────────────────────────────────────────
addWorksheet(wb_out, "District Dataset")
ws1_name <- "District Dataset"
ncol_s1  <- length(ordered_cols)
nrow_s1  <- nrow(master_ordered)

# Row 1: title banner
writeData(wb_out, ws1_name,
          "Free-roaming dog density — Chile 2021  |  Mardones et al.  |  n=25 districts  |  sorted by dogs/km descending",
          startRow = 1, startCol = 1, colNames = FALSE)
addStyle(wb_out, ws1_name, navy, rows = 1, cols = 1:ncol_s1, gridExpand = TRUE)
mergeCells(wb_out, ws1_name, rows = 1, cols = 1:ncol_s1)
setRowHeights(wb_out, ws1_name, rows = 1, heights = 20)

# Row 2: group headers
col_ptr <- 1
for (grp in col_groups) {
  grp_cols <- grp$cols[grp$cols %in% ordered_cols]
  if (length(grp_cols) == 0) next
  span <- length(grp_cols)
  grp_style <- createStyle(
    fgFill = grp$fill, fontName = "Calibri", fontSize = 9,
    textDecoration = "bold", halign = "center", valign = "center",
    border = "TopBottomLeftRight", borderColour = "#666666", borderStyle = "medium"
  )
  writeData(wb_out, ws1_name, grp$name, startRow = 2, startCol = col_ptr, colNames = FALSE)
  if (span > 1) mergeCells(wb_out, ws1_name, rows = 2, cols = col_ptr:(col_ptr+span-1))
  addStyle(wb_out, ws1_name, grp_style, rows = 2, cols = col_ptr:(col_ptr+span-1), gridExpand = TRUE)
  col_ptr <- col_ptr + span
}
setRowHeights(wb_out, ws1_name, rows = 2, heights = 16)

# Row 3: R variable names (italic, small)
key_row <- as.list(ordered_cols)
names(key_row) <- NULL
for (ci in seq_along(ordered_cols)) {
  key_style <- createStyle(
    fgFill = col_groups[[which(sapply(col_groups, function(g) ordered_cols[ci] %in% g$cols))]]$fill,
    fontName = "Calibri", fontSize = 8, textDecoration = "italic",
    fontColour = "#777777", halign = "center", valign = "center",
    border = "Bottom", borderColour = "#AAAAAA", borderStyle = "thin"
  )
  writeData(wb_out, ws1_name, ordered_cols[ci], startRow = 3, startCol = ci, colNames = FALSE)
  addStyle(wb_out, ws1_name, key_style, rows = 3, cols = ci)
}
setRowHeights(wb_out, ws1_name, rows = 3, heights = 12)

# Row 4: display labels
for (ci in seq_along(ordered_cols)) {
  key <- ordered_cols[ci]
  grp_idx <- which(sapply(col_groups, function(g) key %in% g$cols))
  hdr_style <- createStyle(
    fgFill = col_groups[[grp_idx]]$fill,
    fontName = "Calibri", fontSize = 9, textDecoration = "bold",
    halign = "center", valign = "center", wrapText = TRUE,
    border = "Bottom", borderColour = "#333333", borderStyle = "medium"
  )
  lbl <- if (key %in% names(col_labels)) col_labels[key] else key
  writeData(wb_out, ws1_name, lbl, startRow = 4, startCol = ci, colNames = FALSE)
  addStyle(wb_out, ws1_name, hdr_style, rows = 4, cols = ci)
}
setRowHeights(wb_out, ws1_name, rows = 4, heights = 36)

# Data (rows 5 onward)
text_keys <- c("comuna","macrozona","team","effort_flag")
low_rows  <- which(master_ordered$effort_flag == "low-effort")

writeData(wb_out, ws1_name, master_ordered, startRow = 5, startCol = 1, colNames = FALSE)

# Apply alternating row styles
for (ri in seq_len(nrow_s1)) {
  excel_row <- ri + 4
  is_low    <- ri %in% low_rows
  for (ci in seq_along(ordered_cols)) {
    key <- ordered_cols[ci]
    is_text <- key %in% text_keys
    sty <- if (is_low) (if(is_text) low_s else low_n) else {
      base <- if (ri %% 2 == 0) "#FAFAFA" else "#FFFFFF"
      if (is_text)
        createStyle(fontName="Calibri", fontSize=9, fgFill=base, halign="left",
                    border="BottomRight", borderColour="#DDDDDD", borderStyle="hair")
      else
        createStyle(fontName="Calibri", fontSize=9, fgFill=base, halign="center",
                    border="BottomRight", borderColour="#DDDDDD", borderStyle="hair")
    }
    addStyle(wb_out, ws1_name, sty, rows = excel_row, cols = ci)
  }
  setRowHeights(wb_out, ws1_name, rows = excel_row, heights = 13)
}

# Number formats
num_fmt_map <- c(
  dogs_per_km="0.00", dogs_per_km_lcl="0.00", dogs_per_km_ucl="0.00",
  encounter_rate="0.000", walked_km="0.00", raw_km="0.00",
  human_density_km2="0.0", owned_dog_density_km2="0.00",
  human_dog_ratio="0.000", pct_urban_pop="0.0",
  dogs_per_urban_dwell="0.000", dogs_per_rural_dwell="0.000",
  idc2="0.000", bienestar="0.000", economia="0.000",
  educacion="0.000", idh2005="0.000",
  neutered="0.000", always_inside="0.000", roam="0.000",
  chip_per_dog="0.000", vets_per_1000dogs="0.00",
  tasab2019="0.0", tasa_b_sum="0.0",
  PC1="0.000", PC2="0.000", PC3="0.000",
  pop_total="#,##0", hogares="#,##0", owned_dogs_total="#,##0",
  chip_total="#,##0", enc_perro="#,##0", enc_total="#,##0", rm_total="#,##0"
)
for (key in names(num_fmt_map)) {
  ci <- which(ordered_cols == key)
  if (length(ci) == 0) next
  fmt_style <- createStyle(numFmt = num_fmt_map[key])
  addStyle(wb_out, ws1_name, fmt_style, rows = 5:(nrow_s1+4), cols = ci,
           gridExpand = TRUE, stack = TRUE)
}

# Column widths
col_w_map <- c(
  comuna=22, macrozona=11, team=13, effort_flag=9,
  n_survey_days=8, n_records=9, n_dog_records=10, total_dogs=9,
  walked_km=8, raw_km=8, n_jumps=7,
  dogs_per_km=9, dogs_per_km_lcl=10, dogs_per_km_ucl=10,
  encounter_rate=11
)
for (ci in seq_along(ordered_cols)) {
  key <- ordered_cols[ci]
  w   <- if (key %in% names(col_w_map)) col_w_map[key] else 10
  setColWidths(wb_out, ws1_name, cols = ci, widths = w)
}

# Conditional colour scale on dogs_per_km
dpk_ci <- which(ordered_cols == "dogs_per_km")
if (length(dpk_ci) > 0) {
  conditionalFormatting(wb_out, ws1_name,
    cols = dpk_ci, rows = 5:(nrow_s1+4),
    type = "colorScale",
    style = c("#41AB5D","#FEB24C","#D7301F"),
    rule  = c(0, 0.5, 1)
  )
}

# Freeze panes at row 5 after district name (col 2)
freezePane(wb_out, ws1_name, firstActiveRow = 5, firstActiveCol = 2)

# ── SHEET 2: Variable Dictionary ──────────────────────────────────────────────
addWorksheet(wb_out, "Variable Dictionary")
ws2_name <- "Variable Dictionary"

dict_df <- data.frame(
  Variable  = ordered_cols,
  Label     = sapply(ordered_cols, function(k) if(k %in% names(col_labels)) col_labels[k] else k),
  Group     = sapply(ordered_cols, function(k)
    col_groups[[which(sapply(col_groups, function(g) k %in% g$cols))]]$name),
  n_valid   = sapply(ordered_cols, function(k)
    paste0(sum(!is.na(master_ordered[[k]])), "/25")),
  Notes     = sapply(ordered_cols, function(k)
    if(k %in% names(col_notes)) col_notes[k] else ""),
  stringsAsFactors = FALSE
)
row.names(dict_df) <- NULL

hdr_style2 <- createStyle(
  fgFill = "#1F4E79", fontColour = "#FFFFFF", fontName = "Calibri",
  fontSize = 10, textDecoration = "bold", halign = "center", valign = "center",
  border = "Bottom", borderColour = "#000000", borderStyle = "medium"
)
writeData(wb_out, ws2_name, dict_df, startRow = 1, startCol = 1, colNames = TRUE)
addStyle(wb_out, ws2_name, hdr_style2, rows = 1, cols = 1:5, gridExpand = TRUE)
setRowHeights(wb_out, ws2_name, rows = 1, heights = 16)

for (ri in seq_len(nrow(dict_df))) {
  grp_idx <- which(sapply(col_groups, function(g) dict_df$Variable[ri] %in% g$cols))
  fhex    <- if (length(grp_idx) > 0) col_groups[[grp_idx]]$fill else "#FFFFFF"
  alt     <- if (ri %% 2 == 0) "#FAFAFA" else "#FFFFFF"
  for (ci in 1:5) {
    s <- createStyle(
      fgFill = if(ci==3) fhex else alt,
      fontName = "Calibri", fontSize = 9,
      halign = if(ci %in% c(4)) "center" else "left",
      valign = "center", wrapText = (ci == 5),
      border = "Bottom", borderColour = "#DDDDDD", borderStyle = "hair"
    )
    addStyle(wb_out, ws2_name, s, rows = ri+1, cols = ci)
  }
  setRowHeights(wb_out, ws2_name, rows = ri+1, heights = 14)
}

setColWidths(wb_out, ws2_name, cols = 1:5, widths = c(24,32,30,8,55))
freezePane(wb_out, ws2_name, firstActiveRow = 2, firstActiveCol = 1)

# ── SHEET 3: Summary Statistics ───────────────────────────────────────────────
addWorksheet(wb_out, "Summary Statistics")
ws3_name <- "Summary Statistics"

stat_keys <- c("dogs_per_km","encounter_rate","idc2","economia","educacion",
               "idh2005","human_density_km2","pct_urban_pop","hogares","vets",
               "n_ong","chip_total","owned_dog_density_km2","human_dog_ratio",
               "neutered","always_inside","roam","rm_total","tasab2019")

stat_df <- do.call(rbind, lapply(stat_keys, function(k) {
  x <- as.numeric(master_ordered[[k]])
  x <- x[!is.na(x)]
  grp_idx <- which(sapply(col_groups, function(g) k %in% g$cols))
  data.frame(
    Variable = k,
    Label    = if(k %in% names(col_labels)) col_labels[k] else k,
    Group    = if(length(grp_idx)>0) col_groups[[grp_idx]]$name else "",
    n        = length(x),
    Mean     = round(mean(x), 4),
    SD       = round(sd(x),   4),
    Min      = round(min(x),  4),
    P25      = round(quantile(x, 0.25), 4),
    Median   = round(median(x), 4),
    P75      = round(quantile(x, 0.75), 4),
    Max      = round(max(x),  4),
    stringsAsFactors = FALSE
  )
}))

hdr_style3 <- createStyle(
  fgFill = "#1F4E79", fontColour = "#FFFFFF", fontName = "Calibri",
  fontSize = 10, textDecoration = "bold", halign = "center", valign = "center",
  border = "Bottom", borderColour = "#000000", borderStyle = "medium"
)
writeData(wb_out, ws3_name, stat_df, startRow = 1, colNames = TRUE)
addStyle(wb_out, ws3_name, hdr_style3, rows = 1, cols = 1:11, gridExpand = TRUE)
setRowHeights(wb_out, ws3_name, rows = 1, heights = 16)

for (ri in seq_len(nrow(stat_df))) {
  grp_idx <- which(sapply(col_groups, function(g) stat_df$Variable[ri] %in% g$cols))
  fhex    <- if(length(grp_idx)>0) col_groups[[grp_idx]]$fill else "#FFFFFF"
  alt     <- if(ri%%2==0) "#FAFAFA" else "#FFFFFF"
  for (ci in 1:11) {
    s <- createStyle(
      fgFill = if(ci==3) fhex else alt,
      fontName = "Calibri", fontSize = 9,
      halign = if(ci<=3) "left" else "center",
      valign = "center",
      border = "Bottom", borderColour = "#DDDDDD", borderStyle = "hair"
    )
    addStyle(wb_out, ws3_name, s, rows = ri+1, cols = ci)
  }
  num_s3 <- createStyle(numFmt = "0.000")
  addStyle(wb_out, ws3_name, num_s3, rows = ri+1, cols = 5:11, stack = TRUE)
  setRowHeights(wb_out, ws3_name, rows = ri+1, heights = 13)
}

setColWidths(wb_out, ws3_name, cols = 1:11,
             widths = c(22,32,26,5,10,10,10,10,10,10,10))
freezePane(wb_out, ws3_name, firstActiveRow = 2)

# ── F3. Save ─────────────────────────────────────────────────────────────────
saveWorkbook(wb_out, "data/Roaming_dogs_complete_dataset.xlsx", overwrite = TRUE)
cat("\u2713 Dataset exported -> data/Roaming_dogs_complete_dataset.xlsx\n")
cat(sprintf("  Rows: %d districts | Columns: %d variables | Sheets: 3\n",
            nrow(master_ordered), length(ordered_cols)))

# ============================================================================
# SUMMARY
# ============================================================================
cat("\n", strrep("=", 70), "\n")
cat("SUMMARY\n")
cat(strrep("=", 70), "\n")
cat(sprintf(
  "\nA-D: %d districts | %.1f dogs/km (95%% CI %.1f\u2013%.1f) | ICC=%.2f\n",
  nrow(transect), overall$dogs_per_km, overall$lcl, overall$ucl,
  performance::icc(m_nb)$ICC_adjusted))
cat(sprintf("E:   %d PCs (Kaiser, %.1f%% var) | %d districts in PCA\n",
            n_pcs, sum(eig$variance.percent[seq_len(n_pcs)]), nrow(pca_dat)))
cat("\nTop univariate associations (p < 0.15):\n")
print(TableE1_all %>% filter(p < 0.15) %>%
        select(predictor, RR_fmt, p, sig, p_fdr, fdr_sig), n = Inf)

cat("\nMultivariate NB (IDC2 + Neutered + Always indoors):\n")
if (!is.null(TableE5)) print(TableE5, n = Inf)

cat("\nOutputs written to: data/  tables/  figures/\n")
# writeLines(capture.output(sessionInfo()), "tables/sessionInfo.txt")
