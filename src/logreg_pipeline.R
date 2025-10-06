suppressPackageStartupMessages({
  library(readr); library(dplyr); library(rsample); library(yardstick)
  library(tibble); library(yaml)
})

# ---- config ----
cfg <- yaml::read_yaml("config/config.yml")
`%||%` <- function(a,b) if(!is.null(a)) a else b
set.seed(cfg$seed %||% 42)

# ---- read data ----
train_raw <- readr::read_tsv(cfg$train_path, show_col_types = FALSE) %>%
  mutate(label = factor(label, levels = c(0,1)))
test_raw  <- readr::read_tsv(cfg$test_path,  show_col_types = FALSE)

# ---- choose numeric features present in BOTH train & test ----
id_cols <- c("id","ref_id","CNAG_ID","CNAG_id","CNAG","gene_id","transcript_id","evalue","length")
drop_feature_cols <- cfg$features_drop %||% character(0)

numeric_train <- train_raw %>%
  select(-any_of(id_cols), -any_of(drop_feature_cols), -label) %>%
  select(where(is.numeric))

numeric_feature_cols <- intersect(names(numeric_train), names(test_raw))
stopifnot(length(numeric_feature_cols) > 0)

# ---- helpers ----
compute_stats <- function(df) {
  mu  <- vapply(df, function(x) mean(x, na.rm=TRUE),  numeric(1))
  sdv <- vapply(df, function(x) stats::sd(x, na.rm=TRUE), numeric(1))
  sdv[!is.finite(sdv) | sdv == 0] <- 1
  list(mean=mu, sd=sdv)
}
apply_pp <- function(df, st, cols){
  df %>% mutate(across(all_of(cols),
    ~ { col <- cur_column(); (.-st$mean[[col]])/st$sd[[col]] }))
}
mcc_at <- function(truth, prob, thr){
  pred <- factor(ifelse(prob >= thr, 1, 0), levels=c(0,1))
  yardstick::mcc_vec(truth=truth, estimate=pred, event_level="second")
}

# ---- 10× OOF to tune threshold ----
reps <- cfg$cv_repeats %||% 10
model_df <- train_raw %>% select(label, all_of(numeric_feature_cols))
oof <- vector("list", reps)
for(i in seq_len(reps)){
  s  <- rsample::initial_split(model_df, prop=.7, strata=label)
  tr <- rsample::training(s); te <- rsample::testing(s)
  st <- compute_stats(tr %>% select(all_of(numeric_feature_cols)))
  tr_pp <- apply_pp(tr, st, numeric_feature_cols)
  te_pp <- apply_pp(te, st, numeric_feature_cols)
  fit <- glm(label ~ ., data = tr_pp %>% select(label, all_of(numeric_feature_cols)), family = binomial())
  oof[[i]] <- tibble(
    prob  = as.numeric(predict(fit, newdata = te_pp %>% select(all_of(numeric_feature_cols)), type="response")),
    truth = te_pp$label
  )
}
oof <- dplyr::bind_rows(oof)
ths <- seq(0,1,by=.01)
mcc <- vapply(ths, function(t) mcc_at(oof$truth, oof$prob, t), numeric(1))
best_thr <- ths[which.max(mcc)]

# ---- final model & score test ----
full_stats <- compute_stats(model_df %>% select(all_of(numeric_feature_cols)))
train_pp   <- apply_pp(model_df, full_stats, numeric_feature_cols)
final_fit  <- glm(label ~ ., data=train_pp %>% select(label, all_of(numeric_feature_cols)), family=binomial())

test_pp    <- apply_pp(test_raw, full_stats, numeric_feature_cols)
test_prob  <- as.numeric(predict(final_fit, newdata=test_pp %>% select(all_of(numeric_feature_cols)), type="response"))

dir.create("outputs", showWarnings = FALSE)
readr::write_tsv(
  test_raw %>% mutate(prob_good = test_prob,
                      pred_class = as.integer(prob_good >= (cfg$threshold %||% best_thr))) %>%
    arrange(desc(prob_good)),
  "outputs/test_predictions.tsv"
)
readr::write_csv(tibble(threshold=ths, mcc=mcc, best=(ths==best_thr)), "outputs/mcc_vs_threshold.csv")
writeLines(sprintf("Best OOF MCC=%.3f at threshold=%.2f", max(mcc), best_thr), "outputs/summary.txt")
message("Done. See outputs/")
