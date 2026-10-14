# u03：全表型 Cox —— 20 个指数/成分/变体 × 全部合格结局（M1）；12 个不含体型项的指数另做 M2（+BMI）。
# M1：年龄、性别、TDI、吸烟、种族、饮酒，按评估中心分层。性别特异结局（≥90% 事件为同一性别）只在该性别中分析。
suppressMessages({library(data.table); library(survival); library(parallel)})
base <- fread("../cache/u_base.csv")
cen <- fread("../cache/phe_censor.csv"); base <- merge(base, cen, by = "eid")
for (v in c("smk", "eth", "alc", "centre")) base[[v]] <- factor(base[[v]])
IDX <- c("TyG", "TyG_BMI", "TyG_WC", "TyG_WHtR", "TyG_ABSI", "TG_HDL", "AIP", "METS_IR", "neg_eGDR", "neg_eGDR_noHTN", "CHG", "CTI",
         "c_lnTG", "c_lnGlu", "c_negHDL", "c_BMI", "c_WC", "c_HbA1c", "ref_ApoB", "TyG_fs")
IDX2 <- c("TyG", "TG_HDL", "AIP", "neg_eGDR", "neg_eGDR_noHTN", "CHG", "CTI", "c_lnTG", "c_lnGlu", "c_negHDL", "c_HbA1c", "ref_ApoB", "TyG_fs")
fit1 <- function(D, ix, adj_bmi, sexspec) {
  cov <- paste(c("age", if (!sexspec) "sex", "tdi", "smk", "eth", "alc", if (adj_bmi) "bmi"), collapse = " + ")
  f <- as.formula(sprintf("Surv(t, e) ~ z_%s + %s + strata(centre)", ix, cov))
  fit <- tryCatch(suppressWarnings(coxph(f, data = D, ties = "breslow")), error = function(e) NULL)
  if (is.null(fit)) return(NULL)
  s <- summary(fit)$coefficients[paste0("z_", ix), ]
  data.table(index = ix, model = if (adj_bmi) "M2" else "M1", n = fit$n, events = fit$nevent, logHR = s[1], se = s[3], p = s[5])
}
run_outcome <- function(o, days) {
  D <- base[, c("eid", "age", "sex", "tdi", "smk", "eth", "alc", "centre", "bmi", "cday", paste0("z_", IDX)), with = FALSE]
  D[, day := days]
  D <- D[is.na(day) | (day > 0)]                      # 排除基线已患（day<=0）及日期未知（-99999）
  D[, e := as.integer(!is.na(day) & day <= cday)]; D[, t := ifelse(e == 1, day, cday) / 365.25]; D <- D[t > 0]
  sx <- D[e == 1, mean(sex)]; sexspec <- sx > 0.9 | sx < 0.1
  if (sexspec) D <- D[sex == round(sx)]
  if (sum(D$e) < 100) return(NULL)
  r1 <- rbindlist(lapply(IDX, function(ix) fit1(D, ix, FALSE, sexspec)))
  r2 <- rbindlist(lapply(IDX2, function(ix) fit1(D, ix, TRUE, sexspec)))
  r <- rbind(r1, r2); r[, `:=`(outcome = o, sex_specific = sexspec)]; r
}
parts <- sort(Sys.glob("../cache/phe_days_part*.csv")); res <- list(); t0 <- Sys.time()
for (pf in parts) {
  P <- fread(pf); P <- P[match(base$eid, P$eid)]
  outs <- setdiff(names(P), "eid")
  rr <- mclapply(outs, function(o) tryCatch(run_outcome(o, P[[o]]), error = function(e) NULL), mc.cores = 10)
  res[[pf]] <- rbindlist(rr); fwrite(rbindlist(res), "../results/u03_phewas.csv")
  cat(basename(pf), length(outs), "outcomes;", round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1), "min\n")
}
# 死亡与心血管死亡
D0 <- base[, c("eid", "age", "sex", "tdi", "smk", "eth", "alc", "centre", "bmi", "cday", "death_day", "cvdeath", paste0("z_", IDX)), with = FALSE]
for (o in c("death", "cvdeath")) {
  D <- copy(D0); D[, e := as.integer(!is.na(death_day) & (o == "death" | cvdeath == 1))]; D[, t := cday / 365.25]; D <- D[t > 0]
  r <- rbind(rbindlist(lapply(IDX, function(ix) fit1(D, ix, FALSE, FALSE))), rbindlist(lapply(IDX2, function(ix) fit1(D, ix, TRUE, FALSE))))
  r[, `:=`(outcome = o, sex_specific = FALSE)]; res[[o]] <- r
}
R <- rbindlist(res); fwrite(R, "../results/u03_phewas.csv"); cat("done", nrow(R), "rows;", length(unique(R$outcome)), "outcomes\n")
