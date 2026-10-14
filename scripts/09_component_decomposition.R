# u05：成分拆解 —— 每个结局：ln TG、ln 血糖、HDL（取负）、腰围、HbA1c、ApoB 同时入模（M1 协变量），得到各成分的独立效应（每 SD）。
suppressMessages({library(data.table); library(survival); library(parallel)})
base <- fread("../cache/u_base.csv"); cen <- fread("../cache/phe_censor.csv"); base <- merge(base, cen, by = "eid")
for (v in c("smk", "eth", "alc", "centre")) base[[v]] <- factor(base[[v]])
COMP <- c("z_c_lnTG", "z_c_lnGlu", "z_c_negHDL", "z_c_WC", "z_c_HbA1c", "z_ref_ApoB")
run <- function(o, days) {
  D <- base[, c("age", "sex", "tdi", "smk", "eth", "alc", "centre", "cday", COMP), with = FALSE]; D[, day := days]
  D <- D[is.na(day) | (day > 0)]; D[, e := as.integer(!is.na(day) & day <= cday)]; D[, t := ifelse(e == 1, day, cday) / 365.25]; D <- D[t > 0]
  sx <- D[e == 1, mean(sex)]; ss <- sx > 0.9 | sx < 0.1; if (ss) D <- D[sex == round(sx)]
  if (sum(D$e) < 100) return(NULL)
  f <- as.formula(paste("Surv(t, e) ~", paste(COMP, collapse = " + "), "+ age +", if (!ss) "sex +", "tdi + smk + eth + alc + strata(centre)"))
  fit <- tryCatch(suppressWarnings(coxph(f, data = D, ties = "breslow")), error = function(e) NULL); if (is.null(fit)) return(NULL)
  s <- summary(fit)$coefficients[COMP, ]
  data.table(outcome = o, component = sub("^z_(c_|ref_)", "", COMP), events = fit$nevent, logHR = s[, 1], se = s[, 3], p = s[, 5])
}
res <- list()
for (pf in sort(Sys.glob("../cache/phe_days_part*.csv"))) {
  P <- fread(pf); P <- P[match(base$eid, P$eid)]; outs <- setdiff(names(P), "eid")
  res[[pf]] <- rbindlist(mclapply(outs, function(o) tryCatch(run(o, P[[o]]), error = function(e) NULL), mc.cores = 2))
  fwrite(rbindlist(res), "../results/u05_decompose.csv"); cat(basename(pf), "done\n")
}
D0 <- base[, c("age", "sex", "tdi", "smk", "eth", "alc", "centre", "cday", "death_day", "cvdeath", COMP), with = FALSE]
for (o in c("death", "cvdeath")) {
  D <- copy(D0); D[, e := as.integer(!is.na(death_day) & (o == "death" | cvdeath == 1))]; D[, t := cday / 365.25]; D <- D[t > 0]
  fit <- coxph(as.formula(paste("Surv(t, e) ~", paste(COMP, collapse = " + "), "+ age + sex + tdi + smk + eth + alc + strata(centre)")), data = D, ties = "breslow")
  s <- summary(fit)$coefficients[COMP, ]
  res[[o]] <- data.table(outcome = o, component = sub("^z_(c_|ref_)", "", COMP), events = fit$nevent, logHR = s[, 1], se = s[, 3], p = s[, 5])
}
fwrite(rbindlist(res), "../results/u05_decompose.csv"); cat("done\n")
