# u09：TyG 与 ApoB 的不一致分析 —— 按中位数分为 4 组（参照：两者均低），以及两者同时入模（互相独立的效应），M1 协变量。
suppressMessages({library(data.table); library(survival)})
base <- fread("../cache/u_base.csv"); cen <- fread("../cache/phe_censor.csv"); base <- merge(base, cen, by = "eid")
for (v in c("smk", "eth", "alc", "centre")) base[[v]] <- factor(base[[v]])
P <- do.call(cbind, lapply(sort(Sys.glob("../cache/phe_days_part*.csv")), function(f) { x <- fread(f); x <- x[match(base$eid, x$eid)]; x[, eid := NULL]; x }))
OUTS <- list(MI = P$fo_I21, IHD = P$fo_I25, Stroke = P$fo_I63, HF = P$fo_I50, T2D = P$fo_E11, MASLD = P$fo_K76, Gout = P$fo_M10, Gallstones = P$fo_K80, CKD = P$fo_N18,
             Death = ifelse(is.na(base$death_day), NA, base$death_day))
base[, grp := factor(paste0(ifelse(z_TyG > median(z_TyG), "TyG high", "TyG low"), " / ", ifelse(z_ref_ApoB > median(z_ref_ApoB), "ApoB high", "ApoB low")),
                     levels = c("TyG low / ApoB low", "TyG high / ApoB low", "TyG low / ApoB high", "TyG high / ApoB high"))]
COV <- "age + sex + tdi + smk + eth + alc + strata(centre)"; res <- list()
for (o in names(OUTS)) {
  D <- copy(base); D[, day := OUTS[[o]]]; D <- D[is.na(day) | day > 0]; D[, e := as.integer(!is.na(day) & day <= cday)]; D[, t := ifelse(e == 1, day, cday) / 365.25]; D <- D[t > 0]
  f1 <- coxph(as.formula(paste("Surv(t, e) ~ grp +", COV)), data = D, ties = "breslow"); s <- summary(f1)$conf.int[grep("^grp", rownames(summary(f1)$conf.int)), ]
  ev <- D[, .(events = sum(e)), by = grp]
  res[[o]] <- data.table(outcome = o, term = c(levels(D$grp)[-1], "TyG (mutually adjusted, per SD)", "ApoB (mutually adjusted, per SD)"),
                         HR = c(s[, 1], NA, NA), lo = c(s[, 3], NA, NA), hi = c(s[, 4], NA, NA))
  f2 <- coxph(as.formula(paste("Surv(t, e) ~ z_TyG + z_ref_ApoB +", COV)), data = D, ties = "breslow"); s2 <- summary(f2)$conf.int[c("z_TyG", "z_ref_ApoB"), ]
  res[[o]][4:5, `:=`(HR = s2[, 1], lo = s2[, 3], hi = s2[, 4])]
  res[[o]][, events := sum(D$e)]
}
R <- rbindlist(res); fwrite(R, "../results/u09_apob_discordance.csv"); print(R[, .(outcome, term, HR = round(HR, 2), lo = round(lo, 2), hi = round(hi, 2))])
