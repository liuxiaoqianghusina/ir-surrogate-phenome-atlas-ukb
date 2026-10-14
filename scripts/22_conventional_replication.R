# u15：复现以往文献的常规做法 —— TyG 四分位（Q4 vs Q1）与心衰、房颤、卒中、全因死亡，以及两个阴性对照（龋齿、急性扁桃体炎）和阳性对照（T2D、心梗）。
# 模型 A：年龄、性别、种族、TDI、吸烟、饮酒（常见调整，不含 BMI）；模型 B：A + 收缩压、降压药、降脂药、LDL-C（“充分调整”但仍不含 BMI）；模型 C：B + BMI。
suppressMessages({library(data.table); library(survival)})
base <- fread("../cache/u_base.csv"); cen <- fread("../cache/phe_censor.csv"); base <- merge(base, cen, by = "eid")
for (v in c("smk", "eth", "alc", "centre")) base[[v]] <- factor(base[[v]])
base <- copy(base); base[, tyg := log(tg * 88.57 * glu * 18 / 2)]; base[, q := cut(tyg, quantile(tyg, 0:4 / 4), include.lowest = TRUE, labels = paste0("Q", 1:4))]
P <- do.call(cbind, lapply(sort(Sys.glob("../cache/phe_days_part*.csv")), function(f) { x <- fread(f); x <- x[match(base$eid, x$eid)]; x[, eid := NULL]; x }))
OUTS <- list(HF = P$fo_I50, AF = P$fo_I48, Stroke = P$fo_I63, Death = ifelse(is.na(base$death_day), NA, base$death_day), T2D = P$fo_E11, MI = P$fo_I21,
             `Dental caries (NC)` = P$fo_K02, `Acute tonsillitis (NC)` = P$fo_J03)
A <- "age + sex + eth + tdi + smk + alc + strata(centre)"; B <- paste(A, "+ sbp + bpmed + lipmed + ldl"); Cm <- paste(B, "+ bmi")
res <- list()
for (o in names(OUTS)) {
  D <- copy(base); D[, day := OUTS[[o]]]; D <- D[is.na(day) | day > 0]; D[, e := as.integer(!is.na(day) & day <= cday)]; D[, t := ifelse(e == 1, day, cday) / 365.25]
  D <- D[t > 0 & !is.na(sbp) & !is.na(ldl)]
  for (m in c("A", "B", "C")) {
    f <- coxph(as.formula(paste("Surv(t, e) ~ q +", get(if (m == "C") "Cm" else m))), data = D, ties = "breslow"); s <- summary(f)$conf.int["qQ4", ]
    res[[length(res) + 1]] <- data.table(outcome = o, model = m, n = f$n, events = f$nevent, HR_Q4 = s[1], lo = s[3], hi = s[4], p = summary(f)$coefficients["qQ4", 5])
  }
}
R <- rbindlist(res); fwrite(R, "../results/u15_conventional.csv")
R[, s := sprintf("%.2f (%.2f-%.2f) p=%.1g", HR_Q4, lo, hi, p)]; print(dcast(R, outcome ~ model, value.var = "s"))
