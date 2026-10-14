# u06：临床增量预测 —— 关键结局，在人口学基线（A）与临床基线（B：+收缩压、降压药、总胆固醇、HDL、降脂药、BMI）上分别加入各指数，
# 计算 Harrell C 及其相对基线的增量（95% CI，来自 survival::concordance 的协方差矩阵）。
suppressMessages({library(data.table); library(survival); library(parallel)})
base <- fread("../cache/u_base.csv"); cen <- fread("../cache/phe_censor.csv"); base <- merge(base, cen, by = "eid")
for (v in c("smk", "eth", "alc", "centre")) base[[v]] <- factor(base[[v]])
P <- do.call(cbind, lapply(sort(Sys.glob("../cache/phe_days_part*.csv")), function(f) { x <- fread(f); x <- x[match(base$eid, x$eid)]; x[, eid := NULL]; x }))
mace_day <- pmin(ifelse(P$fo_I21 > 0, P$fo_I21, NA), ifelse(P$fo_I63 > 0, P$fo_I63, NA), ifelse(base$cvdeath == 1, base$death_day, NA), na.rm = TRUE)
mace_prev <- (P$fo_I21 <= 0 & !is.na(P$fo_I21)) | (P$fo_I63 <= 0 & !is.na(P$fo_I63))
OUTS <- list(T2D = P$fo_E11, MI = P$fo_I21, Stroke = P$fo_I63, HF = P$fo_I50, AF = P$fo_I48, CKD = P$fo_N18, MASLD = P$fo_K76, Gout = P$fo_M10,
             Gallstones = P$fo_K80, MACE = ifelse(mace_prev, 0L, mace_day), Death = ifelse(is.na(base$death_day), NA, base$death_day))
IDX <- c("TyG", "TyG_fs", "TyG_BMI", "TyG_WC", "TyG_WHtR", "TyG_ABSI", "TG_HDL", "AIP", "METS_IR", "neg_eGDR", "neg_eGDR_noHTN", "CHG", "CTI", "c_HbA1c", "c_lnTG", "ref_ApoB")
A <- "age + sex + eth + tdi + smk"; Bc <- paste(A, "+ sbp + bpmed + tc + hdl + lipmed + bmi")
run <- function(o) {
  D <- copy(base); D[, day := OUTS[[o]]]; D <- D[is.na(day) | day > 0]
  D[, e := as.integer(!is.na(day) & day <= cday)]; D[, t := ifelse(e == 1, day, cday) / 365.25]; D <- D[t > 0 & !is.na(sbp) & !is.na(tc) & !is.na(hdl)]
  out <- list()
  for (bn in c("A", "B")) {
    bf <- if (bn == "A") A else Bc
    f0 <- coxph(as.formula(paste("Surv(t, e) ~", bf)), data = D, ties = "breslow")
    for (ix in IDX) {
      f1 <- coxph(as.formula(paste("Surv(t, e) ~", bf, "+", paste0("z_", ix))), data = D, ties = "breslow")
      cc <- concordance(f0, f1); v <- cc$var; d <- cc$concordance[2] - cc$concordance[1]; sd <- sqrt(v[1, 1] + v[2, 2] - 2 * v[1, 2])
      hr <- exp(coef(f1)[paste0("z_", ix)])
      out[[length(out) + 1]] <- data.table(outcome = o, base = bn, index = ix, n = f1$n, events = f1$nevent, C_base = cc$concordance[1], C_new = cc$concordance[2],
                                           dC = d, dC_lo = d - 1.96 * sd, dC_hi = d + 1.96 * sd, HR_adj = hr)
    }
  }
  rbindlist(out)
}
R <- rbindlist(mclapply(names(OUTS), function(o) tryCatch(run(o), error = function(e) {message(o, ": ", conditionMessage(e)); NULL}), mc.cores = 6))
fwrite(R, "../results/u06_prediction.csv")
print(dcast(R[base == "B"], index ~ outcome, value.var = "dC", fun.aggregate = function(x) round(x * 1000, 1)))
