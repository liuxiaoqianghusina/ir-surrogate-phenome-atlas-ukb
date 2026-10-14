# u07：剂量反应（自然样条，以 0 为参照）与亚组（性别、年龄、种族、BMI、空腹 ≥8 h），关键指数 × 关键结局，M1。
suppressMessages({library(data.table); library(survival); library(splines); library(parallel)})
base <- fread("../cache/u_base.csv"); cen <- fread("../cache/phe_censor.csv"); base <- merge(base, cen, by = "eid")
for (v in c("smk", "eth", "alc", "centre")) base[[v]] <- factor(base[[v]])
P <- do.call(cbind, lapply(sort(Sys.glob("../cache/phe_days_part*.csv")), function(f) { x <- fread(f); x <- x[match(base$eid, x$eid)]; x[, eid := NULL]; x }))
OUTS <- list(T2D = P$fo_E11, MI = P$fo_I21, Stroke = P$fo_I63, HF = P$fo_I50, AF = P$fo_I48, CKD = P$fo_N18, MASLD = P$fo_K76, Gout = P$fo_M10, Gallstones = P$fo_K80,
             Death = ifelse(is.na(base$death_day), NA, base$death_day))
IDX <- c("TyG", "TyG_fs", "neg_eGDR", "neg_eGDR_noHTN", "TyG_WC", "METS_IR", "TG_HDL", "CHG")
COV <- "age + sex + tdi + smk + eth + alc + strata(centre)"
prep <- function(o) { D <- copy(base); D[, day := OUTS[[o]]]; D <- D[is.na(day) | day > 0]
  D[, e := as.integer(!is.na(day) & day <= cday)]; D[, t := ifelse(e == 1, day, cday) / 365.25]; D[t > 0] }
base[, agegrp := ifelse(age < 60, "<60", ">=60")]; base[, bmigrp := cut(bmi, c(0, 25, 30, 100), labels = c("<25", "25-30", ">=30"), right = FALSE)]
base[, ethgrp := ifelse(eth %in% c("White", "South Asian", "Black"), as.character(eth), NA)]; base[, fastgrp := ifelse(fast >= 8, ">=8h", "<8h")]
base[, sexgrp := ifelse(sex == 1, "Men", "Women")]
run <- function(o) {
  D <- prep(o); curves <- list(); sub <- list()
  for (ix in IDX) {
    z <- D[[paste0("z_", ix)]]; D[, xx := pmin(pmax(z, quantile(z, 0.005)), quantile(z, 0.995))]
    kn <- quantile(D$xx, c(0.05, 0.35, 0.65, 0.95))
    fs <- coxph(as.formula(paste("Surv(t, e) ~ ns(xx, knots = kn[2:3], Boundary.knots = kn[c(1,4)]) +", COV)), data = D, ties = "breslow")
    fl <- coxph(as.formula(paste("Surv(t, e) ~ xx +", COV)), data = D, ties = "breslow")
    g <- seq(quantile(D$xx, 0.01), quantile(D$xx, 0.99), length.out = 80)
    Bm <- ns(g, knots = kn[2:3], Boundary.knots = kn[c(1, 4)]); B0 <- ns(0, knots = kn[2:3], Boundary.knots = kn[c(1, 4)]); Bd <- sweep(Bm, 2, B0)
    be <- coef(fs)[1:3]; V <- vcov(fs)[1:3, 1:3]; lp <- as.vector(Bd %*% be); se <- sqrt(rowSums((Bd %*% V) * Bd))
    curves[[ix]] <- data.table(outcome = o, index = ix, x = g, HR = exp(lp), lo = exp(lp - 1.96 * se), hi = exp(lp + 1.96 * se),
                               p_nl = pchisq(2 * (fs$loglik[2] - fl$loglik[2]), 2, lower.tail = FALSE))
    for (gv in c("sexgrp", "agegrp", "ethgrp", "bmigrp", "fastgrp")) {
      Dg <- D[!is.na(get(gv))]; Dg[, grp := factor(get(gv))]
      f0 <- coxph(as.formula(sprintf("Surv(t, e) ~ z_%s + grp + %s", ix, COV)), data = Dg, ties = "breslow")
      f1 <- coxph(as.formula(sprintf("Surv(t, e) ~ z_%s * grp + %s", ix, COV)), data = Dg, ties = "breslow")
      p_int <- pchisq(2 * (f1$loglik[2] - f0$loglik[2]), nlevels(Dg$grp) - 1, lower.tail = FALSE)
      for (lv in levels(Dg$grp)) {
        Dl <- Dg[grp == lv]; if (sum(Dl$e) < 30) next
        cv <- if (gv == "sexgrp") sub("sex \\+ ", "", COV) else COV
        fl2 <- tryCatch(coxph(as.formula(sprintf("Surv(t, e) ~ z_%s + %s", ix, cv)), data = Dl, ties = "breslow"), error = function(e) NULL); if (is.null(fl2)) next
        s <- summary(fl2)$conf.int[paste0("z_", ix), ]
        sub[[length(sub) + 1]] <- data.table(outcome = o, index = ix, subgroup = gv, level = lv, events = fl2$nevent, HR = s[1], lo = s[3], hi = s[4], p_int = p_int)
      }
    }
  }
  list(rbindlist(curves), rbindlist(sub))
}
RR <- mclapply(names(OUTS), function(o) tryCatch(run(o), error = function(e) {message(o, ": ", conditionMessage(e)); NULL}), mc.cores = 5)
fwrite(rbindlist(lapply(RR, `[[`, 1)), "../results/u07_rcs.csv"); fwrite(rbindlist(lapply(RR, `[[`, 2)), "../results/u07_subgroups.csv")
S <- rbindlist(lapply(RR, `[[`, 2)); print(S[index %in% c("TyG", "neg_eGDR") & subgroup == "ethgrp", .(outcome, index, level, events, HR = round(HR, 2), p_int = signif(p_int, 2))])
