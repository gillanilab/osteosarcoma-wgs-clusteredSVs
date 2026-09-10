# ---------------------------------------------------------------------------
# Figure 5A -- GLMM: truncal ~ hotspot + (1 | patient)
#
# This scripts require the reg_df_for_glmm.csv input.
# It is generated in the "Fig5A_final.ipynb". 
# This script reads that CSV and fits the GLMM:
#   glmer(truncal ~ hotspot + (1 | patient), family = binomial)
# ---------------------------------------------------------------------------

library(lme4)
library(ggplot2)

reg_df <- read.csv("reg_df_for_glmm.csv")

# Sanity checks to catch column/types
stopifnot(all(c("patient", "hotspot", "truncal") %in% colnames(reg_df)))
stopifnot(all(reg_df$hotspot %in% c(0, 1)))
stopifnot(all(reg_df$truncal %in% c(0, 1)))

reg_df$patient <- as.factor(reg_df$patient)

cat(sprintf("Rows (unique clustered SVs): %d\n", nrow(reg_df)))
cat(sprintf("Unique patients: %d\n", length(unique(reg_df$patient))))
cat("SVs per patient:\n")
print(summary(as.numeric(table(reg_df$patient))))


# Fit the model
model <- glmer(truncal ~ hotspot + (1 | patient), data = reg_df, family = binomial)
print(summary(model))

# Singularity check for patient variance 

# isSingular() flags whether the estimated variance of the patient-level random intercept 
# has collapsed to (or is numerically indistinguishable from) zero
singular <- isSingular(model)
cat(sprintf("\nisSingular(model): %s\n", singular))
if (singular) {
  cat("WARNING: singular fit -- random-effect variance is at/near zero, do not trust this model's p-value without addressing it.\n")
}


# Effect size (odds ratio) with both Wald and profile-likelihood CIs, plus a likelihood-ratio test

# Pulls intercept and coefficient (hotspot) from model
print(exp(fixef(model)))

# Wald for approximation of CIs and compare to ci-profile
print(exp(confint(model, method = "Wald")))

# Models CI profile (should be same or similar to Wald)
ci_profile <- exp(confint(model, method = "profile"))
print(ci_profile)

# Null model for glmer (with hotspot removed)
model_null <- glmer(truncal ~ 1 + (1 | patient), data = reg_df, family = binomial)
# Anova test to compare the output of both models
lrt <- anova(model_null, model)
print(lrt)


# Leave-one-patient-out sensitivity analysis

# Refit the model once per patient, excluding that patient each time, and track how much the 
# hotspot OR and p-value move. This checks whether the result is a broad, patient-general 
#effect or is still being driven disproportionately by one or two patients.
patients <- levels(reg_df$patient)
lopo_results <- data.frame(
  excluded_patient = character(), n_svs_excluded = integer(),
  OR = numeric(), ci_low = numeric(), ci_high = numeric(),
  p_value = numeric(), singular = logical()
)

for (pid in patients) {
  sub_df <- reg_df[reg_df$patient != pid, ]
  sub_df$patient <- droplevels(sub_df$patient)

  fit <- tryCatch(
    glmer(truncal ~ hotspot + (1 | patient), data = sub_df, family = binomial),
    warning = function(w) w, error = function(e) e
  )

  if (inherits(fit, "warning") || inherits(fit, "error")) {
    lopo_results <- rbind(lopo_results, data.frame(
      excluded_patient = pid,
      n_svs_excluded = sum(reg_df$patient == pid),
      OR = NA, ci_low = NA, ci_high = NA, p_value = NA, singular = NA
    ))
    next
  }

  coefs <- summary(fit)$coefficients
  or_val <- exp(coefs["hotspot", "Estimate"])
  se_val <- coefs["hotspot", "Std. Error"]
  p_val <- coefs["hotspot", "Pr(>|z|)"]
  ci <- exp(coefs["hotspot", "Estimate"] + c(-1.96, 1.96) * se_val)

  lopo_results <- rbind(lopo_results, data.frame(
    excluded_patient = pid,
    n_svs_excluded = sum(reg_df$patient == pid),
    OR = or_val, ci_low = ci[1], ci_high = ci[2],
    p_value = p_val, singular = isSingular(fit)
  ))
}

# Leave-one-patient-out results (sorted by SVs excluded, descending)
print(lopo_results[order(-lopo_results$n_svs_excluded), ], row.names = FALSE)

cat(sprintf(
  "\nFull-data OR = %.3f. LOPO OR range: %.3f to %.3f.\n",
  exp(fixef(model))["hotspot"], min(lopo_results$OR, na.rm = TRUE), max(lopo_results$OR, na.rm = TRUE)
))
cat(sprintf(
  "Patients where excluding them flips significance (p crosses 0.05): %d of %d\n",
  sum(lopo_results$p_value > 0.05, na.rm = TRUE), nrow(lopo_results)
))
# If the OR stays fairly stable and p stays well under 0.05 across all leave-one-out fits,
# the effect is not being driven by one patient. 

# Removing patient with most SVs to test if OR moves but does not
largest_patient <- lopo_results$excluded_patient[which.max(lopo_results$n_svs_excluded)]
cat(sprintf("\nLargest-SV patient in the cohort: %s (%d SVs). Its LOPO row:\n",
            largest_patient, max(lopo_results$n_svs_excluded)))
print(lopo_results[lopo_results$excluded_patient == largest_patient, ], row.names = FALSE)

