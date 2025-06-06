# Health Data Debugging Prompts

## Prompt 1

**Plain Request**: 

"Broken Script: ```r, library(survival), fit <- survfi, (Surv(time, status), data = patients), plot(fit, xlab = "Days", ylab = "Survival Probability")"

**Prompt Template**:

"Fix this R script performing a Kaplan-Meier survival analysis. Make sure the syntax is correct, required packages are loaded, and the plot is properly rendered."

---

## Prompt 2

**Plain Request**: 

"My SQL query returns duplicates in patient records."

**Prompt Template**:

"Diagnose and correct an SQL query returning duplicates from a patient records table when joining to visit logs. The corrected query should return one row per patient with their most recent visit date."