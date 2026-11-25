# dev/

Exploratory scripts that are not part of the pipeline. Run them from the repository root; they read the frozen results in `reference/` and write nothing.

| Script | Question |
|---|---|
| `review_WTY_match_effects_tangent.R` | What did the legacy `WYT_` grouping of HOR_CHP release groups capture? Compares release-group random effects between the legacy and fixed fits |
| `review_group_covariates_tangent.R` | Does release date, or another release-group covariate in `x.df`, explain the release-group pattern? |

## Legacy `WYT_` grouping vs the fixed grouping (HOR_CHP)

Compared `reference/HOR_CHP_Rxtemp_Tmsd_2026-11-24/` (legacy `WYT_`: release-group lookup that mislabels 6 of 19 groups) with `reference/HOR_CHP_Rxtemp_Tmsd_WYTfixed_2026-11-24/` (`WYT_` from `water.yr`). The two differ only in the `WYT_drought` column.

- **Legacy fit:** the release-group random effects are tiny (SD 0.07 in the best legacy candidate, against 0.27 in the fixed fit). Across all 1,024 candidates the fixed fit's random-effect SD is 1.6 to 14 times larger (median 3.4). The legacy water-year dummies absorbed almost all of the release-group-level variation.
- **The pattern in survival by release group:** high for 2015 S2, 2015 S3, 2016 S1 to S3 (plus 2011 S4 and 2013 S3); low for 2012 S1 to S3, 2014 S2, 2014 S3 and 2015 S1.
- **Why the legacy coding fit better:** it put the high groups (2016, 2015 S2 and S3, 2013) together as `drought` and the low groups (2012, 2014, 2015 S1) together as `dry`, a partition that follows the group means. The real water-year types do not: 2012 and 2016 are both `dry` but sit at opposite ends. In the fixed fit, random effects correlate 0.63 with the legacy `drought` label and 0.01 with the true `drought` label.
- **Coefficients:** the legacy fit has large positive effects relative to `dry` (`wet` +0.33 to +0.49, `drought` +0.69 to +0.77). The fixed fit has `wet` at about -0.2 to -0.5 and `drought` near zero (-0.01 to -0.21). The legacy water-year "effect" was largely an artefact of the grouping.
- **Fit:** the fixed fit is 5.6 to 11.2 log-likelihood units worse per candidate (11 to 22 in AIC), with the same number of parameters.

The unexplained release-group variation is therefore real, not explained by water-year type, and in the fixed fit it sits in the random effect.

## Do release date or other covariates explain the group pattern?

The group effect is the fixed fit's water-year effect plus release-group random effect for the best legacy candidate (`111111111110000001011`), one value per release group (n = 19). Covariates are group means of columns of `x.df`.

| Covariate | Spearman correlation with the group effect |
|---|---|
| Release date (mean day of year) | -0.31 |
| Release order within the season (S1 to S5) | 0.11 |
| Fish length (mean) | -0.11 |
| Release year (linear) | 0.63 (R-squared 0.49) |
| Release year as a factor | R-squared 0.69 |

- **Within-season timing does not explain it.** Release date and release order within the season are only weakly related to the group effect. The groups within one year differ (for example 2015 S1 is low and 2015 S2 and S3 are high), so part of the variation is within-year.
- **Most of the pattern is between years.** Year as a factor explains 69% of the variance in the group effect (5 degrees of freedom for 19 groups, so this is optimistic). The effect rises from 2011 to 2016 with dips in 2012 and 2014.
- **Scan of 331 covariates** (group means of numeric `x.df` columns, excluding detection histories, travel times and other outcome-type columns): the largest correlations are with San Joaquin gauge flow at the HOR site (`Qsjg.hor.*`, Spearman about -0.6 to -0.7; higher flow, lower group effect). A permutation test of the largest correlation across all 331 columns, which accounts for the scan, gives p = 0.063 (observed 0.70; the 95th percentile under permutation is 0.71). That is not clearly distinguishable from chance. The columns are strongly collinear.
- **Caveats:** 19 groups, one target (the best legacy candidate), no uncertainty on the random-effect estimates, and correlations of group means with a group-level quantity. This is a description of where the variation sits, not a test of a causal covariate.

## Notes

- `review_WTY_match_effects_tangent.R` keeps "WTY" as the file name was given.
- Both scripts draw plots at the end; nothing is saved.
