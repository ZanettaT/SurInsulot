# Internal helpers (not exported).

# Coerce a "male" argument to a 0/1 indicator.
# Accepts logical (TRUE/FALSE) or the numeric codes 1 (man) / 0 (woman).
# Anything else (e.g. a raw sex code of 2 for women) becomes NA, forcing the
# caller to pass an explicit indicator rather than a study-specific code.
.male_indicator <- function(male) {
  ifelse(male == 1, 1, ifelse(male == 0, 0, NA_real_))
}
