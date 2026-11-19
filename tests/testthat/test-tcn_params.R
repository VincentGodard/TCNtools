test_that("spallation rates are those of Borchers et al. (2016), Table 7", {
  table7 <- rbind(
    Be10 = c(St = 4.01, Lm = 4.00, Sf = 4.09, Sa = 3.92),
    Al26 = c(St = 27.93, Lm = 27.93, Sf = 28.61, Sa = 28.54),
    C14 = c(St = 12.24, Lm = 12.22, Sf = 12.72, Sa = 12.76)
  )
  for (n in rownames(table7)) {
    for (s in colnames(table7)) {
      p <- suppressWarnings(tcn_params(n, scaling = s))
      expect_equal(
        p$production$value[p$production$pathway == "spallation"],
        table7[n, s],
        label = paste(n, s)
      )
    }
  }
})

test_that("muon rates are those of Braucher et al. (2011) and Lupker et al. (2015)", {
  muons <- function(n) {
    p <- suppressWarnings(tcn_params(n, scaling = "St"))$production
    setNames(p$value, p$pathway)[c("muon_stop", "muon_fast")]
  }
  expect_equal(muons("Be10"), c(muon_stop = 0.012, muon_fast = 0.039))
  expect_equal(muons("Al26"), c(muon_stop = 0.84, muon_fast = 0.081))
  expect_equal(muons("C14"), c(muon_stop = 3.31, muon_fast = 0))
})

test_that("muon rates and attenuation lengths do not depend on the scheme", {
  st <- tcn_params("Be10", scaling = "St")
  sa <- tcn_params("Be10", scaling = "Sa")
  expect_equal(st$attenuation, sa$attenuation)
  expect_equal(st$production[-1, ], sa$production[-1, ], ignore_attr = TRUE)
  expect_equal(st$attenuation$value, c(160, 1500, 4320))
})

test_that("decay constant follows from the half-life", {
  be <- tcn_params("Be10", scaling = "St")
  expect_equal(be$nuclide$lambda[["value"]], log(2) / 1.387e6)
  expect_equal(be$nuclide$lambda[["sd"]], log(2) * 0.012e6 / 1.387e6^2)
  c14 <- tcn_params("C14", scaling = "St")
  expect_equal(c14$nuclide$lambda[["value"]], log(2) / 5700)
  al <- tcn_params("Al26", scaling = "St")
  expect_equal(al$nuclide$half_life[["value"]], 0.705e6)
  expect_equal(al$nuclide$lambda[["value"]], log(2) / 0.705e6)
})

test_that("a stable nuclide has a zero decay constant", {
  ne <- tcn_params("Ne21", scaling = "St")
  expect_identical(ne$nuclide$lambda[["value"]], 0)
  expect_identical(ne$nuclide$half_life[["value"]], Inf)
})

test_that("every supported nuclide and scheme has complete parameters", {
  for (n in c("Be10", "Al26", "C14", "Ne21")) {
    for (s in c("St", "Lm", "Sf", "Sa")) {
      expect_no_warning(tcn_params(n, scaling = s))
    }
  }
})

test_that("21Ne spallation rates are those of CRONUScalc", {
  cronuscalc <- c(St = 16.63, Lm = 16.59, Sf = 16.96, Sa = 16.26)
  for (s in names(cronuscalc)) {
    p <- tcn_params("Ne21", scaling = s)$production
    expect_equal(p$value[p$pathway == "spallation"], cronuscalc[[s]],
                 label = paste("Ne21", s))
  }
  p <- tcn_params("Ne21", scaling = "St")$production
  expect_equal(setNames(p$value, p$pathway)[c("muon_stop", "muon_fast")],
               c(muon_stop = 0, muon_fast = 0))
})

test_that("21Ne St rate agrees with independent calibrations", {
  p <- tcn_params("Ne21", scaling = "St")$production
  st <- p$value[p$pathway == "spallation"]
  # SPICE (Fenton et al. 2019): 17.0 +/- 1.1 at/g/yr, 2 sigma
  expect_lte(abs(st - 17.0), 1.1)
  # Balco v3 online calculator (make_consts_v3.m): 16.896 +/- 1.033
  expect_lte(abs(st - 16.896), 1.033)
})

test_that("overrides replace the value and record the source", {
  p <- tcn_params("Be10", scaling = "St", production_spallation = c(4.2, 0.2))
  expect_equal(p$production$value[1], 4.2)
  expect_equal(p$production$sd[1], 0.2)
  s <- as.data.frame(p)
  expect_identical(s$reference[s$parameter == "production_spallation"], "user")

  p <- tcn_params("Be10", scaling = "St", attenuation_spallation = 150)
  expect_equal(p$attenuation$value[1], 150)
  expect_true(is.na(p$attenuation$sd[1]))
})

test_that("invalid input is rejected with an informative error", {
  expect_error(tcn_params("Be10"), "scaling")
  expect_error(tcn_params("Be7", scaling = "St"), "nuclide")
  expect_error(tcn_params("Be10", scaling = "LSDn"), "scaling")
  expect_error(tcn_params(c("Be10", "Al26"), scaling = "St"), "nuclide")
  expect_error(tcn_params("Be10", scaling = "St", target = "olivine"), "target")
  expect_error(tcn_params("Be10", scaling = "St", 4.2), "named")
  expect_error(tcn_params("Be10", scaling = "St", foo = 1), "Unknown")
  expect_error(tcn_params("Be10", scaling = "St", half_life = -1), "non-negative")
  expect_error(tcn_params("Be10", scaling = "St", half_life = 1:3), "length")
})

test_that("every value in the parameter table has a source", {
  p <- nuclide_parameters
  ok <- is.na(p$value) | !is.na(p$reference) | is.infinite(p$value)
  expect_true(all(ok))
})

test_that("print and as.data.frame methods work", {
  p <- tcn_params("Be10", scaling = "Sa")
  expect_output(print(p), "Be10 in quartz, scaling scheme Sa")
  d <- as.data.frame(p)
  expect_s3_class(d, "data.frame")
  expect_true(all(c("parameter", "value", "reference") %in% names(d)))
})

test_that("10Be and 26Al spallation uncertainties are those of Phillips et al. (2016)", {
  sd <- rbind(
    Be10 = c(St = 0.33, Lm = 0.32, Sf = 0.35, Sa = 0.31),
    Al26 = c(St = 2.8, Lm = 2.7, Sf = 3.3, Sa = 3.1)
  )
  for (n in rownames(sd)) {
    for (s in colnames(sd)) {
      p <- tcn_params(n, scaling = s)$production
      expect_equal(p$sd[p$pathway == "spallation"], sd[n, s],
                   label = paste(n, s))
    }
  }
  p <- tcn_params("C14", scaling = "St")$production
  expect_true(is.na(p$sd[p$pathway == "spallation"]))
})

test_that("print shows uncertainties with three significant digits", {
  expect_output(print(tcn_params("Ne21", scaling = "Sa")), "1.96")
})
