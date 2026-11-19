test_that("returns one tcn_params per nuclide, named and in order", {
  p <- tcn_params_list(c("Al26", "Be10"), scaling = "Sf")
  expect_type(p, "list")
  expect_named(p, c("Al26", "Be10"))
  expect_true(all(vapply(p, inherits, logical(1), what = "tcn_params")))
  expect_identical(p$Be10, tcn_params("Be10", scaling = "Sf"))
  expect_identical(p$Al26, tcn_params("Al26", scaling = "Sf"))
})

test_that("all elements share the scaling scheme and target", {
  p <- tcn_params_list(c("Be10", "Al26", "C14", "Ne21"), scaling = "Lm")
  expect_true(all(vapply(p, function(x) x$scaling, character(1)) == "Lm"))
  expect_true(all(vapply(p, function(x) x$target, character(1)) == "quartz"))
})

test_that("shared overrides apply to every nuclide", {
  p <- tcn_params_list(c("Be10", "Al26"), scaling = "St",
                       attenuation_spallation = 150)
  for (x in p) {
    expect_equal(x$attenuation$value[x$attenuation$pathway == "spallation"],
                 150)
  }
})

test_that("a single nuclide gives a list of length one", {
  p <- tcn_params_list("Be10", scaling = "St")
  expect_length(p, 1)
  expect_named(p, "Be10")
})

test_that("invalid input is rejected", {
  expect_error(tcn_params_list(character(0), scaling = "St"), "non-empty")
  expect_error(tcn_params_list(c("Be10", NA), scaling = "St"), "NA")
  expect_error(tcn_params_list(1, scaling = "St"), "character")
  expect_error(tcn_params_list(c("Be10", "Be10"), scaling = "St"),
               "duplicates")
  expect_error(tcn_params_list(c("Be10", "Xx99"), scaling = "St"), "nuclide")
  expect_error(tcn_params_list(c("Be10", "Al26"), scaling = c("St", "Sa")),
               "scaling")
})

test_that("nuclide-specific overrides are rejected", {
  expect_error(
    tcn_params_list(c("Be10", "Al26"), scaling = "St",
                    production_spallation = 4.2),
    "production_spallation"
  )
  expect_error(
    tcn_params_list(c("Be10", "Al26"), scaling = "St", half_life = 1e6),
    "shared by all nuclides"
  )
  expect_error(
    tcn_params_list(c("Be10", "Al26"), scaling = "St",
                    attenuation_spallation = 150, production_muon_fast = 0.1),
    "production_muon_fast"
  )
})

test_that("the replacement pattern for a nuclide-specific value works", {
  p <- tcn_params_list(c("Be10", "Al26"), scaling = "St")
  p$Be10 <- tcn_params("Be10", scaling = "St", production_spallation = 4.2)
  sp <- function(x) x$production$value[x$production$pathway == "spallation"]
  expect_equal(sp(p$Be10), 4.2)
  expect_equal(sp(p$Al26), 27.93)
})
