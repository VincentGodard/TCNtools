test_that("without filters, returns the whole internal table", {
  tab <- tcn_parameters()
  expect_s3_class(tab, "data.frame")
  expect_equal(nrow(tab), nrow(nuclide_parameters))
  expect_named(tab, c("nuclide", "target", "parameter", "pathway", "element",
                      "scaling", "value", "sd", "unit", "reference",
                      "checked", "note"))
})

test_that("filters keep the selected rows and the shared ones", {
  tab <- tcn_parameters("Be10", scaling = "St")
  expect_true(all(tab$nuclide %in% c("Be10", NA)))
  expect_true(all(tab$scaling %in% c("St", NA)))
  # one half-life, one spallation rate, two muon rates, three attenuation lengths
  expect_equal(nrow(tab), 7)
  expect_equal(sum(is.na(tab$nuclide)), 3)
})

test_that("several nuclides and schemes can be selected", {
  tab <- tcn_parameters(c("Be10", "Al26"), scaling = c("St", "Sa"))
  sp <- tab[tab$parameter == "production" & tab$pathway %in% "spallation", ]
  expect_equal(nrow(sp), 4)
  expect_setequal(sp$nuclide, c("Be10", "Al26"))
  expect_setequal(sp$scaling, c("St", "Sa"))
})

test_that("values agree with tcn_params()", {
  for (n in c("Be10", "Al26", "C14", "Ne21")) {
    for (s in c("St", "Lm", "Sf", "Sa")) {
      tab <- tcn_parameters(n, scaling = s)
      p <- tcn_params(n, scaling = s)
      sp <- tab$value[tab$parameter == "production" & tab$pathway %in% "spallation"]
      expect_equal(sp, p$production$value[p$production$pathway == "spallation"],
                   label = paste(n, s))
    }
  }
})

test_that("row names are reset", {
  tab <- tcn_parameters("C14")
  expect_identical(rownames(tab), as.character(seq_len(nrow(tab))))
})

test_that("invalid input is rejected", {
  expect_error(tcn_parameters("Xx99"), "nuclides")
  expect_error(tcn_parameters(scaling = "LSD"), "scaling")
  expect_error(tcn_parameters(target = "olivine"), "target")
  expect_error(tcn_parameters(character(0)), "nuclides")
  expect_error(tcn_parameters(c("Be10", NA)), "nuclides")
  expect_error(tcn_parameters(1), "nuclides")
})
