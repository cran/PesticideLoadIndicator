test_that("load indicator produces known result", {
  library(readxl)
  substances <- read_excel("./Table_R_substances_example.xlsx")
  products <- read_excel("./Table_R_products_example.xlsx")

  result <- compute_pesticide_load_indicator(substances, products)

  expected_col_names <- c(
    "product", "crop",
    "sum.risk.score", "reference.sum.risk.scores",
    "formula", "amount.applied",
    "standard.dosage", "HL",
    "TL", "FL",
    "L", "STI",
    "LI"
  )

  expect_equal(names(result), expected_col_names)
  expect_equal(dim(result), c(4, 13))

  expect_equal(result$HL, c(0.6429, 0.1071, 0.0857, 0.557), tolerance=1e-3)
  expect_equal(result$TL, c(0.274, 0.1424, 0.0173, 0.216), tolerance=1e-3)
  expect_equal(result$L, c(1.087, 0.278, 0.231, 0.918), tolerance=1e-3)
  expect_equal(result$STI, c(0.5, 1.667, NA, NA), tolerance=1e-3)
  expect_equal(result$LI, c(0.543, 0.4638, NA, NA), tolerance=1e-3)

})
