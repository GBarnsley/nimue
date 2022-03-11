context("Deterministic model")
# library(nimue)
library(squire)

test_that("compare deterministic vaccine model to SEEIR model", {
  pop <- get_population("Angola")
  mm <- get_mixing_matrix("Angola")

  # Reference model
  m1 <- run_deterministic_SEIR_model(
    population = pop$n,
    contact_matrix_set = mm,
    hosp_bed_capacity = 100000,
    ICU_bed_capacity = 1000000,
    seed = 1,
    day_return = TRUE,
    replicates = 1,
    seeding_cases = 20,
    time_period = 100,
  )
  oi1 <- odin_index(m1$model)

  # Vaccine model, no vaccine
  m2 <- run(
    population = pop$n,
    contact_matrix_set = mm,
    hosp_bed_capacity = 100000,
    ICU_bed_capacity = 1000000,
    dur_R = Inf,
    first_doses = 0,
    seed = 1,
    replicates = 1,
    seeding_cases = 20,
    time_period = 100,
  )
  oi2 <- odin_index(m2$model)

  # Compare shared compartments (ones we expect to be equal)
  compare_compartments <- names(oi1)[names(oi1) %in% names(oi2)][2:7]

  for(i in seq_along(compare_compartments)){
    # Isolare squire output
    t1 <- m1$output[,unlist(oi1[compare_compartments[i]]),]
    # Isolate nimue output and collapse vaccine compartments if needed
    if(is.matrix(oi2[[compare_compartments[i]]])){
      t2 <- apply(oi2[[compare_compartments[i]]], 1, function(x, y){
        rowSums(y[,x,1])
      }, y = m2$output)
    } else {
      t2 <- m2$output[,unlist(oi2[compare_compartments[i]]),]
    }
    # Clear attributes
    attributes(t1) <- NULL
    attributes(t2) <- NULL

    expect_equal(t1, t2, tol = 0.00001)
  }

  # Check all vaccine-related compartments are 0
  expect_equal(sum(m2$output[,unlist(oi2[c("vaccinated_first_dose", "vaccinated_second_dose", "vaccines", "vaccinated_waned")]),]), 0)

  # Check population size is constant at specified level
  expect_equal(format(m2, "N", NULL)$value,
               rep(sum(pop$n), 100))

})


test_that("Vaccine on works", {
  pop <- get_population("Angola")
  mm <- get_mixing_matrix("Angola")

  # Vaccine model 100% efficacy against infection
  m1 <- run(
    population = pop$n,
    contact_matrix_set = mm,
    hosp_bed_capacity = 100000,
    ICU_bed_capacity = 1000000,
    dur_R = Inf,
    first_doses = 10000,
    seed = 1,
    replicates = 1,
    time_period = 100,
    seeding_cases = 20
  )

  # Check individuals reaching V
  expect_gt(sum(format(m1, "vaccinated_first_dose", NULL)$value), 0)
  expect_gt(sum(format(m1, "vaccinated_second_dose", NULL)$value), 0)
  expect_gt(sum(format(m1, "vaccinated_waned", NULL)$value), 0)

  # Check population size is constant at specified level
  expect_equal(format(m1, "N", NULL)$value,
               rep(sum(pop$n), 100))


  # Vaccine model 100% efficacy against infection (infinite vaccine duration)
  m2 <- run(
    population = pop$n,
    contact_matrix_set = mm,
    hosp_bed_capacity = 100000,
    ICU_bed_capacity = 1000000,
    dur_R = Inf,
    first_doses = 10000,
    seed = 1,
    replicates = 1,
    time_period = 100,
    seeding_cases = 20,
    dur_V = c(Inf, Inf, Inf)
  )

  # Check individuals reaching V
  expect_gt(sum(format(m2, "vaccinated_first_dose", NULL)$value), 0)
  expect_gt(sum(format(m2, "vaccinated_second_dose", NULL)$value), 0)
  # But not leaving v
  expect_equal(sum(format(m2, "vaccinated_waned", NULL)$value), 0)

  # Check population size is constant at specified level
  expect_equal(format(m2, "N", NULL)$value,
               rep(sum(pop$n), 100))

  # Vaccine model 100% efficacy against infection (no second doses)
  m_3 <- run(
    population = pop$n,
    contact_matrix_set = mm,
    hosp_bed_capacity = 100000,
    ICU_bed_capacity = 1000000,
    dur_R = Inf,
    first_doses = 10000,
    second_doses = 0,
    seed = 1,
    replicates = 1,
    time_period = 100,
    seeding_cases = 20
  )

  # Check individuals reaching V
  expect_gt(sum(format(m_3, "vaccinated_first_dose", NULL)$value), 0)
  #but no second
  expect_equal(sum(format(m_3, "vaccinated_second_dose", NULL)$value), 0)
  # and not leaving v
  expect_equal(sum(format(m_3, "vaccinated_waned", NULL)$value), 0)

  # Check population size is constant at specified level
  expect_equal(format(m_3, "N", NULL)$value,
               rep(sum(pop$n), 100))

  # Vaccine model 100% efficacy against infection (no boosters)
  m_4 <- run(
    population = pop$n,
    contact_matrix_set = mm,
    hosp_bed_capacity = 100000,
    ICU_bed_capacity = 1000000,
    dur_R = Inf,
    first_doses = 10000,
    booster_doses = 0,
    seed = 1,
    replicates = 1,
    time_period = 100,
    seeding_cases = 20
  )

  # Check individuals reaching V
  expect_gt(sum(format(m_4, "vaccinated_first_dose", NULL)$value), 0)
  expect_gt(sum(format(m_4, "vaccinated_second_dose", NULL)$value), 0)
  # and less waned that in baseline
  expect_gt(sum(format(m_4, "vaccinated_waned", NULL)$value), sum(format(m1, "vaccinated_waned", NULL)$value))

  # Check population size is constant at specified level
  expect_equal(format(m_4, "N", NULL)$value,
               rep(sum(pop$n), 100))
})

test_that("Time-varying works", {
  pop <- get_population("Angola")
  mm <- get_mixing_matrix("Angola")

  # Vaccine model time varying
  m1 <- run(
    population = pop$n,
    contact_matrix_set = mm,
    hosp_bed_capacity = 100000,
    ICU_bed_capacity = 1000000,
    dur_R = Inf,
    first_doses = c(0, 1000, 0),
    tt_first_doses = c(0, 10, 20),
    second_doses = c(0, 1000, 0),
    tt_second_doses = c(0, 20, 30),
    booster_doses = 0,
    seed = 1,
    replicates = 1,
    seeding_cases = 20,
    time_period = 100
  )

  # Check individuals in youngest age group reaching V
  t_v <- format(m1, NULL, "vaccines")
  expect_equal(sum(dplyr::filter(t_v, t < 10)$value, na.rm = TRUE), 0)
  expect_gt(sum(dplyr::filter(t_v, t >= 10, t <30)$value, na.rm = TRUE), 0)
  expect_equal(sum(dplyr::filter(t_v, t > 31)$value, na.rm = TRUE), 0)
})

test_that("Efficacy against infection works", {
  pop <- get_population("Angola")
  mm <- get_mixing_matrix("Angola")

  # No vaccine
  m1 <- run(
    population = pop$n,
    contact_matrix_set = mm,
    hosp_bed_capacity = 100000,
    ICU_bed_capacity = 1000000,
    dur_R = Inf,
    seed = 1,
    replicates = 1,
    seeding_cases = 20,
    time_period = 100
  )
  # Vaccine 50% efficacy against infection
  m2 <- run(
    population = pop$n,
    contact_matrix_set = mm,
    hosp_bed_capacity = 100000,
    ICU_bed_capacity = 1000000,
    dur_R = Inf,
    seed = 1,
    replicates = 1,
    seeding_cases = 20,
    time_period = 100,
    first_doses = 10000,
    vaccine_efficacy_infection = rep(0.5, 4),
    vaccine_efficacy_disease = rep(0, 4)
  )
  # Vaccine 100% efficacy against infection
  m3 <- run(
    population = pop$n,
    contact_matrix_set = mm,
    hosp_bed_capacity = 100000,
    ICU_bed_capacity = 1000000,
    dur_R = Inf,
    seed = 1,
    replicates = 1,
    seeding_cases = 20,
    time_period = 100,
    first_doses = 10000,
    vaccine_efficacy_infection = rep(1, 4),
    vaccine_efficacy_disease = rep(0, 4)
  )

  i1 <- sum(format(m1, NULL, "infections")$value, na.rm = TRUE)
  i2 <- sum(format(m2, NULL, "infections")$value, na.rm = TRUE)
  i3 <- sum(format(m3, NULL, "infections")$value, na.rm = TRUE)

  expect_gt(i1, i2)
  expect_gt(i2, i3)
})

test_that("Efficacy against disease works", {
  pop <- get_population("Angola")
  mm <- get_mixing_matrix("Angola")

  # No vaccine
  m1 <- run(
    population = pop$n,
    contact_matrix_set = mm,
    hosp_bed_capacity = 100000,
    ICU_bed_capacity = 1000000,
    dur_R = Inf,
    seed = 1,
    replicates = 1,
    seeding_cases = 20,
    time_period = 100
  )
  # Vaccine 50% efficacy against disease
  m2 <- run(
    population = pop$n,
    contact_matrix_set = mm,
    hosp_bed_capacity = 100000,
    ICU_bed_capacity = 1000000,
    dur_R = Inf,
    seed = 1,
    replicates = 1,
    seeding_cases = 20,
    time_period = 100,
    first_doses = 10000,
    vaccine_efficacy_disease = rep(0.5, 4),
    vaccine_efficacy_infection = rep(0, 4)
  )
  # Vaccine 100% efficacy against disease
  m3 <- run(
    population = pop$n,
    contact_matrix_set = mm,
    hosp_bed_capacity = 100000,
    ICU_bed_capacity = 1000000,
    dur_R = Inf,
    seed = 1,
    replicates = 1,
    seeding_cases = 20,
    time_period = 100,
    first_doses = 10000,
    vaccine_efficacy_disease = rep(1, 4),
    vaccine_efficacy_infection = rep(0, 4)
  )

  i1 <- sum(format(m1, NULL, "deaths")$value, na.rm = TRUE)
  i2 <- sum(format(m2, NULL, "deaths")$value, na.rm = TRUE)
  i3 <- sum(format(m3, NULL, "deaths")$value, na.rm = TRUE)

  expect_gt(i1, i2)
  expect_gt(i2, i3)
})


test_that("Efficacy parameterisation options work", {
  run1 <- run(country = "France", first_doses = 5000000, vaccine_efficacy_infection = rep(0.9, 4), seed = 1)
  run2 <- run(country = "France", first_doses = 5000000, vaccine_efficacy_infection = matrix(0.9, nrow = 4, ncol = 17), seed = 1)
  expect_identical(run1$output, run2$output)

  run3 <- run(country = "France", first_doses = 5000000, vaccine_efficacy_disease = rep(0.9, 4), seed = 1)
  run4 <- run(country = "France", first_doses = 5000000, vaccine_efficacy_disease = matrix(0.9, nrow = 4, ncol = 17), seed = 1)
  expect_identical(run3$output, run4$output)

  expect_error(run(country = "France", vaccine_efficacy_infection = rep(0.9, 2)), "If element of vaccine_efficacy_infection is a vector, it must have 4 values corresponding to first dose, second dose, and the two waning levels")
  expect_error(run(country = "France", vaccine_efficacy_disease = rep(0.9, 2)), "If element of vaccine_efficacy_disease is a vector, it must have 4 values corresponding to first dose, second dose, and the two waning levels")
})
