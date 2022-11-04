context("SOI are generated successfully")

test_that("A SOI param class can be created",{
  s <- RHermes:::SOIParam()
  expect_s4_class(s, "SOIParam")
  d <- getSOIpar()
  expect_s4_class(d, "SOIParam")
})

test_that("SOI generation works",{
    myHermes <- readRDS(system.file("extdata",
                                    "exampleObject.rds",
                                    package = "RHermes"))
    myHermes <- findSOI(myHermes, getSOIpar(), 1)
    expect_equal(nrow(myHermes@data@SOI[[2]]@SOIList), 11)
})

test_that("Blank substraction is configured and works", {
    library(data.table)
    set.seed(1234)
    #Does the model work as intended?
    blank <- data.table(rt = seq(0,20,0.2),
                        rtiv = rnorm(101, 10, 3),
                        formv = "foo", isov = "M0", mz = 100)
    setkeyv(blank, "formv")
    group <- dplyr::tibble(start = c(0,10), end = c(10,20),
                        peaks = list(
                            dplyr::tibble(
                                rt = seq(0,10,0.2),
                                rtiv = rnorm(51, 10, 3)
                            ),
                            dplyr::tibble(
                                rt = seq(10,20,0.2),
                                rtiv = rnorm(51, 2e5, 3)
                            )

                        ), formula = c("foo", "foo"))

    #Check that heuristics work
    expect_false(RHermes:::firstCleaning(1, group, blank)) #Blank-like
    expect_true(RHermes:::firstCleaning(2, group, blank)) #Totally different

})

test_that("SOI plot works", {
    myHermes <- RHermesExp()
    myHermes@metadata@ExpParam@adlist <- data.frame(adduct = "M+H")
    myHermes@data@SOI <- list(RHermesSOI(
        SOIList = data.table(),
        PlotDF = data.table(rt = seq(0,10,0.2),
                            rtiv = rnorm(51, 10, 3),
                            form = rep("[C6H13O6]+", 51),
                            isov = rep("M0", 51)),
        filename = "foo"
        ))

    myHermes@data@PL <- list(RHermesPL(
        peaklist = data.table(rt = seq(0,10,0.2),
            rtiv = rnorm(51, 10, 3),
            formv = rep("[C6H13O6]+", 51),
            isov = rep("M0", 51)
        ),
        raw = data.table(), filename = "foo"
        )
    )
    myHermes@metadata@ExpParam@ionF <- list(list(),
                                            data.table(f = "C6H12O6",
                                                    ion = "[C6H13O6]+",
                                                    an = "M+H"))
    p <- RHermes::plotSOI(myHermes, 1, "C6H12O6",
                            rtrange = c(0,1500), dynamicaxis = TRUE,
                            ads = "M+H")
  expect_true(is(p, "plotly"))
})


context("SOI cleanup works")
test_that("SOI are filtered correctly", {
    myHermes <- readRDS(system.file("extdata",
                                    "exampleObject.rds",
                                    package = "RHermes"))
    myHermes <- filterSOI(myHermes, 1, 20000, TRUE)
    #Performs equal to the precalculated version
    expect_equal(nrow(myHermes@data@SOI[[1]]@SOIList), 8)
})


