context("General methods")
test_that("RHermesExp methods work", {
    myHermes <- readRDS(system.file("extdata",
                                    "exampleObject.rds",
                                    package = "RHermes"))
    show(myHermes)
    expect_message(readTime(myHermes), regexp = ".")

    myHermes <- setExpParam(myHermes, template = "orbi-pos")
    expect_warning(setExpParam(myHermes, template = c("orbi-pos", "orbi-neg")))
    expect_error(setExpParam(myHermes, template = c("invalidTemplate")))

    myHermes <- setExpParam(myHermes, params = ExpParam(ppm = 10, ion = "-"))
    expect_warning(remAd(myHermes, "M+H"))

    myHermes <- setDB(myHermes)
    myHermes <- addAd(myHermes, "M+Cs", deltam = 132.905, ch = 1,  mult = 1,
                      toadd = "Cs")
    expect_warning(addAd(myHermes, "M+Cs", deltam = 132.905, ch = 0,  mult = 1,
                      toadd = "Cs"))
    remAd(myHermes, "[M-H]-")
    expect_warning(remAd(myHermes, "M+2Rb"))
})

context("MS1 processing methods")

test_that("SOI methods work", {
    myHermes <- readRDS(system.file("extdata",
                                    "exampleObject.rds",
                                    package = "RHermes"))
    myHermes <- removeSOI(myHermes, 1)
    succeed()
})
