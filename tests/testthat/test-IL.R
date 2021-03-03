context("Inclusion List generation")

test_that("Regular IL generation works",{
    myHermes <- readRDS(system.file("extdata",
                                    "exampleObject.rds",
                                    package = "RHermes"))

    myHermes <- generateIL(myHermes, 1, ILParam())
    expect_equal(nrow(myHermes@data@MS2Exp[[2]]@IL@IL),
                 nrow(myHermes@data@MS2Exp[[1]]@IL@IL))

})

test_that("Prioritized IL generation works",{
    myHermes <- readRDS(system.file("extdata",
                                    "exampleObject.rds",
                                    package = "RHermes"))
    myHermes <- generateIL(myHermes, 1, ILParam(filtermz = 0.1,
                                        priorization = "yes", ad = "M+H"))
    expect_equal(nrow(myHermes@data@MS2Exp[[2]]@IL@IL), 5)

})

test_that("IL can be exported", {
    #Reason: Generates local files on the computer (that are removed afterwards,
    #but still)
    skip_on_bioc()
    myHermes <- readRDS(system.file("extdata",
                                    "exampleObject.rds",
                                    package = "RHermes"))
    myHermes <- generateIL(myHermes, 1, ILParam())
    exportIL(myHermes, id = 1, folder = getwd(), maxOver = 5, sepFiles = FALSE)
    exportIL(myHermes, id = 1, folder = getwd(), maxOver = 5, sepFiles = TRUE)
    file.remove(c("./ExportedIL.csv",
                paste0(paste("./Injection", seq(1), sep = "_"), ".csv")))
    succeed()
})


