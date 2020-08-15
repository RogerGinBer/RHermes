app <- ShinyDriver$new("../../", loadTimeout = 1e+05, shinyOptions = list(port = 1234),
                       phantomTimeout = 1e+05)
app$snapshotInit("mytest")

app$setInputs(shutdown = "click")
app$snapshot()
