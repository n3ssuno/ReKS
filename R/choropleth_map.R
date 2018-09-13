# https://deanattali.com/2015/04/21/r-package-shiny-app/
choropleth_map <- function(){
    appDir <- system.file("Choropleth_map", package = "ReKS")
    if (appDir == "") {
        stop(paste0("Could not find Choropleth_map directory. ",
                    "Try re-installing `ReKS`."), call. = FALSE)
    }

    shiny::runApp(appDir, display.mode = "normal")
}
