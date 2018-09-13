# get_eurostat_geospatial <- function(output_class = "sf", resolution = "60",
#                                     nuts_level = "all", year = "2016",
#                                     cache = TRUE, update_cache = FALSE,
#                                     cache_dir = NULL) {
#     eurostat_geodata_60_2016 <- NULL
#     LEVL_CODE <- NULL
#     data("eurostat_geodata_60_2013",
#          envir = environment(),
#          package = "eurostat")
#     resolution <- as.character(resolution)
#     resolution <- gsub("^0+", "", resolution)
#     if (!as.numeric(resolution) %in% c(1, 3, 10, 20, 60)) {
#         stop("Resolution should be one of 01, 1, 03, 3, 10, 20, 60")
#     }
#     resolution <- gsub("^1$", "01", resolution)
#     resolution <- gsub("^3$", "03", resolution)
#     year <- as.character(year)
#     if (!as.numeric(year) %in% c(2003, 2006, 2010, 2013, 2016)) {
#         stop("Year should be one of 2003, 2006, 2010, 2013, 2016")
#     }
#     if (as.numeric(year) == 2003 & as.numeric(resolution) ==
#         60) {
#         stop(paste0("NUTS 2003 is not provided at 1:60 million resolution. ",
#                     "Try 1:1 million, 1:3 million, 1:10 million or ",
#                     "1:20 million"))
#     }
#     message(paste0("\nCOPYRIGHT NOTICE\n\nWhen data downloaded from this",
#                    " page \n<http://ec.europa.eu/eurostat/web/gisco/geodata",
#                    "/reference-data/administrative-units-statistical-units>\n",
#                    "is used in any printed or electronic publication, \n",
#                    "in addition to any other provisions \napplicable to the ",
#                    "whole Eurostat website, \ndata source will have to be ",
#                    "acknowledged \nin the legend of the map and \nin the ",
#                    "introductory page of the publication \nwith the following ",
#                    "copyright notice:\n\n- EN: (C) EuroGeographics for the ",
#                    "administrative boundaries\n- FR: (C) EuroGeographics pour ",
#                    "les limites administratives\n- DE: (C) EuroGeographics ",
#                    "bezuglich der Verwaltungsgrenzen\n\nFor publications in ",
#                    "languages other than \nEnglish, French or German, \n",
#                    "the translation of the copyright notice \nin the language ",
#                    "of the publication shall be used.\n\nIf you intend to use ",
#                    "the data commercially, \nplease contact EuroGeographics ",
#                    "for \ninformation regarding their licence agreements.",
#                    "\n          "))
#     if (resolution == "60" & year == 2013) {
#         if (nuts_level %in% c("all")) {
#             shp <- eurostat_geodata_60_2013
#         }
#         if (nuts_level == "0")
#             shp <- subset(eurostat_geodata_60_2013, LEVL_CODE == 0)
#         if (nuts_level == "1")
#             shp <- subset(eurostat_geodata_60_2013, LEVL_CODE == 1)
#         if (nuts_level == "2")
#             shp <- subset(eurostat_geodata_60_2013, LEVL_CODE == 2)
#         if (nuts_level == "3")
#             shp <- subset(eurostat_geodata_60_2013, LEVL_CODE == 3)
#         if (output_class == "df") {
#             nuts_sp <- as(shp, "Spatial")
#             nuts_sp$id <- row.names(nuts_sp)
#             nuts_ff <- broom::tidy(nuts_sp)
#             shp <- left_join(nuts_ff, nuts_sp@data)
#         }
#         if (output_class == "spdf") {
#             shp <- as(shp, "sp::Spatial")
#         }
#     }
#     else {
#         if (cache) {
#             update_cache <- update_cache | getOption("eurostat_update",
#                                                      FALSE)
#             if (is.null(cache_dir)) {
#                 cache_dir <- getOption("eurostat_cache_dir",
#                                        NULL)
#                 if (is.null(cache_dir)) {
#                     cache_dir <- file.path(tempdir(), "eurostat")
#                     if (!file.exists(cache_dir))
#                         dir.create(cache_dir)
#                 }
#             }
#             else {
#                 if (!file.exists(cache_dir)) {
#                     stop("The folder ", cache_dir, " does not exist")
#                 }
#             }
#             cache_file <- file.path(cache_dir, paste0(output_class, resolution,
#                                                       nuts_level, year,
#                                                       ".RData"))
#         }
#         if (!cache || update_cache || !file.exists(cache_file)) {
#             if (nuts_level %in% c("0", "all")) {
#                 resp <- httr::GET(paste0("http://ec.europa.eu/eurostat/cache/",
#                                          "GISCO/distribution/v2/nuts/geojson/",
#                                          "NUTS_RG_", resolution, "M_", year,
#                                          "_4258_LEVL_0.geojson"))
#                 nuts0 <- sf::st_read(content(resp, as = "text"),
#                                      stringsAsFactors = FALSE, quiet = TRUE)
#             }
#             if (nuts_level %in% c("1", "all")) {
#                 resp <- httr::GET(paste0("http://ec.europa.eu/eurostat/cache/",
#                                          "GISCO/distribution/v2/nuts/geojson/",
#                                          "NUTS_RG_", resolution, "M_", year,
#                                          "_4258_LEVL_1.geojson"))
#                 nuts1 <- sf::st_read(content(resp, as = "text"),
#                                      stringsAsFactors = FALSE, quiet = TRUE)
#             }
#             if (nuts_level %in% c("2", "all")) {
#                 resp <- httr::GET(paste0("http://ec.europa.eu/eurostat/cache/",
#                                          "GISCO/distribution/v2/nuts/geojson/",
#                                          "NUTS_RG_", resolution, "M_", year,
#                                          "_4258_LEVL_2.geojson"))
#                 nuts2 <- sf::st_read(content(resp, as = "text"),
#                                      stringsAsFactors = FALSE, quiet = TRUE)
#             }
#             if (nuts_level %in% c("3", "all")) {
#                 resp <- httr::GET(paste0("http://ec.europa.eu/eurostat/cache/",
#                                          "GISCO/distribution/v2/nuts/geojson/",
#                                          "NUTS_RG_", resolution, "M_", year,
#                                          "_4258_LEVL_3.geojson"))
#                 nuts3 <- sf::st_read(content(resp, as = "text"),
#                                      stringsAsFactors = FALSE, quiet = TRUE)
#             }
#             if (nuts_level %in% c("all")) {
#                 shp <- rbind(nuts0, nuts1, nuts2, nuts3)
#             }
#             if (nuts_level == "0")
#                 shp <- nuts0
#             if (nuts_level == "1")
#                 shp <- nuts1
#             if (nuts_level == "2")
#                 shp <- nuts2
#             if (nuts_level == "3")
#                 shp <- nuts3
#             if (output_class == "df") {
#                 nuts_sp <- as(shp, "Spatial")
#                 nuts_sp$id <- row.names(nuts_sp)
#                 nuts_ff <- broom::tidy(nuts_sp)
#                 shp <- left_join(nuts_ff, nuts_sp@data)
#             }
#             if (output_class == "spdf") {
#                 shp <- as(shp, "Spatial")
#             }
#         }
#     }
#     if (resolution != "60" & year != 2013) {
#         if (cache & file.exists(cache_file)) {
#             cf <- path.expand(cache_file)
#             message(paste("Reading cache file", cf))
#             load(file = cache_file)
#             if (output_class == "sf")
#                 message(paste("sf at resolution 1:",
#                               resolution, " from year ", year,
#                               " read from cache file: ", cf))
#             if (output_class == "df")
#                 message(paste("data_frame at resolution 1:",
#                               resolution, " from year ", year,
#                               " read from cache file: ", cf))
#             if (output_class == "spdf")
#                 message(paste("SpatialPolygonDataFrame at resolution 1:",
#                               resolution, " from year ", year,
#                               " read from cache file: ", cf))
#         }
#         if (cache && (update_cache || !file.exists(cache_file))) {
#             save(shp, file = cache_file)
#             if (output_class == "sf")
#                 message(paste("sf at resolution 1:",
#                               resolution, " cached at: ",
#                               path.expand(cache_file)))
#             if (output_class == "df")
#                 message(paste("data_frame at resolution 1:",
#                               resolution, " cached at: ",
#                               path.expand(cache_file)))
#             if (output_class == "spdf")
#                 message(paste("SpatialPolygonDataFrame at resolution 1:",
#                               resolution, " cached at: ",
#                               path.expand(cache_file)))
#         }
#     }
#     if (resolution == "60" & year == 2013) {
#         if (output_class == "sf")
#             message(paste("sf at resolution 1:60 read from local file"))
#         if (output_class == "df")
#             message(paste("data_frame at resolution 1:60 read from local file"))
#         if (output_class == "spdf")
#             message(paste("SpatialPolygonDataFrame at resolution 1:60 read ",
#                           "from local file"))
#     }
#     message("\n# --------------------------\nHEADS UP!!\n\nFunction now ",
#             "returns the data in 'sf'-class (simple features) \nby default ",
#             "which is different \nfrom previous behaviour's ",
#             "'SpatialPolygonDataFrame'. \n\nIf you prefer either ",
#             "'SpatialPolygonDataFrame' or \nfortified 'data_frame' ",
#             "(for ggplot2::geom_polygon), \nplease specify it explicitly to ",
#             "'output_class'-argument!\n\n# --------------------------",
#             "          \n          ")
#     shp$geo <- shp$NUTS_ID
#     return(shp)
# }

msrs <- list("Variety (total)"   = "Variety.total.1_4",
             "Related variety"   = "Variety.related.1_4",
             "Unrelated variety" = "Variety.unrelated.1_4",
             "Coherence"         = "Coherence.4",
             "Reg. Complexity"   = "Complexity.4",
             "Reg. Fitness (log)"= "Fitness.4")

# Define UI for application that draws a histogram
ui <- shiny::fluidPage(

    # Application title
    shiny::titlePanel("EU Regional Knowledge Space"),

    # Sidebar with a slider input for number of bins
    shiny::sidebarLayout(
        shiny::sidebarPanel(
            shiny::selectInput(inputId = "geo_region",
                               label = h3("Geographical region"),
                               choices = list("European Union"   = "eu")),
            shiny::sliderInput(inputId = "geo_level",
                               label = h3("Geographical level of aggregation"),
                               min = 2,
                               max = 3,
                               value = 2,
                               step = 1),
            shiny::sliderInput(inputId = "year",
                               label = h3("Year"),
                               min = 2000,
                               max = 2012,
                               value = 13,
                               sep = ""),
            shiny::selectInput(inputId = "measure",
                               label = h3("Measure"),
                               choices = msrs)
        ),

        # Show a plot of the generated distribution
        shiny::mainPanel(
            shiny::plotOutput(outputId = "map",
                              height = "600px")
        )
    )
)

server <- function(input, output) {

    cache_map <- FALSE

    map_dta <- reactive({
        if (input$geo_region == "eu") {
            geo_dta <-
                eurostat::get_eurostat_geospatial(resolution = "60",
                                                  nuts_level = input$geo_level,
                                                  year = 2013,
                                                  output_class = "sf")
        }

        basic_src <- paste0("~/Documents/UniTo/Progetto_di_Ricerca/",
                            "Resilience_EU_Regions/Relatedness_measures/")
        if (input$geo_level == 2) {
            load(paste0(basic_src, "Relatedness.NUTS2.RData"))
            dta <- as.data.frame(relatedness_nuts2)
        } else{
            if (input$geo_level == 3) {
                load(paste0(basic_src, "Relatedness.NUTS3.RData"))
                dta <- as.data.frame(relatedness_nuts3)
            }
        }
        whch <- which(dta[, "Year"] == input$year)
        dta_subset <- dta[whch, c(paste0("NUTS", input$geo_level),
                                  input$measure)]
        colnames(dta_subset) <- c("geo", "measure")

        if (input$measure == "Fitness.4") {
            dta_subset[, "measure"] <- log(dta_subset[, "measure"])
        }
        dta_subset[, "measure"] <- scale(dta_subset[, "measure"])
        dta_subset[, "measure"] <- scales::rescale(dta_subset[, "measure"],
                                                   c(0, 1))

        map_dta <- merge(geo_dta, dta_subset, all.x = TRUE)

        return(map_dta)
    })

    output$map <- shiny::renderPlot({
        clrs <- RColorBrewer::brewer.pal(n = 5, name = "YlGn")
        cptn <- paste("(C) EuroGeographics for the administrative",
                      "boundaries.\n Map produced in R with a help",
                      "from Eurostat-package <github.com/ropengov/",
                      "eurostat/>")
        plt_ttl <- paste(input$measure, "in", input$year)
        lgnd_ttl <- input$measure

        m <- ggplot2::ggplot(data = map_dta()) +
            ggplot2::geom_sf(ggplot2::aes(fill = measure),
                             color = "dim grey", size = .1) +
            ggplot2::scale_fill_gradientn(colours = clrs,
                                          na.value = "gray90") +
            ggplot2::guides(fill = ggplot2::guide_legend(reverse=T,
                                                         title = lgnd_ttl)) +
            ggplot2::labs(title = plt_ttl,
                          caption = cptn) +
            ggplot2::theme_light() +
            ggplot2::theme(legend.position = c(.85, .85)) +
            ggplot2::coord_sf(xlim = c(-12, 44), ylim = c(35, 70))

        return(m)
    })
}

# Run the application
shiny::shinyApp(ui = ui, server = server)
