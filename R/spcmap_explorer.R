
#' Interactive spectral map explorer
#'
#' @description Shiny app to interactively explore spatially resolved \code{hyperSpec} objects.
#' @param hyperspec_obj A \code{hyperSpec} object with variables "x" and "y" in \code{@data}.
#' @param fixed_y logical or numeric; if TRUE the initial y axis range is set to the intensity range of the complete dataset until the axis is manually modified or autoscale activated; alternatively provide a numeric vector of length 2 to set a fixed y range; if FALSE autoscale is activated.
#' @param flip logical; if TRUE the x and y axes are transposed.
#' @param startband numeric; band for the displayed image at startup.
#' @param metavar character; a variable in \code{@data} can be selected as an alternative to the intensity map (e.g. to display cluster results).
#' @importFrom magrittr "%>%"
#' @export

spcmap_explorer <- function(hyperspec_obj, fixed_y = TRUE, flip = FALSE, startband = 1600, metavar = NULL) {

  # input validation
  if (!is_hyperSpecMap(hyperspec_obj)) {
    stop("Please provide a hyperSpec object with variables 'x' and 'y' in @data")
  }
  stopifnot(is.logical(flip))
  stopifnot(is.logical(fixed_y) | (is.numeric(fixed_y) & (length(fixed_y) == 2)))
  stopifnot(is.numeric(startband) & (length(startband) == 1))


  # todo
  # - get rid of the event_register warning
  # - consistent naming (use spc, wl for specplot)
  # - in far future: add region selections instead of single bands; band height ratios (right click?)


  if (!is.null(metavar)) {
    hyperspec_obj@data[, metavar] <- as.numeric(hyperspec_obj@data[, metavar])
  }
  removebuttons <- c("sendDataToCloud", "toggleSpikelines", "resetScale2d", "hoverClosestCartesian", "hoverCompareCartesian")

  shiny::shinyApp(
    ui = shiny::fluidPage(
      shiny::tags$style(shiny::HTML("#big-heading{background-color: lightgray; padding: 20px 40px 20px 20px; border-radius: 0px 40px 10px 0px; text-align: right; font-size: 30px;}")),
      shiny::fluidRow(
        shiny::column(6, "Spectral map explorer", id = "big-heading"),
        shiny::column(3,
          shiny::sliderInput("col_slider", "colorramp range",
            min = 0, max = 100,
            value = c(0, 100), round = 0
          ),
          offset = 1
        ),
        shiny::column(2, shiny::div(shiny::actionButton("close", "Close app"),
          style = "position: relative; top: 50%; -webkit-transform: translateY(50%); -ms-transform: translateY(50%); transform: translateY(50%);"
        ))
      ),
      plotly::plotlyOutput("plot"),
      plotly::plotlyOutput("plotspc") # ,
      # verbatimTextOutput("doubleclick"),
    ),
    server = function(input, output, session) {

      # fix labels -------------------------------------------------------------
      if (any(vapply(hyperspec_obj@label, is.expression, logical(1)))) {
        warning("One or multiple labels of the provided of hyperSpec object contained an expression which plotly can't handle - it is now replaced by a default label which may be not applicable")
      }
      xlabel <- ifelse((is.expression(hyperspec_obj@label$x) | !any("x" %in% names(hyperspec_obj@label)) | is.null(hyperspec_obj@label$x)),
        "x",
        hyperspec_obj@label$x
      )
      ylabel <- ifelse((is.expression(hyperspec_obj@label$y) | !any("y" %in% names(hyperspec_obj@label)) | is.null(hyperspec_obj@label$y)),
        "y",
        hyperspec_obj@label$y
      )
      wllabel <- ifelse((is.expression(hyperspec_obj@label$.wavelength) | !any(".wavelength" %in% names(hyperspec_obj@label)) | is.null(hyperspec_obj@label$.wavelength)),
        "wavenumber",
        hyperspec_obj@label$.wavelength
      )
      spclabel <- ifelse((is.expression(hyperspec_obj@label$spc) | !any("spc" %in% names(hyperspec_obj@label)) | is.null(hyperspec_obj@label$spc)),
        "intensity",
        hyperspec_obj@label$spc
      )

      hyperspec_obj@label$x <- xlabel
      hyperspec_obj@label$y <- ylabel
      hyperspec_obj@label$.wavelength <- wllabel
      hyperspec_obj@label$spc <- spclabel

      # init reactive variables ------------------------------------------------

      xselect <- shiny::reactiveVal() # derived x value from click()
      yselect <- shiny::reactiveVal() # derived y value from click()

      layout_var <- shiny::reactiveVal() # variable to store dynamic layout changes on spectrumplot
      layout_var_raster <- shiny::reactiveVal() # variable to store dynamic layout changes on raster map plot

      xrange_start <- range(hyperspec_obj@wavelength)
      global_y_limits <- range(hyperspec_obj[[]])

      if (is.numeric(fixed_y) && (length(fixed_y) == 2)) {
        yrange_start <- fixed_y
      } else if (fixed_y) {
        yrange_start <- global_y_limits
      } else {
        yrange_start <- NULL
      }
      xrange <- shiny::reactiveVal(value = xrange_start) # start values for spectral plot (before dynamic update through "plotly_relayout")
      yrange <- shiny::reactiveVal(value = yrange_start) # start values for spectral plot (before dynamic update through "plotly_relayout")


      if (flip) {
        xrange_raster_start <- range(hyperspec_obj@data$y)
        yrange_raster_start <- range(hyperspec_obj@data$x)
      } else {
        xrange_raster_start <- range(hyperspec_obj@data$x)
        yrange_raster_start <- range(hyperspec_obj@data$y)
      }
      xrange_raster <- shiny::reactiveVal(value = xrange_raster_start) # start values for raster map plot (before dynamic update through "plotly_relayout")
      yrange_raster <- shiny::reactiveVal(value = yrange_raster_start) # start values for raster map plot (before dynamic update through "plotly_relayout")




      # raster map plot ---------------------------------------------------------

      output$plot <- plotly::renderPlotly({

        # update axis ranges if map is zoomed of panned
        if (!is.null(layout_var_raster())) {
          xrange_raster(c(layout_var_raster()$`xaxis.range[0]`, layout_var_raster()$`xaxis.range[1]`))
          yrange_raster(c(layout_var_raster()$`yaxis.range[0]`, layout_var_raster()$`yaxis.range[1]`))
        }

        # show only one wavenumber either selected by argument or click action (band())
        if (is.null(band())) {
          bandselect <- startband
        } else {
          bandselect <- band()[["x"]]
        }

        # subset hyperspec object
        plotdat <- hyperspec_obj[, , bandselect]


        # if metavar is provided as an argument use it instead of band. Assign x y and z to the respective vectors
        if (!is.null(metavar)) {
          z <- stats::as.formula(paste0("~", metavar))
          raster_plt_title <- paste0("map for variable: ", metavar)
        } else {
          z <- ~ spc[, 1]
          raster_plt_title <- paste0("intensity map for ", wllabel, " ", bandselect)
        }
        if (flip) {
          x <- ~y
          y <- ~x
          template <- paste0("<b>", ylabel, "</b>: %{x:.2f}<br><b>", xlabel, "</b>: %{y:.2f}<br><b>", spclabel, "</b>: %{z:.2f}<extra></extra>")
        } else {
          x <- ~x
          y <- ~y
          template <- paste0("<b>", xlabel, "</b>: %{x:.2f}<br><b>", ylabel, "</b>: %{y:.2f}<br><b>", spclabel, "</b>: %{z:.2f}<extra></extra>")
        }

        colsliderpos <- input$col_slider
        colramp <- c(
          rep("#440154FF", colsliderpos[1]),
          viridis::viridis(diff(colsliderpos)),
          rep("#FDE725FF", 100 - colsliderpos[2])
        )

        # plot
        raster_plt <- plotdat %>%
          as.data.frame() %>%
          plotly::plot_ly(
            x = x, y = y, z = z, source = "A", type = "heatmap",
            hovertemplate = template, colors = colramp
          ) %>%
          plotly::layout(
            xaxis = list(
              range = xrange_raster(), autorange = "FALSE",
              scaleanchor = "y", title = ifelse(flip, ylabel, xlabel)
            ),
            yaxis = list(
              range = yrange_raster(), autorange = "FALSE",
              title = ifelse(flip, xlabel, ylabel)
            ),
            title = list(
              text = raster_plt_title,
              x = 0.1, font = list(size = 15)
            )
          ) %>%
          plotly::config(
            scrollZoom = TRUE, doubleClick = FALSE, displaylogo = FALSE,
            modeBarButtonsToRemove = removebuttons
          ) %>%
          plotly::colorbar(title = spclabel)
      })


      # spectrum plot -----------------------------------------------------------

      output$plotspc <- plotly::renderPlotly({
        if (!is.null(click())) {

          # on click plot selected spectrum (x and y coordinates have to be equal up to a set tolerance)
          tolerance <- 0.00001
          selection <- which((abs(hyperspec_obj@data$x - xselect()) <= tolerance) &
            (abs(hyperspec_obj@data$y - yselect()) <= tolerance))



          # static plot using ggplot
          p_tmp <- hyperspec_obj[selection, , ] %>%
            hyperSpec::qplotspc() +
            ggplot2::theme_minimal()

          # update x and y range
          if (!is.null(layout_var()$`xaxis.range[0]`)) {
            xrange(c(layout_var()$`xaxis.range[0]`, layout_var()$`xaxis.range[1]`))
          }
          if (!is.null(layout_var()$`yaxis.range[0]`)) {
            yrange(c(layout_var()$`yaxis.range[0]`, layout_var()$`yaxis.range[1]`))
          }

          # autorange: when autoscale button is pressed, autoscale is applied to all spectra until a scaling is done (zoom, pan) - then the altered scale applied
          if (!is.null(layout_var()$`xaxis.autorange`)) {
            xrange(range(p_tmp$data$.wavelength))
          }
          if (!is.null(layout_var()$`yaxis.autorange`)) {
            yrange(range(p_tmp$data$spc))
          }

          # create plotly object from static plot
          out <- plotly::ggplotly(p_tmp, source = "B", tooltip = "text") %>%
            plotly::style(text = paste0(
              wllabel, ": ", round(p_tmp$data$.wavelength, 2),
              "<br>", spclabel, ": ", round(p_tmp$data$spc, 2)
            )) %>%
            plotly::layout(
              xaxis = list(range = xrange(), autorange = "FALSE"),
              yaxis = list(
                range = yrange(), autorange = "FALSE",
                fixedrange = FALSE
              ),
              title = list(
                text = paste0(
                  "spectrum at position ", xlabel, ": ",
                  round(xselect(), 2),
                  ", ", ylabel, ": ",
                  round(yselect(), 2)
                ),
                x = 0.1, font = list(size = 15)
              )
            ) %>%
            plotly::config(
              scrollZoom = TRUE, doubleClick = FALSE, displaylogo = FALSE,
              modeBarButtonsToRemove = removebuttons
            )

          # draw vertical line ###bugfix necessary - update when autorange is pressed (otherwise too short if autorange zooms out)
          if (!is.null(band()[["x"]])) {
            out <- out %>%
              plotly::add_segments(
                x = band()[["x"]], xend = band()[["x"]],
                y = yrange()[1], yend = yrange()[2]
              )
          }

          # return plot
          out
        }
      })

      # update reactive variables when map is clicked or band selected
      click <- shiny::reactive({
        plotly::event_data("plotly_click", source = "A")
      })

      shiny::observeEvent(click(), {
        layout_var(plotly::event_data("plotly_relayout", source = "B"))
        clicked <- c(as.numeric(click()[["x"]]), as.numeric(click()[["y"]]))
        if (flip) {
          xselect(clicked[[2]])
          yselect(clicked[[1]])
        } else {
          xselect(clicked[[1]])
          yselect(clicked[[2]])
        }
      })

      # as the event_data call on source B is happening here before the plot is called a warning is thrown
      band <- shiny::reactive({
        plotly::event_data("plotly_click", source = "B")
      })

      shiny::observeEvent(band(), {
        layout_var(plotly::event_data("plotly_relayout", source = "B"))
        layout_var_raster(plotly::event_data("plotly_relayout", source = "A"))
        shiny::updateSliderInput(session, "col_slider", value = c(0, 100))
      })

      # stop app if browser is closed
      session$onSessionEnded(shiny::stopApp)

      # close button
      shiny::observeEvent(input$close, {
        shiny::stopApp()
      })


      # output$info <- renderPrint({
      #   d <- event_data("plotly_selected")
      #   if (is.null(d)) "Infos appear here (double-click to clear)" else d
      # })
    }
  )
}


# testdat <- readr::read_rds('data/shiny-testdata.rds')
# testdat[[(testdat$y %in% c(-794, -795) & (testdat$x %in% 9:11)), ]] <- NA
# spcmap_explorer(testdat, fixed_y = c(-10, 500), flip = TRUE)
# spcmap_explorer(testdat, fixed_y = c(-10, 500), flip = TRUE, metavar = "cluster")
# spcmap_explorer(testdat)
#
# #plotly cant deal with expressions - therefore i change the labels in chondro
# chondro2 <- chondro
# chondro2@label$.wavelength <- "wavenumber"
# chondro2@label$spc <- "intensity"
# chondro2@label$clusters <- "clusters"
# spcmap_explorer(chondro2, startband = 1450)
# spcmap_explorer(chondro2, metavar = 'clusters')
