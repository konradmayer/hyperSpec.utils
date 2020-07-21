
#' Interactive spectral map explorer
#'
#' @description Shiny app to interactively explore spatially resolved \code{hyperSpec} objects.
#' @param hyperspec_obj A \code{hyperSpec} object with variables "x" and "y" in \code{@data}.
#' @param fixed_y logical or numeric; if TRUE the initial y axis range is set to the intensity range of the complete dataset until the axis is manually modified or autoscale activated; alternatively provide a numeric vector of length 2 to set a fixed y range; if FALSE autoscale is activated.
#' @param flip logical; if TRUE the x and y axes are transposed.
#' @param startband numeric; band for the displayed image at startup.
#' @param metavar character; a variable in \code{@data} can be selected as an alternative to the intensity map (e.g. to display cluster results).
#' @export

spcmap_explorer <- function(hyperspec_obj, fixed_y = TRUE, flip = FALSE, startband = 1600, metavar = NULL) {

  # input validation
  if(!is_hyperSpecMap(hyperspec_obj)) {stop("Please provide a hyperSpec object with variables 'x' and 'y' in @data")}
  stopifnot(is.logical(flip))
  stopifnot(is.logical(fixed_y) | (is.numeric(fixed_y) & (length(fixed_y) == 2)))
  stopifnot(is.numeric(startband) & (length(startband) == 1))


  # todo
  # - explicitely call functions from namespace and remove require calls
  # - get rid of double click action
  # - use labels stored in hyperspec object
  # - in far future: add region selections instead of single bands; band height ratios (right click?)


  require(tidyverse)
  require(shiny)
  require(plotly)
  require(hyperSpec)

  global_y_limits <- range(hyperspec_obj[[]])
  if (!is.null(metavar)) {hyperspec_obj@data[,metavar] <- as.numeric(hyperspec_obj@data[,metavar])}
  removebuttons <- c('sendDataToCloud', 'toggleSpikelines', 'resetScale2d', 'hoverClosestCartesian', 'hoverCompareCartesian')

  shinyApp(
    ui = fluidPage(
      tags$style(HTML("#big-heading{background-color: lightgray; padding: 20px; border-radius: 5px; text-align: right; font-size: 30px;}")),
      fluidRow(
        column(6,"Spectral map explorer", id = "big-heading"),
        column(3,
               sliderInput('col_slider', 'intensity range',
                           min = 0, max = 100,
                           value = c(0, 100), round = 0),
               offset = 1
        ),
        column(2, div(actionButton("close", "Close app"),
                      style = "position: relative; top: 50%; -webkit-transform: translateY(50%); -ms-transform: translateY(50%); transform: translateY(50%);")
        )
      ),
      plotlyOutput("plot"),
      plotlyOutput("plotspc")#,
      #verbatimTextOutput("doubleclick"),
    ),
    server = function(input, output, session) {



      # init reactive variables -------------------------------------------------
      click <- reactiveVal() #click on x,y position on map to extract spectrum
      xselect <- reactiveVal() #derived x value from click()
      yselect <- reactiveVal() #derived y value from click()

      band <- reactiveVal() #click on position on spectrum to extract wavenumber

      layout_var <- reactiveVal() # variable to store dynamic layout changes on spectrumplot
      xrange_start <- range(hyperspec_obj@wavelength)
      if (is.numeric(fixed_y) && (length(fixed_y) == 2)){
        yrange_start <- fixed_y
      } else if (fixed_y) {
        yrange_start <- global_y_limits
      } else {
        yrange_start <- NULL
      }
      xrange <- reactiveVal(xrange_start) # start values for spectral plot (before dynamic update through "plotly_relayout")
      yrange <- reactiveVal(value = yrange_start) # start values for spectral plot (before dynamic update through "plotly_relayout")

      layout_var_raster <- reactiveVal() # variable to store dynamic layout changes on raster map plot

      if (flip) {
        xrange_raster_start <- range(hyperspec_obj@data$y)
        yrange_raster_start <- range(hyperspec_obj@data$x)
      } else {
        xrange_raster_start <- range(hyperspec_obj@data$x)
        yrange_raster_start <- range(hyperspec_obj@data$y)
      }
      xrange_raster <- reactiveVal(value = xrange_raster_start) # start values for raster map plot (before dynamic update through "plotly_relayout")
      yrange_raster <- reactiveVal(value = yrange_raster_start) # start values for raster map plot (before dynamic update through "plotly_relayout")




      # raster map plot ---------------------------------------------------------

      output$plot <- renderPlotly({

        # update axis ranges if map is zoomed of panned
        if (!is.null(layout_var_raster())) {
          xrange_raster(c(layout_var_raster()$`xaxis.range[0]`,layout_var_raster()$`xaxis.range[1]`))
          yrange_raster(c(layout_var_raster()$`yaxis.range[0]`,layout_var_raster()$`yaxis.range[1]`))
        }

        # show only one wavenumber either selected by argument or click action (band())
        if (is.null(band())) {
          bandselect <- startband
        } else {
          bandselect <- band()[['x']]
        }

        # subset hyperspec object
        plotdat <- hyperspec_obj[, , bandselect]


        # if metavar is provided as an argument use it instead of band. Assign x y and z to the respective vectors
        if (!is.null(metavar)) {
          z = as.formula(paste0("~", metavar))
          raster_plt_title <- paste0('map for variable: ', metavar)
        } else {
          z = ~spc[, 1]
          raster_plt_title <- paste0('intensity map for wavenumber: ', bandselect)
        }
        if (flip) {
          x = ~y; y = ~x; template <- '<b>X</b>: %{y:.2f}<br><b>Y</b>: %{x:.2f}<br><b>Intensity</b>: %{z:.2f}<extra></extra>'
        } else {
          x = ~x; y = ~y; template <- '<b>X</b>: %{x:.2f}<br><b>Y</b>: %{y:.2f}<br><b>Intensity</b>: %{z:.2f}<extra></extra>'}

        colsliderpos <- input$col_slider
        colramp <- c(rep("#440154FF", colsliderpos[1]),
                     viridis::viridis(diff(colsliderpos)),
                     rep("#FDE725FF", 100 - colsliderpos[2]))

        # plot
        raster_plt <- plotdat %>%
          as.data.frame() %>%
          plotly::plot_ly(x = x, y = y, z = z, source = 'A', type = 'heatmap',
                          hovertemplate = template, colors = colramp) %>%
          plotly::layout(xaxis = list(range = xrange_raster(), autorange = 'FALSE',
                                      scaleanchor = "y"),
                         yaxis = list(range = yrange_raster(), autorange = 'FALSE'),
                         title = list(text = raster_plt_title,
                                      x = 0.1, font = list(size = 15))) %>%
          config(scrollZoom = TRUE, doubleClick = FALSE, displaylogo = FALSE,
                 modeBarButtonsToRemove = removebuttons) %>%
          colorbar(title = "Intensity")


      })


      # spectrum plot -----------------------------------------------------------

      output$plotspc <- renderPlotly({

        if (!is.null(click())) {

          # on click plot selected spectrum (x and y coordinates have to be equal up to a set tolerance)
          tolerance <- 0.00001
          selection <- which((abs(hyperspec_obj@data$x - xselect()) <= tolerance) &
                               (abs(hyperspec_obj@data$y - yselect()) <= tolerance))



          # static plot using ggplot
          p_tmp <- hyperspec_obj[selection, , ] %>%
            qplotspc() +
            theme_minimal()

          # update x and y range
          if (!is.null(layout_var()$`xaxis.range[0]`)) {
            xrange(c(layout_var()$`xaxis.range[0]`,layout_var()$`xaxis.range[1]`))
          }
          if (!is.null(layout_var()$`yaxis.range[0]`)) {
            yrange(c(layout_var()$`yaxis.range[0]`,layout_var()$`yaxis.range[1]`))
          }

          # autorange: when autoscale button is pressed, autoscale is applied to all spectra until a scaling is done (zoom, pan) - then the altered scale applied
          if (!is.null(layout_var()$`xaxis.autorange`)) {
            xrange(range(p_tmp$data$.wavelength))
          }
          if (!is.null(layout_var()$`yaxis.autorange`)) {
            yrange(range(p_tmp$data$spc))
          }

          # create plotly object from static plot
          out <- ggplotly(p_tmp, source = 'B', tooltip = "text") %>%
            style(text = paste0('Wavenumber: ', round(p_tmp$data$.wavelength, 2),
                                '<br>', 'Intensity: ', round(p_tmp$data$spc, 2))) %>%
            plotly::layout(xaxis = list(range = xrange(), autorange = 'FALSE'),
                           yaxis = list(range = yrange(), autorange = 'FALSE',
                                        fixedrange = FALSE),
                           title = list(text = paste0('spectrum at position x: ',
                                                      round(xselect(), 2),
                                                      ', y: ',
                                                      round(yselect(), 2)),
                                        x = 0.1, font = list(size = 15))) %>%
            config(scrollZoom = TRUE, doubleClick = FALSE, displaylogo = FALSE,
                   modeBarButtonsToRemove = removebuttons)

          # draw vertical line ###bugfix necessary - update when autorange is pressed (otherwise too short if autorange zooms out)
          if (!is.null(band()[['x']])) {
            out <- out %>%
              add_segments(x = band()[['x']], xend = band()[['x']],
                           y = yrange()[1], yend = yrange()[2])
          }

          # return plot
          out

        }
      })

      # update reactive variables when map is clicked or band selected
      click <- reactive({event_data("plotly_click", source = "A")})
      observeEvent(click(), {
        layout_var(event_data("plotly_relayout", source = "B"))
        clicked <- c(as.numeric(click()[['x']]),  as.numeric(click()[['y']]))
        if (flip) {
          xselect(clicked[[2]]); yselect(clicked[[1]])
        } else {
          xselect(clicked[[1]]); yselect(clicked[[2]])
        }
      })

      band <- reactive({event_data("plotly_click", source = "B")})
      observeEvent(band(), {
        layout_var(event_data("plotly_relayout", source = "B"))
        layout_var_raster(event_data("plotly_relayout", source = "A"))
        updateSliderInput(session, 'col_slider', value = c(0, 100))
      })

      # stop app if browser is closed
      session$onSessionEnded(stopApp)

      # close button
      observeEvent(input$close, {
        stopApp()
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
