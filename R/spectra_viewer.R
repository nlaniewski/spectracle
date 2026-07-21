spectra_viewer <- function(instrument = c("Aurora")){
  ##
  traces <- readRDS(list.files(
    system.file("extdata/spectra_library", package = "spectracle"),
    full.names = T,
    pattern = paste0(match.arg(instrument), ".*long")
  ))
  ##
  title.instrument <- traces[,unique(instrument)]
  title.configuration <- traces[, unique(configuration)]
  subtitle <- sprintf("%s::%s", title.instrument, title.configuration)
  ##
  ui <- shiny::fluidPage(
    shiny::titlePanel("Spectracle Spectra Viewer", windowTitle = "Spectracle Spectra Viewer"),
    ##
    # shinyWidgets::pickerInput(
    #   inputId = "selected_laser",
    #   label = "Laser:",
    #   choices = traces[,levels(laser)],
    #   options = list(
    #     `actions-box` = TRUE,
    #     `live-search` = TRUE
    #   ),
    #   multiple = TRUE
    # ),
    ##
    shinyWidgets::pickerInput(
      inputId = "selected_detector",
      label = "Detector:",
      choices = traces[, unique(detector)],
      selected = NULL,
      options = list(
        `actions-box` = TRUE,
        `live-search` = TRUE
      ),
      multiple = TRUE
    ),
    ##
    shinyWidgets::pickerInput(
      inputId = "selected_fluor",
      label = "Fluorochrome:",
      choices = traces[, unique(fluorochrome)],
      selected = traces[, unique(fluorochrome)][1],
      options = list(
        `actions-box` = TRUE,
        `live-search` = TRUE
      ),
      multiple = TRUE
    ),
    ##
    shiny::mainPanel(
      plotly::plotlyOutput(outputId = "spectral_trace")
    )
  )
  ##
  server <- function(input, output, session) {
    ##
    shiny::observeEvent(input$selected_fluor, {
      shiny::req(input$selected_fluor)
      shinyWidgets::updatePickerInput(
        session, "selected_laser", selected = character(0)
      )
      shinyWidgets::updatePickerInput(
        session, "selected_detector", selected = character(0)
      )
    }, ignoreNULL = TRUE)
    ##
    shiny::observeEvent(input$selected_detector, {
      shiny::req(input$selected_detector)
      shinyWidgets::updatePickerInput(
        session, "selected_laser", selected = character(0)
      )
      shinyWidgets::updatePickerInput(
        session, "selected_fluor", selected = character(0)
      )
    }, ignoreNULL = TRUE)
    ##
    # shiny::observeEvent(input$selected_laser, {
    #   shiny::req(input$selected_laser)
    #   shinyWidgets::updatePickerInput(
    #     session, "selected_detector", selected = character(0)
    #   )
    #   shinyWidgets::updatePickerInput(
    #     session, "selected_fluor", selected = character(0)
    #   )
    # }, ignoreNULL = TRUE)
    ##
    output$spectral_trace <- plotly::renderPlotly({
      ##
      if(shiny::isTruthy(input$selected_fluor)){
        ##
        n <- traces[
          i = fluorochrome %in% input$selected_fluor,
          j = data.table::uniqueN(fluorochrome)
        ]
        ##
        p <- spectral_trace_plotly(
          .data = traces[fluorochrome %in% input$selected_fluor],
          n,
          subtitle
        )
      }else if(shiny::isTruthy(input$selected_detector)){
        ##
        n <- traces[
          i = detector %in% input$selected_detector,
          j = data.table::uniqueN(fluorochrome)
        ]
        ##
        p <- spectral_trace_plotly(
          .data = traces[detector %in% input$selected_detector],
          n,
          subtitle
        )
      }else if(shiny::isTruthy(input$selected_laser)){
        ##
        n <- traces[
          i = laser %in% input$selected_laser,
          j = data.table::uniqueN(fluorochrome)
        ]
        ##
        p <- spectral_trace_plotly(
          .data = traces[laser %in% input$selected_laser],
          n,
          subtitle
        )
      }
    })
  }
  ##
  shiny::shinyApp(ui, server)
  ##
}
##
spectral_trace_plotly <- function(.data, n, subtitle){
  p <- plotly::plot_ly(
    data = .data,
    x = ~Detector,
    y = ~value,
    color = ~fluorochrome,
    colors = colors.kelly[seq_len(n)],
    type = "scatter",
    mode = "lines"
  )
  p <- plotly::layout(p, xaxis = list(type = 'category'))
  p <- plotly::layout(
    p,
    title = sprintf("<b>Spectral Trace(s)</b><br><sup>%s</sup>", subtitle),
    xaxis = list(tickmode = 'linear', dtick = 1, tickangle = 270),
    yaxis = list(title = "Emission (Normalized [0,1])"),
    plot_bgcolor = "white",
    legend = list(
      orientation = "h",
      xanchor = "center",
      x = 0.5,
      y = -0.2
    ),
    showlegend = T,
    margin = list(t = 50)
  )
  ##
  return(p)
}
