


HLIRApp <- function(){
  HLIR_func<-function(times,y,parms) {
    dH<- -parms["beta"]*y["H"]*y["I"]
    dL<- parms["beta"]*y["H"]*y["I"] - parms["omega"]*y["L"]
    dI<- parms["omega"]*y["L"] - parms["mu"]*y["I"]
    dR<- parms["mu"]*y["I"]
    list(c(dH,dL,dI,dR)) }


  test1_fun = function(beta, omega, mu, lat, inf, rem, time){
    parms<-c("beta"=beta, "omega"=omega, "mu"=mu)
    state<-c("H"=1000 - (lat+inf+rem), "L"=lat, "I"=inf, "R"=rem)
    times<-seq(0,time,1)

    deSolve::lsoda(y=state, times=times, func=HLIR_func, parms=parms) %>%
      as.data.frame(.) %>%
      mutate(To_I = L+I+R) %>%
      tidyr::pivot_longer(cols = -time, names_to = "var", values_to = "value") %>%
      mutate(var = factor(var, levels = c("H",  "L", "I","R", "To_I"),
                          labels = c("Healthy", "Latently", "Infectious",  "Removed", "Total Infected"))) %>%
      rename("Time" = time)             }


  HLIR_plot = function(data){
    ggplot2::ggplot(data, ggplot2::aes(x = Time, y = value, color = var)) +
      ggplot2::geom_line(size = 1) +
      ggplot2::labs(y = "Density", x = "Time", color = "Individuals") +

      ggplot2::scale_color_manual(values = c(
        "#00A087FF",
        "#E64B35FF",
        "#8491B4FF",
        "#3C5488FF",
        "#DC0000FF"
      )) +
      ggplot2::theme(
        panel.background = ggplot2::element_blank() ,
        panel.border = ggplot2::element_blank(),
        panel.grid.major = ggplot2::element_line(colour = "grey75"),

        axis.title.y = ggplot2::element_text(face = "bold", size = 14, hjust = .9),
        axis.title.x = ggplot2::element_text(face = "bold", size = 14, hjust = .1),

        axis.text  = ggplot2::element_text(colour = "black"),

        legend.key = ggplot2::element_rect(fill = "white"),
        legend.position = "top",
        legend.justification = "right",
        legend.title =  ggplot2::element_text(face = "bold", size = 12, hjust = .9),
        legend.text = ggplot2::element_text(colour = "black", size = 12)
      ) +
      ggplot2::guides(colour = ggplot2::guide_legend(title.position = "top", title.hjust = 0.1)) }

  # User interface
  ui <- fluidPage(

    theme = shinytheme("readable"),
    # Application title
    titlePanel("HLIR (SEIR) model"),

    # Sidebar with a shiny::slider input for number of bins
    sidebarLayout(
      sidebarPanel(
        numericInput(inputId = "inf", label = "# Infectious individuals",
                     min = 0, max = 999, value = 1),

        numericInput(inputId = "lat", label = "# Latently infected individuals",
                     min = 0, max = 999, value = 0),

        numericInput(inputId = "rem", label = "# Removed individuals",
                     min = 0, max = 999, value = 0),


        shiny::sliderInput(inputId = "beta", label = HTML("Infection rate (&beta;)"),
                           min = 0.00000, max = 0.001, value = 0.00015, step = 0.00001),

        shiny::sliderInput(inputId = "omega", label = HTML("Latently infected to infectious (&omega;)"),
                           min = 0, max = 1, value = 0.1, step = 0.05),

        shiny::sliderInput(inputId = "mu", label = HTML("Infectious to removed (&mu;)"),
                           min = 0, max = 1, value = 0.1, step = 0.05),


        shiny::sliderInput(inputId = "time", label = "Simulation length",
                           min = 0, max = 1500, value = 450, step = 50),

        shiny::actionButton("refresh", label = "RUN")
      ),

      shiny::mainPanel(
        shiny::tabsetPanel(
          shiny::tabPanel(title = "Plots", plotOutput("plot1")),

          shiny::tabPanel(title = "Data",
                          shiny::DTOutput("table1"),
                          shiny::downloadButton("downloadData", "Download Results")
          )

        )
      )
    )
  )


  # server
  server <- function(input, output, session) {

    datasetInput =  shiny::eventReactive(input$refresh,{
      test1_fun(input$beta, input$omega, input$mu, input$lat, input$inf, input$rem, input$time)
    }
    )

    output$plot1 =  shiny::renderPlot({
      HLIR_plot(datasetInput())
    })

    output$table1 <- DT::renderDT({(
      datasetInput()%>%
        tidyr::pivot_wider(names_from = var, values_from = value)) %>%
        DT::datatable(rownames = FALSE,
                      options = list(searching=FALSE,pageLength = 20)) %>%
        DT::formatRound( columns = c(2:6), digits= 0)
    })


    output$downloadData <-
      downloadHandler(
        filename = "HILR_results.csv",
        content = function(file){

          write.csv(datasetInput()%>% pivot_wider(names_from = var, values_from = value), file)
        })
  }

  # Run the application
  shinyApp(ui = ui, server = server)

}

