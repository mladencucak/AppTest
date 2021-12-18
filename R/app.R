

###########################################################
### Growth models app
###########################################################
GrowModAPP <- function(){

  slider_disease <- "Initial disease(y0):"
  slider_rate <- "The rate of disease progress(r):"

  min_slider_yo <-  1e-10
  max_slider_yo <- .02
  default_yo <- 1e-3

  min_slider_r <-  1e-10
  max_slider_r <- .02
  default_r <- .005

  # Define UI for slider demo app
  ui <- shiny::fluidPage(

    # titlePanel("Logistic model"), #add title

    # Set sliders to be on the side
    shiny::sidebarLayout(
      shiny::sidebarPanel(
        # Add slider for initial disease and set parameters
        shiny::sliderInput("initial", slider_disease,
                    min = min_slider_yo, max = max_slider_yo, value = default_yo
        ),
        shiny::sliderInput("rate", slider_rate,
                    min = min_slider_r, max = max_slider_r, value = default_r
        ),

        shiny::radioButtons("plot_choice",
                     label = "Plot",
                     choices = c("Together", "Separate"),
                     selected = c("Together"),
                     inline = TRUE),
        shiny::actionButton("go", "Take a screenshot")


      ),
      # Main panel for displaying outputs ----
      shiny::mainPanel(
        # Output:
        shiny::plotOutput("myplot")
      )
    )
  )

  # Define server logic for slider examples
  server <-
    function(input, output, session) {
      output$myplot <-
        shiny::renderPlot({
          r = input$rate
          time <-  seq(0, 80, by = 2)

          y <- ex <- gom <- mon <- numeric(length(time))
          y[1] <- ex[1] <- gom[1] <-  mon[1] <- input$initial



          logi_fun <- function (t, y, par) {
            y <- y[1]
            r <- par$r
            dy <- y * r * (1 - y)
            return(list(c(dy)))
          }

          mono_fun <- function (t, y, par) {
            y <- y[1]
            r <- par$r
            dy <- r * (1 - y)
            return(list(c(dy)))
          }

          expo_fun <- function (t, y, par) {
            y <- y[1]
            r <- par$r
            dy <- y * r
            return(list(c(dy)))
          }

          gompi_fun <- function (t, y, par) {
            y <- y[1]
            r <- par$r
            dy <- y * r * (log(1) - log(y))
            return(list(c(dy)))
          }

          for (k in 1:(length(time) - 1)) {
            r[k + 1] <- r[k]
            InitCond <- c(y[k])
            parms <- list(r = r[k])
            # logistic
            ode_logi <- deSolve::ode(c(y[k]), time, logi_fun, parms)
            y[k + 1] <- ode_logi[length(ode_logi[, 2]), 2]
            # exp
            ode_ex <- deSolve::ode(c(ex[k]), time, expo_fun, parms)
            ex[k + 1] <- ode_ex[length(ode_ex[, 2]), 2]
            # mon
            ode_mon <- deSolve::ode(c(mon[k]), time, mono_fun, parms)
            mon[k + 1] <- ode_mon[length(ode_mon[, 2]), 2]

            #gomp
            ode_g <- deSolve::ode(c(gom[k]), time, gompi_fun, parms)
            gom[k + 1] <- ode_g[length(ode_g[, 2]), 2]
          }

          # data.frame(time = time, log = y, exp = ex, gom = gom, mon = mon)
          dta <-
            data.frame(
              model = rep(
                c("Monomolecular", "Exponential", "Logistic", "Gomperz"),
                each = length(time)
              ),
              time = rep(time, 4),
              dis = c(mon, ex, y, gom)
            )

          # Fix exponential curve to go all the way up to 1
          if (!is.na(dta[which(dta$model == "Exponential" &
                               dta$dis > .9999)[1], "dis"])) {
            dta[which(dta$model == "Exponential" &
                        dta$dis > .9999)[1], "dis"] <- .99999
          }
          # dta <-  dta[- which(dta$model== "Exponential"&dta$dis > .9999999),]

          plot_dis <-
            ggplot2::ggplot(dta) +
            ggplot2::geom_line(ggplot2::aes(time, dis, color = model)) +
            ggplot2::ylim(0, 1) +
            ggplot2::theme_bw() +
            ggplot2::scale_y_continuous(
              limits = c(0, 1),
              expand = c(0, 0),
              breaks = seq(0, 1, 0.2),
              name = "Disease"
            ) +
            ggplot2::scale_x_continuous(
              expand = c(0, 0),
              breaks = seq(0, 80, 20),
              name = "Time"
            ) +
            ggplot2::theme(legend.position = "bottom",
                  text = ggplot2::element_text(size = 12))

          # plot_dis
          if (input$plot_choice == "Separate") {
            plot_dis +
              ggplot2::facet_wrap( ~ model) +
              ggplot2::theme(legend.position = "none")
          }

          if (input$plot_choice == "Together") {
            plot_dis
          }

        })
      shiny::observeEvent(input$go, {
        shinyscreenshot::screenshot(
          # selector="#myplot",
          filename = "combined"
        )
      })


    }
  # Create Shiny app
  shiny::shinyApp(ui = ui, server = server)

}




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

