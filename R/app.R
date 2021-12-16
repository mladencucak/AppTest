library("ggplot2")
library("shinyscreenshot")
library("shiny")
library("deSolve")

###########################################################
### Logistic model app
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
  ui <- fluidPage(

    # titlePanel("Logistic model"), #add title

    # Set sliders to be on the side
    sidebarLayout(
      sidebarPanel(
        # Add slider for initial disease and set parameters
        sliderInput("initial", slider_disease,
                    min = min_slider_yo, max = max_slider_yo, value = default_yo
        ),
        sliderInput("rate", slider_rate,
                    min = min_slider_r, max = max_slider_r, value = default_r
        ),

        radioButtons("plot_choice",
                     label = "Plot",
                     choices = c("Together", "Separate"),
                     selected = c("Together"),
                     inline = TRUE),
        actionButton("go", "Take a screenshot")


      ),
      # Main panel for displaying outputs ----
      mainPanel(
        # Output:
        plotOutput("myplot")
      )
    )
  )

  # Define server logic for slider examples
  server <-
    function(input, output, session) {
      output$myplot <-
        renderPlot({
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
            ggplot(dta) +
            geom_line(aes(time, dis, color = model)) +
            ylim(0, 1) +
            theme_bw() +
            scale_y_continuous(
              limits = c(0, 1),
              expand = c(0, 0),
              breaks = seq(0, 1, 0.2),
              name = "Disease"
            ) +
            scale_x_continuous(
              expand = c(0, 0),
              breaks = seq(0, 80, 20),
              name = "Time"
            ) +
            theme(legend.position = "bottom",
                  text = element_text(size = 12))

          # plot_dis
          if (input$plot_choice == "Separate") {
            plot_dis +
              facet_wrap( ~ model) +
              theme(legend.position = "none")
          }

          if (input$plot_choice == "Together") {
            plot_dis
          }

        })
      observeEvent(input$go, {
        screenshot(
          # selector="#myplot",
          filename = "combined"
        )
      })


    }
  # Create Shiny app
  shinyApp(ui = ui, server = server)

}



