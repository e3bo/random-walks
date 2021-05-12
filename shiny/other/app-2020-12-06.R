# Load packages
library(dplyr)
library(tidyselect)
library(tidyr)
library(readr)
library(stringr)
library(here)
library(shiny)
library(shinyWidgets)
library(ggplot2)
library(plotly)
library(RColorBrewer)
library(httr) #for reading from github

#################################
#define some variables 
#################################

#starting date for date slider and default starting date to show
mindate = as.Date("2020-02-01","%Y-%m-%d")
defaultdate = as.Date("2020-03-01","%Y-%m-%d")
#needs to follow order of scenarios
scenarionames = c("Increase social distancing", "Return to normal", "Maintain social distancing")

#################################
# functions
#################################

#function to capatalize the first letter in a string-used to make make tidy labels for the y-axis
capitalize_first <- function(str) {
  substr(str, 1, 1) <- toupper(substr(str, 1, 1))
  return(str)
}

#function to build tidy legend labels by scenario 
build_legend <- function(p_dat, scenario_var){
  #make reactive legend names based on picker inputs
  make_legend <- match(p_dat$scenario, scenario_var)
  make_legend <- recode(make_legend,
                        "1" = "Increase social distancing",
                        "2" = "Return to normal",
                        "3" = "Maintain social distancing")
  return(make_legend)
}

#################################
# Load all data
# should be online so things update automatically
# for speed, we only get data from the online scenario if the data is old, otherwise we load locally
#################################

# to ensure data gets refreshed on server, we need this
get_data <- function()
{
  filename = here("data",paste0("clean_data_",Sys.Date(),'.rds')) #if the data file for today is here, load then return from function
  if (file.exists(filename)) 
  {
     all_data <- readRDS(file = filename)    
     return(all_data)  
  }
  #if data file is not here, go through all of the below
  
  all_data = list() #will save and return all datasets as list
  
  #################################
  #pull data from repo, combine and process
  #################################
  #message('starting data processing')
  
  us_popsize <- readRDS(here("data","us_popsize.rds")) %>% rename(state_abr = state, location = state_full, pop_size = total_pop)
  us_dat_raw <- readr::read_csv("https://raw.githubusercontent.com/CEIDatUGA/COVID-stochastic-fitting/master/output/us_current_results.csv")
  
  #Notes on current data format:
  #sim_type contains 3 different future scenarios, has NA for data
  #period is Past or Future, based on date up to which data is fit. Only applies to model results
  #variable includes latent_trend and mobility_trend and combined_trend. Those are stored in mean_value. Data are stored in mean_value as well.
  us_dat <- us_dat_raw %>% 
            filter(variable != c("daily_hosps","cumulative_hosps")) %>% #not using hosp right now
            # mutate(variable = recode(variable,daily_cases = "Daily_Cases", 
            #                                   daily_deaths = "Daily_Deaths",
            #                                   daily_all_infections = "Daily_Allinfected", 
            #                                   actual_daily_deaths = "Actual_Daily_Deaths",
            #                                   actual_daily_cases = "Actual_Daily_Cases",
            #                                   cumulative_cases = "Cumulative_Cases",                         
            #                                   cumulative_deaths = "Cumulative_Deaths",
            #                                   cumulative_all_infections = "Cumulative_Allinfected", #we are missing this in the raw data
            #                                   actual_cumula_deaths = "Actual_Daily_Deaths",
            #                                   actual_daily_cases = "Actual_Daily_Cases",
            #                          
            #                                   combined_trend = "Transmissionstrength"
            #                         )) %>%
            select(-c(lower_95,lower_90,upper_95)) %>% #currently only using 80CI
            left_join(us_popsize, by = "location") %>%
            rename(populationsize = pop_size, scenario = sim_type)
  
  #combine data in list  
  #currently only US, but set up for future use
  all_data$us_dat = us_dat
  
  message('Data cleaning done.')
  
  #save the data
  saveRDS(all_data, filename)    
  return(all_data)
} # end the get-data function which pulls data from the various online scenarios and processes/saves  


###########################################
# Define server functions
###########################################
server <- function(input, output, session) 
{

  ###########################################
  # function that re-reads the data every so often
  ###########################################
  all_data <- reactivePoll(intervalMillis = 1000*60*60*3, # pull new data every N hours
                           session = NULL,
                           checkFunc = function() {Sys.time()}, #this will always return a different value, which means at intervals specified by intervalMillis the new data will be pulled
                           valueFunc = function() {
                             withProgress(message = 'Loading data, please be patient.',
                                          detail = "", value = 0,
                                          {
                                            get_data()
                                          } ) #end with-progress wrapper
                             }) #end reactive poll
  
  #read data is reactive, doesn't work for rest below 
  all_dat = isolate(all_data())
  # pull data out of list 
  us_dat = all_dat$us_dat 
  #define variables for location and scenario selectors
  state_var = sort(unique(us_dat$location))  
  scenario_var = sort(unique(us_dat$scenario))

  names(scenario_var) <- scenarionames

  #getting the date at which fitting was done, i.e. switch from 
  #past to future. Might not be same as today's date
  #done a bit clunky, could likely be done better (and tidy)
  #assigned to a global variable
  x <- us_dat %>% filter(location == "Alabama", variable == "actual_daily_cases")
  rundate <<- x$date[max(which(x$period == "Past"))]
  maxdate <<- as.Date(rundate,"%Y-%m-%d") + 42 #this works but is dependant on modeldate working and assumes the data always uses a 42 day prediction
  
  
  ################################################################################################
  #create the following UI elements on server and then add to UI since they depend on the variables above 
  #those variables are only defined on server
  ################################################################################################
  output$state_selector = renderUI({
    shinyWidgets::pickerInput("state_selector", "Select State(s)", state_var, multiple = FALSE, options = list(`actions-box` = TRUE), selected = c("Georgia") )
  })
  output$scenario_selector = renderUI({
    checkboxGroupInput("scenario_selector", label = h4("Select scenario(s)"), choices = scenario_var, selected = c("status_quo","return_normal","linear_increase_sd") )
  })
  output$x_limit = renderUI({
    sliderInput(inputId = "x_limit", "Select a range of dates to be plotted.", min = mindate,  max = maxdate, value = c(defaultdate, maxdate))
  })
  
  
  

  ###########################################
  # function that takes data and makes plots with plotly
  ###########################################
  make_plotly <- function(all_dat, location_selector, scenario_selector, daily_tot,
                          xscale, yscale, absolute_scaled, x_limit, conf_int, ylabel, outtype)
  {
    
    #create name for outcome
    outcome = paste0(daily_tot,"_",outtype)
    #daily/cumulative does not exist for combined_trend, so adjust for that
    if (outtype == 'combined_trend') {outcome = outtype}
    
    #create y-axis names
    y_tag = c(" Cases"," Deaths"," All Infected")
    ylabel = paste0(daily_tot,y_tag,sep=" ")
    ylabel = capitalize_first(ylabel)
    #apply names to plots by outtype
    if (outtype == 'cases') {ylabel = ylabel[1]}
    if (outtype == 'deaths') {ylabel = ylabel[2]}
    if (outtype == 'all_infections') {ylabel = ylabel[3]}
    
    #filter data based on user selections
    #keep all outcomes/variables for now so we can do x-axis adjustment
    #filtering of only the outcome to plot is done after x-scale adjustment
    #not filter scenario here otherwise the data gets lost
    plot_dat <- all_dat %>%   filter(location %in% location_selector) %>%      
                                      filter(date >= input$x_limit[1] &  date <= input$x_limit[2])
    
    # if we want scaling by 100K, do extra scaling 
    # this needs to be applied to all values, including CI, so need to figure out how to do-should be working now
    # don't apply to the transmissionstrength since it's an input
    if ((absolute_scaled == 'scaled') && (outtype != "combined_trend"))
    {
      plot_dat <- plot_dat %>% mutate(median_value = median_value / populationsize * 100000) %>%
                               mutate(lower_80 = lower_80 / populationsize * 100000) %>%
                               mutate(upper_80 = upper_80 / populationsize * 100000) 
    } #end scaling function
 
    linesize = 1.5
    ncols = max(3,length(unique(plot_dat$location))) #number of colors for plotting
    #define a specified colorset for all scenarios
    colorset = c('linear_increase_sd'='#5798d1','return_normal'='#a11c3e','status_quo'='#319045')
    
    # make plot
    if(outtype != "combined_trend")
    {
      p_dat <- plot_dat %>% 
        filter(variable == outcome) %>%
        filter(scenario %in% scenario_selector) %>%
        group_by(scenario,location) %>%
        arrange(date)
      
      #make reactive legend names based on picker inputs
      legend_name <- build_legend(p_dat, scenario_var)

      pl <- p_dat %>%
          plotly::plot_ly() %>% 
          plotly::add_trace(x = ~date, y = ~median_value, type = 'scatter',
                                 mode = 'lines', 
                                 linetype = ~location, 
                                 line = list(width = linesize), 
                                 name = ~legend_name, 
                                 showlegend = FALSE, #not needed to run but indicating for reference
                                 color = ~scenario, colors = colorset) %>%
                          layout(xaxis = list(title = "Date")) %>%
                          layout(yaxis = list(title=ylabel, type = yscale, size = 18)) %>%
                          layout(legend = list(orientation = "h", x = 0.2, y = -0.3))
      
      maxy = max(p_dat$median_value, na.rm = TRUE)
      
      #Adds prediction interval ribbons and/or current date bar
      if(conf_int == "Yes")
      {
        #add prediction interval ranges
        pl <- pl %>% add_ribbons(x = ~date, ymin = ~lower_80, ymax = ~upper_80, 
                                 name = "80% Prediction Interval", 
                                 color = ~scenario, showlegend = FALSE, opacity = 0.5) 
        maxy = max(p_dat$upper_80, na.rm = TRUE)
        
      }
      
      #add actual data on top of model data
      actual_data <- plot_dat %>%  
                     filter(variable == paste0("actual_",outcome)) %>%
                     group_by(scenario,location) %>%
                     arrange(date) 
      
      pl <- pl %>% plotly::add_trace(x = ~date, y = ~median_value, type = 'scatter',
                                     mode = 'lines+markers',  data = actual_data,
                                     line = list(color = "grey"), opacity = 0.5, name = "Reported Data", 
                                     showlegend = FALSE,
                                     marker = list(size = 5, color = "black", opacity = 0.5))
    } #end non-transmissions strength plots
    
    if(outtype == "combined_trend")
    {
      p_dat <- plot_dat %>% 
        filter(variable == outcome) %>%
        filter(scenario %in% scenario_selector) %>%
        group_by(scenario,location) %>%
        arrange(date)
      
      #make reactive legend names based on picker inputs
      legend_name <- build_legend(p_dat, scenario_var)

      pl <- p_dat %>%
        plotly::plot_ly() %>% 
        plotly::add_trace(x = ~date, y = ~mean_value, type = 'scatter', 
                          mode = 'lines', 
                          linetype = ~location, 
                          line = list(width = linesize), 
                          name = ~legend_name,
                          showlegend = FALSE, 
                          color = ~scenario, colors = colorset) %>%
        layout(xaxis = list(title = "Date")) %>%
        layout(yaxis = list(title="Relative Transmission Strength", type = yscale, size = 18)) %>%
        layout(legend = list(orientation = "h", x = 0.2, y = -0.3))
      
      maxy = max(p_dat$mean_value, na.rm = TRUE)
      
    }

    #add date marker
    pl <- pl %>% plotly::add_segments(x = rundate, xend = rundate, 
                                      y = 0, yend = maxy, name = rundate,
                                      color = I("black"), alpha = 0.75,
                                      showlegend = FALSE)  %>%
            layout(annotations = list(x = rundate, y = maxy, text = paste0("Last model run: ", rundate),
                                      xref = "x", yref = "y",
                                      showarrow = TRUE, arrowhead = 3,
                                      ax = 20, ay = -40))
    
    
    return(pl)
  }
  
  ###########################################
  #make plots
  output$all_plots <- renderPlotly({
    
  if (!is.null(input$scenario_selector))
  {
    #transmission strength plot
    p1 <- make_plotly(us_dat, input$state_selector, input$scenario_selector, "Daily",
                      input$xscale, input$yscale, "absolute", input$x_limit,  input$conf_int, ylabel = 1, outtype = "combined_trend")
    #create case plot
    p2 <- make_plotly(us_dat, input$state_selector, input$scenario_selector, input$daily_tot,
                      input$xscale, input$yscale, input$absolute_scaled, input$x_limit, input$conf_int, ylabel = 1, outtype = "cases")
    #death plot
    p3 <- make_plotly(us_dat, input$state_selector, input$scenario_selector, input$daily_tot,
                      input$xscale, input$yscale, input$absolute_scaled, input$x_limit, input$conf_int, ylabel = 1, outtype = "deaths")
    
    #create all  infections plot
    p4 <- make_plotly(us_dat, input$state_selector, input$scenario_selector, input$daily_tot,
                      input$xscale, input$yscale, input$absolute_scaled, input$x_limit, input$conf_int, ylabel = 1, outtype = "all_infections")
  
  all_plots <-  plotly::subplot(p1, p2, p3, p4, 
                    nrows = 4, heights = c(.25,.25,.25,.25), 
                    shareX = TRUE, titleY = TRUE
    )
  }
  }) #end function making plots
  
  
} #end server function


#################################
# Define UI
#################################
ui <- fluidPage(
  tags$head(includeHTML(here("www","google-analytics.html"))), #this is for Google analytics tracking.
  includeCSS(here("www","appstyle.css")),
  #main tabs
  sidebarLayout(
      sidebarPanel(
        uiOutput('state_selector'),
        uiOutput('scenario_selector'),
        #br(),
        #shiny::checkboxInput("conf_int", label = h4("Show 80% Prediction interval"), value = TRUE),
        shiny::radioButtons("conf_int", label = h4("Show 80% prediction interval"), choices = list("Yes" = "Yes", "No" = "No"), selected = "Yes"),
        #shiny::div("Show 80% confidence interval for forecast data."),
        shiny::radioButtons("daily_tot", label = h4("Daily or cumulative numbers"), choices = list("Daily" = "daily", "Cumulative" = "cumulative"), selected = "daily"),
       #shiny::selectInput("daily_tot", "Daily or cumulative numbers", c("Daily" = "daily", "Cumulative" = "cumulative" )),
        shiny::div("Modify all plots to show daily or cumulative data."),
        br(),
        shiny::radioButtons("absolute_scaled", label = h4("Absolute or scaled values"), choices = list("Absolute Number" = "absolute", "Per 100,000 persons" = "scaled"), selected = "absolute"),
        #shiny::selectInput( "absolute_scaled","Absolute or scaled values",c("Absolute Number" = "absolute", "Per 100,000 persons" = "scaled") ),
        shiny::div("Modify the bottom three plots to display values scaled by the state population size."),
        br(),
        uiOutput('x_limit'),
        shiny::div("Modify plots to display a specified date range."),
        br(),
        shiny::radioButtons("yscale", label = h4("Y-scale"), choices = list("Linear" = "lin", "Logarithmic" = "log"), selected = "lin"),
        #shiny::selectInput(  "yscale", "Y-scale", c("Linear" = "lin", "Logarithmic" = "log")),
        shiny::div("Modify outcome plots to show data on a linear or logarithmic scale."),
       
      ),         #end sidebar panel
      # Output:
      mainPanel(
        #change to plotOutput if using static ggplot object
        plotlyOutput(outputId = "all_plots", height = "1000px"),
        img(src='legend.png', align = "bottom", height = "100px", width = "250px"),
         tags$div(
          "These interactive plots are part of CEID @ UGA's work on COVID-19 modeling. If you ended up on this site outside our main page, see", a("the link to this website", href = "https://www.covid19.uga.edu", target = "_blank" ),
          "for more information and explanations."
        )        
      ) #end main panel
    ) #end sidebar layout     

) #end fluidpage and UI part of shiny app
#end UI of shiny app
###########################################



# Create Shiny object
shinyApp(ui = ui, server = server)

