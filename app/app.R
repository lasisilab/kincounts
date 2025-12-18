#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    https://shiny.posit.co/
#

# Load required libraries
library(shiny)        # Web application framework
library(bslib)        # Bootstrap themes and components
library(pscl)         # Zero-inflated models
library(DT)           # Interactive data tables
library(dplyr)        # Data manipulation
library(tidyr)        # Data reshaping
library(viridisLite)  # Color palettes
library(ggplot2)      # Data visualization

# Load utility functions for fertility estimation, simulation, and probability calculations
source("utils.R")

# Define UI for application (Tabset with Tab Nav)
ui <- page_fluid(
    "Long range familial search",

    # Tab panels for different sections of the application
    navset_tab(
      # Introduction tab: provides overview and documentation
      nav_panel(
        title = "Introduction",
        card(
          card_header("Long Range Familial Search: Estimating Identity Probabilities"),
          card_body(
            h3("Overview"),
            p("This application helps estimate the probability of identifying individuals through genetic databases. 
              It estimates the fertility and sibling parameters based on fertility count data, 
              simulates the distribution of distant relatives (cousins),
              and calculates the likelihood of finding genetic matches in a database."),
            
            h3("How It Works"),
            tags$ol(
              tags$li(strong("Fertility Distribution:"), " Upload fertility count data to estimate parameters for the 
                      zero-inflated negative binomial distribution that models fertility and sibling patterns."),
              tags$li(strong("Simulation:"), " Generate simulated relative counts based on fertility parameters, 
                      including first, second, and third cousins."),
              tags$li(strong("Identity Probability:"), " Calculate the probability of finding at least one genetic match 
                      in a database of a given size.")
            ),
            
            h3("Getting Started"),
            p("Navigate through the tabs to:"),
            tags$ul(
              tags$li(strong("Fertility distribution:"), " Estimate fertility parameters from your data"),
              tags$li(strong("Simulation:"), " Adjust parameters and simulate relative counts"),
              tags$li(strong("Identity probability:"), " Upload simulated data and calculate match probabilities")
            ),
            
            h3("Key Parameters"),
            tags$ul(
              tags$li(strong("Pi (π):"), " Zero-inflation probability - the excess proportion of individuals with no children"),
              tags$li(strong("Size:"), " Negative binomial dispersion parameter"),
              tags$li(strong("Mu (μ):"), " Negative binomial mean parameter"),
              tags$li(strong("Population size:"), " Total population being considered"),
              tags$li(strong("Database size:"), " Number of individuals in the genetic database")
            ),
            
            h3("Methods"),
            p("The application uses a zero-inflated negative binomial (ZINB) model to capture fertility distributions. 
              Genetic match probabilities are calculated for cousins of varying degrees.")
          )
        )
      ),

      # Fertility distribution tab: upload data and estimate parameters
      nav_panel(
        title = "Fertility distribution",
        sidebarLayout(
          sidebarPanel(
            # File input for fertility count data (must contain 'num_children' column)
            fileInput("fert_count", "Upload fertility count data (CSV)", accept = ".csv")
          ),
          mainPanel(
            # Display estimated parameters for fertility and sibling distributions
            h4("Estimated Fertility and Sibling Distribution Parameters"),
            dataTableOutput("fert_params")
          )
        )
      ),

      # Simulation tab: configure parameters and generate simulated relative counts
      nav_panel(
        title = "Simulation",
        sidebarLayout(
          sidebarPanel(
            # Zero-inflated negative binomial parameters
            h4("Fertility distribution parameters:"),
            numericInput("pi", "Pi (Zero-inflation probability)", value = 0.1, min = 0, max = 1, step = 0.05),
            numericInput("size", "Size (Negative Binomial size)", value = 3, min = 0),
            numericInput("mu", "Mu (Negative Binomial mean)", value = 3, min = 0),
            # Simulation configuration
            h4("Simulation parameters:"),
            numericInput("n_sim", "Number of simulations", value = 10000, min = 1),
            numericInput("max_kids", "Maximum number of children to consider", value = 20, min = 1),
            # Download button for simulated data
            h4("Download simulated relative counts:"),
            downloadButton("download_relative_counts", "Download")
          ),
          mainPanel(
            # Preview of simulated data
            h4("Head of simulated relative counts:"),
            dataTableOutput("relative_counts_head"),
            # Summary statistics for all relationship types
            h4("Summary of simulated relative counts:"),
            dataTableOutput("relative_counts_summary"),
            # Visualizations of relative count distributions
            h4("Boxplots of simulated relative counts:"),
            plotOutput("relative_counts_boxplot_1"),  # First through second degree + first cousins
            plotOutput("relative_counts_boxplot_2"),  # First through third cousins comparison
            # Probability mass functions for children and siblings
            h4("PMF for children, and siblings"),
            plotOutput("pmf_plot_children"),
            plotOutput("pmf_plot_siblings")
          )
        )
      ),

      # Identity probability tab: calculate match probabilities given population and database sizes
      nav_panel(
        title = "Identity probability",
        sidebarLayout(
          sidebarPanel(
            # Upload simulated relative counts from the Simulation tab
            fileInput("relative_counts", "Upload simulated relative count data (CSV)", accept = ".csv"),
            # Population and database parameters
            numericInput("pop_size", "Population size", value = 3e8, min = 0),
            numericInput("db_size", "Database size", value = 1e7, min = 0)
          ),
          mainPanel(
            # Display probability of finding at least one match for each cousin degree
            h4("Estimated probability of finding at least one cousin:"),
            verbatimTextOutput("id_prob")
          )
        )
      )
    ),

    id = "tabs"
)

# Define server logic
server <- function(input, output) {
  # ========== Fertility Distribution Tab ==========
  # Estimate fertility and sibling distribution parameters from uploaded data
  output$fert_params <- renderDataTable({
    req(input$fert_count)  # Require file upload
    fert_count <- read.csv(input$fert_count$datapath)
    # Validate that required column exists
    validate(need("num_children" %in% colnames(fert_count),
                  "The uploaded CSV must contain a column named 'num_children'."))
    # Estimate parameters using zero-inflated negative binomial model
    fert_parameters <- estimate_fertility_parameters(fert_count)
    # Display results in interactive table
    datatable(fert_parameters |> mutate(across(where(is.numeric), ~ round(.x, digits = 3))),
              colnames = c("Distribution", "Pi", "Size", "Prob", "Mu", "Mean", "Variance"),
              caption = "Estimated Fertility and Sibling Distribution Parameters")
  })

  # ========== Simulation Tab ==========
  # Reactive expression to generate simulated relative counts based on user parameters
  relative_counts <- reactive({
    # Extract parameters from UI inputs
    pi = input$pi; size = input$size; mu = input$mu; n_sim = input$n_sim; max_kids = input$max_kids
    # Run simulation using ZINB distribution
    simulate_relatives_zinb(mu, size, pi, n_sim, max_kids)
  })

  # Download handler for simulated relative counts
  output$download_relative_counts <- downloadHandler(
    filename = function() {
      paste("simulated_relative_counts.csv")
    },
    content = function(file) {
      # Write the reactive data to CSV file
      write.csv(relative_counts(), file, row.names = FALSE)
    }
  )

  # Display first few rows of simulated data
  output$relative_counts_head <- renderDataTable({
    datatable(head(relative_counts()),
              caption = "Head of Simulated Relative Counts")
  })

  # Calculate and display summary statistics for all relationship types
  output$relative_counts_summary <- renderDataTable({
    rel_counts <- relative_counts()
    # Transform data: remove 'num_' prefix, pivot to long format, and calculate statistics
    rel_counts_summary <- rel_counts %>% rename_with(~ sub("^num_?", "", .x)) %>%
      {
        col_order <- names(.)                                # Save current column order
        pivot_longer(., everything(),
                     names_to = "relationship", values_to = "value") %>%
          mutate(relationship = factor(relationship, levels = col_order))  # Preserve order
      } %>%
      group_by(relationship) %>%
      summarise(
        mean   = mean(value, na.rm = TRUE),
        sd     = sd(value, na.rm = TRUE),
        median = median(value, na.rm = TRUE),
        max = max(value, na.rm = TRUE),
        .groups = "drop"
      )

    # Display summary table with rounded values
    datatable(rel_counts_summary |> mutate(across(where(is.numeric), ~ round(.x, digits = 2))),
              caption = "Summary Statistics of Simulated Relative Counts")
  })

  # Reactive expression to transform data to long format for plotting
  relative_counts_long <- reactive({
    rel_counts <- relative_counts()
    # Transform to long format while preserving relationship order
    rel_counts %>%
      rename_with(~ sub("^num_?", "", .x)) %>%
      {
        col_order <- names(.)                                # Save current column order
        pivot_longer(., everything(),
                     names_to = "relationship", values_to = "count") %>%
          mutate(relationship = factor(relationship, levels = col_order))  # Preserve order
      }
  })

  # Boxplot for close relationships (first and second degree + first cousins)
  output$relative_counts_boxplot_1 <- renderPlot({
    rel_counts_long <- relative_counts_long()
    # Filter for specified relationship types
    ggplot(rel_counts_long %>% filter(relationship %in% c("parents", "siblings", "children", "first_degree", "grandparents", "aunts_uncles", "nephews_nieces", "second_degree", "first_cousins")),
           aes(x = relationship, y = count, fill = relationship)) +
      geom_boxplot() +
      labs(title = "Boxplot of Simulated Relative Counts",
           x = "Relationship",
           y = "Relative count",
           fill = "Relationship") +
      theme_minimal() +
      theme(legend.position = "top", axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
  })

  # Boxplot comparing cousin degrees (first through third cousins)
  output$relative_counts_boxplot_2 <- renderPlot({
    rel_counts_long <- relative_counts_long()
    # Focus on cousin relationships for comparison
    ggplot(rel_counts_long %>% filter(relationship %in% c("first_cousins", "second_cousins", "third_cousins")),
           aes(x = relationship, y = count, fill = relationship)) +
      geom_boxplot() +
      labs(title = "Boxplot of Simulated Relative Counts",
           x = "Relationship",
           y = "Relative count",
           fill = "Relationship") +
      theme_minimal() +
      theme(legend.position = "top", axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
  })

  # Probability mass function plot for children
  output$pmf_plot_children <- renderPlot({
    rel_counts <- relative_counts()
    plot_pmf(rel_counts, count_var = "num_children", category_label = "Children")
  })

  # Probability mass function plot for siblings
  output$pmf_plot_siblings <- renderPlot({
    rel_counts <- relative_counts()
    plot_pmf(rel_counts, count_var = "num_siblings", category_label = "Siblings")
  })

  # ========== Identity Probability Tab ==========
  # Calculate probability of finding at least one cousin match in the database
  output$id_prob <- renderText({
    req(input$relative_counts)  # Require file upload
    relative_counts <- read.csv(input$relative_counts$datapath)
    # Validate that required cousin count columns exist
    validate(need(all(c("num_first_cousins", "num_second_cousins", "num_third_cousins") %in% colnames(relative_counts)),
                  "The uploaded CSV must contain columns named 'num_first_cousins', 'num_second_cousins', and 'num_third_cousins'."))
    # Get population and database parameters
    pop_size <- input$pop_size
    db_size <- input$db_size

    # Calculate identity probabilities for each cousin degree
    id_prob <- identity_probability(relative_counts, pop_size, db_size)
    id_prob <- round(id_prob, digits = 4)

    # Format output as text
    paste0("first cousin: ", id_prob["first_cousin"],
           "\nsecond cousin: ", id_prob["second_cousin"],
           "\nthird cousin: ", id_prob["third_cousin"])
  })

}

# Run the application
shinyApp(ui = ui, server = server)
