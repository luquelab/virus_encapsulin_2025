# Interactive exploration of Fig. 3D and smoothing curves (USER-UPLOADED INPUT)

library(shiny)
library(dplyr)
library(ggplot2)
library(mgcv)
library(scales)
library(hexbin)

# -------------------------
# UI
# -------------------------
ui <- fluidPage(
  titlePanel("Encapsulin fraction vs. genome length"),
  sidebarLayout(
    sidebarPanel(
      strong("1) Upload your data"),
      fileInput(
        "datafile",
        "Upload CSV file",
        accept = c(".csv")
      ),
      checkboxInput("header", "File has header row", value = TRUE),
      radioButtons(
        "sep",
        "Separator",
        choices = c("Comma (,)" = ",", "Semicolon (;)" = ";", "Tab (\\t)" = "\t"),
        selected = ",",
        inline = TRUE
      ),
      tags$small(
        "Expected: a column for encapsulin fraction (0–1), a column for genome length in kbp,",
        "and optionally a categorical group column."
      ),
      hr(),
      
      strong("2) Map your columns"),
      uiOutput("colSelectors"),
      hr(),
      
      strong("3) Smoother settings"),
      selectInput(
        "smoothMethod",
        "Smoother:",
        choices = c("GAM (mgcv)" = "gam",
                    "LOESS"      = "loess",
                    "Smoothing spline" = "spline"),
        selected = "spline"
      ),
      sliderInput(
        "k", "Spline basis dimension k [GAM]:",
        min = 3, max = 80, value = 10, step = 1
      ),
      sliderInput(
        "loessSpan", "LOESS span [LOESS]:",
        min = 0.1, max = 1.0, value = 0.6, step = 0.05
      ),
      sliderInput(
        "spar", "Spline spar [Smoothing spline]:",
        min = 0.1, max = 1.0, value = 0.5, step = 0.05
      ),
      sliderInput(
        "bins", "Number of hexbin bins:",
        min = 10, max = 60, value = 25, step = 1
      ),
      sliderInput(
        "ciLevel", "CI level for GAM (0 = no ribbon):",
        min = 0, max = 0.95, value = 0.50, step = 0.05
      ),
      checkboxInput(
        "signalOnly",
        "Use only genomes with encapsulin signal (> 0) for fit",
        value = FALSE
      ),
      checkboxInput(
        "showZeroPoints",
        "Show points with encapsulin_fraction = 0 in the scatter",
        value = TRUE
      ),
      checkboxInput(
        "fitLogX",
        "Fit on log10(genome length) instead of raw (x-axis stays in kbp)",
        value = FALSE
      ),
      checkboxInput(
        "clampZero",
        "Clamp curve at 0 (no negative encapsulin fraction)",
        value = TRUE
      ),
      hr(),
      strong("Current model specification:"),
      verbatimTextOutput("modelText", placeholder = TRUE)
    ),
    
    mainPanel(
      plotOutput("fig3D", height = "600px"),
      hr(),
      strong("Loaded data preview:"),
      tableOutput("preview")
    )
  )
)

# -------------------------
# Server
# -------------------------
server <- function(input, output, session) {
  
  # 1) Read raw uploaded data
  raw_data <- reactive({
    req(input$datafile)
    
    tryCatch(
      read.csv(
        input$datafile$datapath,
        header = isTRUE(input$header),
        sep = input$sep,
        stringsAsFactors = FALSE,
        check.names = FALSE
      ),
      error = function(e) {
        showNotification(paste("Error reading CSV:", e$message), type = "error")
        NULL
      }
    )
  })
  
  # 2) Column selectors (appear after upload)
  output$colSelectors <- renderUI({
    df <- raw_data()
    if (is.null(df)) return(tags$em("Upload a CSV to select columns."))
    
    cols <- names(df)
    
    # Try smart defaults if typical names exist
    default_enc <- if ("encapsulin_fraction" %in% cols) "encapsulin_fraction" else cols[1]
    default_len <- if ("query_genome_length_kbp" %in% cols) "query_genome_length_kbp" else cols[min(2, length(cols))]
    default_grp <- if ("group" %in% cols) "group" else "(none)"
    
    tagList(
      selectInput("col_enc", "Encapsulin fraction column:", choices = cols, selected = default_enc),
      selectInput("col_len", "Genome length (kbp) column:", choices = cols, selected = default_len),
      selectInput("col_grp", "Group column (optional):", choices = c("(none)", cols), selected = default_grp)
    )
  })
  
  # 3) Standardize to the expected internal column names
  data_all <- reactive({
    df <- raw_data()
    req(df, input$col_enc, input$col_len, input$col_grp)
    
    validate(
      need(input$col_enc %in% names(df), "Encapsulin fraction column not found."),
      need(input$col_len %in% names(df), "Genome length column not found.")
    )
    
    out <- df %>%
      transmute(
        encapsulin_fraction      = suppressWarnings(as.numeric(.data[[input$col_enc]])),
        query_genome_length_kbp  = suppressWarnings(as.numeric(.data[[input$col_len]])),
        group = if (!is.null(input$col_grp) &&
                    input$col_grp != "(none)" &&
                    input$col_grp %in% names(df)) {
          as.factor(.data[[input$col_grp]])
        } else {
          factor("All")
        }
      )
    
    validate(
      need(any(!is.na(out$encapsulin_fraction)), "Encapsulin fraction column has no numeric values (all NA after conversion)."),
      need(any(!is.na(out$query_genome_length_kbp)), "Genome length column has no numeric values (all NA after conversion).")
    )
    
    out
  })
  
  # Preview table
  output$preview <- renderTable({
    df <- data_all()
    req(df)
    head(df, 10)
  })
  
  # -------------------------
  # Existing logic, now using data_all()
  # -------------------------
  data_signal <- reactive({
    data_all() %>%
      dplyr::filter(!is.na(encapsulin_fraction),
                    !is.na(query_genome_length_kbp),
                    encapsulin_fraction > 0)
  })
  
  data_scatter <- reactive({
    df <- data_all()
    req(df)
    
    if (input$signalOnly) {
      data_signal()
    } else {
      if (input$showZeroPoints) {
        df %>%
          dplyr::filter(!is.na(encapsulin_fraction),
                        !is.na(query_genome_length_kbp))
      } else {
        data_signal()
      }
    }
  })
  
  predictions <- reactive({
    df_fit <- data_scatter() %>%
      dplyr::filter(!is.na(encapsulin_fraction),
                    !is.na(query_genome_length_kbp))
    
    validate(need(nrow(df_fit) >= 5, "Not enough points to fit a curve (need at least ~5)."))
    
    rng <- range(df_fit$query_genome_length_kbp, na.rm = TRUE)
    x_grid <- seq(rng[1], rng[2], length.out = 200)
    
    method <- input$smoothMethod
    clamp  <- isTRUE(input$clampZero)
    
    if (method == "gam") {
      k_val <- input$k
      
      if (input$fitLogX) {
        df_fit2 <- df_fit %>% dplyr::mutate(x_var = log10(query_genome_length_kbp))
        newdata <- data.frame(x_var = log10(x_grid))
      } else {
        df_fit2 <- df_fit %>% dplyr::mutate(x_var = query_genome_length_kbp)
        newdata <- data.frame(x_var = x_grid)
      }
      
      gam_fit <- mgcv::gam(
        encapsulin_fraction ~ s(x_var, k = k_val),
        data   = df_fit2,
        family = gaussian()
      )
      
      pred <- predict(gam_fit, newdata = newdata, se.fit = TRUE)
      y  <- as.numeric(pred$fit)
      se <- as.numeric(pred$se.fit)
      
      level <- input$ciLevel
      if (level > 0 && level < 1) {
        z <- qnorm((1 + level) / 2)
        lower <- y - z * se
        upper <- y + z * se
      } else {
        lower <- upper <- rep(NA_real_, length(y))
      }
      
      if (clamp) {
        y     <- pmax(y, 0)
        lower <- pmax(lower, 0)
        upper <- pmax(upper, 0)
      }
      
      data.frame(
        method  = "gam",
        x       = x_grid,
        y       = y,
        y_lower = lower,
        y_upper = upper
      )
      
    } else if (method == "loess") {
      if (input$fitLogX) {
        df_fit2 <- df_fit %>% dplyr::mutate(x_var = log10(query_genome_length_kbp))
        newx <- log10(x_grid)
      } else {
        df_fit2 <- df_fit %>% dplyr::mutate(x_var = query_genome_length_kbp)
        newx <- x_grid
      }
      
      loess_fit <- loess(
        encapsulin_fraction ~ x_var,
        data    = df_fit2,
        span    = input$loessSpan,
        degree  = 2,
        surface = "direct"
      )
      
      y <- as.numeric(predict(loess_fit, newdata = data.frame(x_var = newx)))
      if (clamp) y <- pmax(y, 0)
      
      data.frame(
        method  = "loess",
        x       = x_grid,
        y       = y,
        y_lower = NA_real_,
        y_upper = NA_real_
      )
      
    } else {
      if (input$fitLogX) {
        x_var <- log10(df_fit$query_genome_length_kbp)
        newx  <- log10(x_grid)
      } else {
        x_var <- df_fit$query_genome_length_kbp
        newx  <- x_grid
      }
      
      ss_fit <- smooth.spline(
        x    = x_var,
        y    = df_fit$encapsulin_fraction,
        spar = input$spar
      )
      
      y <- as.numeric(predict(ss_fit, x = newx)$y)
      if (clamp) y <- pmax(y, 0)
      
      data.frame(
        method  = "spline",
        x       = x_grid,
        y       = y,
        y_lower = NA_real_,
        y_upper = NA_real_
      )
    }
  })
  
  output$fig3D <- renderPlot({
    df_hex <- data_all()
    req(df_hex)
    
    df_scatter <- data_scatter()
    pred_df    <- predictions()
    method     <- unique(pred_df$method)
    
    p <- ggplot() +
      stat_binhex(
        data = df_hex,
        aes(
          x = query_genome_length_kbp,
          y = encapsulin_fraction,
          fill = pmin(..count.., 5)
        ),
        bins = input$bins,
        na.rm = TRUE
      ) +
      scale_fill_gradient(
        name   = "Bin count",
        low    = "grey85",
        high   = "grey20",
        limits = c(1, 5),
        breaks = c(1, 5),
        labels = c("1", "5+")
      ) +
      geom_point(
        data  = df_scatter,
        aes(x = query_genome_length_kbp, y = encapsulin_fraction, shape = group),
        size  = 1.5,
        stroke= 0.3,
        colour= "black",
        fill  = "white"
      )
    
    if (method == "gam" && any(!is.na(pred_df$y_lower))) {
      p <- p + geom_ribbon(
        data = pred_df,
        aes(x = x, ymin = y_lower, ymax = y_upper),
        inherit.aes = FALSE,
        fill  = "grey80",
        alpha = 0.5
      )
    }
    
    p +
      geom_line(
        data = pred_df,
        aes(x = x, y = y),
        inherit.aes = FALSE,
        colour = "black",
        linewidth = 0.6
      ) +
      scale_x_log10(
        name   = "Genome length (kbp)",
        breaks = c(1, 5, 10, 50, 100, 200, 300),
        labels = scales::comma
      ) +
      coord_cartesian(ylim = c(0, 1)) +
      ylab("Encapsulin fraction") +
      theme_minimal() +
      theme(panel.background = element_rect(fill = "white", colour = NA)) +
      guides(shape = guide_legend(title = "Group")) +
      ggtitle("Encapsulin fraction vs. genome length")
  })
  
  output$modelText <- renderText({
    method <- input$smoothMethod
    
    base_txt <- paste0(
      "Smoother: ",
      if (method == "gam") "GAM (mgcv)"
      else if (method == "loess") "LOESS"
      else "Cubic smoothing spline",
      "\n"
    )
    
    model_txt <- if (method == "gam") {
      paste0(
        "Model:\n  encapsulin_fraction ~ s(",
        if (input$fitLogX) "log10(query_genome_length_kbp)" else "query_genome_length_kbp",
        ", k = ", input$k, ")\n"
      )
    } else if (method == "loess") {
      paste0(
        "Model:\n  encapsulin_fraction ~ loess(",
        if (input$fitLogX) "log10(query_genome_length_kbp)" else "query_genome_length_kbp",
        ", span = ", input$loessSpan, ")\n"
      )
    } else {
      paste0(
        "Model:\n  encapsulin_fraction ~ smooth.spline(",
        if (input$fitLogX) "log10(query_genome_length_kbp)" else "query_genome_length_kbp",
        ", spar = ", input$spar, ")\n"
      )
    }
    
    data_txt <- if (input$signalOnly) {
      "Data used for fit:\n  - Only genomes with encapsulin_fraction > 0\n"
    } else if (!input$signalOnly && input$showZeroPoints) {
      "Data used for fit:\n  - All genomes (including encapsulin_fraction = 0)\n"
    } else {
      "Data used for fit:\n  - All genomes for hexbin, but only encapsulin_fraction > 0 for points & curve\n"
    }
    
    other_txt <- paste0(
      "Other parameters:\n",
      "  - Hexbin bins       = ", input$bins, "\n",
      if (method == "gam") paste0("  - CI level (GAM)    = ", input$ciLevel, "\n") else "",
      "  - Clamp at 0?       = ", if (input$clampZero) "YES" else "NO", "\n"
    )
    
    paste0(base_txt, model_txt, "\n", data_txt, "\n", other_txt)
  })
}

shinyApp(ui, server)
# END