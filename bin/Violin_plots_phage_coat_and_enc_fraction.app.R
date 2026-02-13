# app.R
# Interactive violin + jitter with hover tooltips (ggplot2 + plotly + shiny)

library(shiny)
library(tidyverse)
library(ggplot2)
library(viridis)
library(plotly)
library(rlang)

# -------------------------
# Settings
# -------------------------
default_group_order <- c("Micro", "Mini", "Control", "Fringe", "Refseq", "PICIs",
                         "ENC_1", "ENC_2", "ENC_3", "ENC_4")

special_trim_groups <- c("ENC_2", "ENC_3", "ENC_4")

# -------------------------
# Helpers
# -------------------------

build_tooltip <- function(df, cols) {
  cols <- cols[cols %in% names(df)]
  if (length(cols) == 0) return(rep("", nrow(df)))
  
  purrr::pmap_chr(df[, cols, drop = FALSE], function(...) {
    vals <- list(...)
    nm <- names(vals)
    parts <- map2_chr(
      nm, vals,
      ~ paste0("<b>", .x, ":</b> ", ifelse(is.na(.y), "NA", as.character(.y)))
    )
    paste(parts, collapse = "<br>")
  })
}

make_violin_gg <- function(df, value_col, title_text,
                           group_order,
                           positive_only,
                           show_points,
                           show_violin,
                           violin_width,
                           violin_width_factor,   # shrink max width only when scale == "width"
                           violin_adjust,
                           violin_scale,
                           violin_trim_ui,
                           violin_alpha,
                           point_size,
                           point_alpha,
                           jitter_width,
                           jitter_height,
                           y_max,
                           show_counts_in_x,
                           tooltip_cols,
                           point_color_mode) {
  
  value_sym <- if (is.character(value_col)) rlang::sym(value_col) else rlang::ensym(value_col)
  metric_name <- rlang::as_string(value_sym)
  
  if (!(metric_name %in% names(df))) stop("Column not found in dataframe: ", metric_name)
  if (!("group" %in% names(df))) stop("Missing column: group")
  
  df2 <- df %>%
    mutate(group = factor(group, levels = group_order)) %>%
    filter(group %in% group_order) %>%
    filter(!is.na(!!value_sym))
  
  if (positive_only) df2 <- df2 %>% filter(!!value_sym > 0)
  
  # x labels with counts
  if (show_counts_in_x) {
    counts <- df2 %>%
      count(group, name = "n") %>%
      tidyr::complete(group = factor(group_order, levels = group_order), fill = list(n = 0)) %>%
      arrange(group)
    
    group_labels <- setNames(
      paste0(as.character(counts$group), "\n(n=", counts$n, ")"),
      as.character(counts$group)
    )
  } else {
    group_labels <- setNames(as.character(group_order), as.character(group_order))
  }
  
  # Tooltip text for points
  df2 <- df2 %>% mutate(.tooltip = build_tooltip(., tooltip_cols))
  
  # SECRET RULE #1: always truncate violin data above 100 (density only)
  df_violin_base <- df2 %>% filter(!!value_sym <= 100)
  
  # effective width: shrink only when scale == "width"
  eff_violin_width <- violin_width * ifelse(violin_scale == "width", violin_width_factor, 1)
  
  p <- ggplot(df2, aes(x = group, y = !!value_sym))
  
  # -------------------------
  # Violin layer(s)
  # -------------------------
  if (show_violin) {
    
    # SECRET RULE #2: if metric is enc_perc_fraction, force trim=TRUE for ENC_2/ENC_3/ENC_4 only
    if (metric_name == "enc_perc_fraction") {
      
      df_violin_special <- df_violin_base %>% filter(group %in% special_trim_groups)
      df_violin_other   <- df_violin_base %>% filter(!(group %in% special_trim_groups))
      
      if (nrow(df_violin_other) > 0) {
        p <- p +
          geom_violin(
            data = df_violin_other,
            aes(fill = group),
            trim = violin_trim_ui,
            adjust = violin_adjust,
            scale = violin_scale,
            width = eff_violin_width,
            alpha = violin_alpha
          )
      }
      
      if (nrow(df_violin_special) > 0) {
        p <- p +
          geom_violin(
            data = df_violin_special,
            aes(fill = group),
            trim = TRUE,
            adjust = violin_adjust,
            scale = violin_scale,
            width = eff_violin_width,
            alpha = violin_alpha
          )
      }
      
    } else {
      p <- p +
        geom_violin(
          data = df_violin_base,
          aes(fill = group),
          trim = violin_trim_ui,
          adjust = violin_adjust,
          scale = violin_scale,
          width = eff_violin_width,
          alpha = violin_alpha
        )
    }
    
    p <- p + scale_fill_viridis(discrete = TRUE, option = "A", drop = FALSE)
  }
  
  # -------------------------
  # Points layer (color toggle)
  # -------------------------
  if (show_points) {
    
    if (identical(point_color_mode, "Group colors")) {
      p <- p +
        geom_jitter(
          aes(color = group, text = .tooltip),
          width = jitter_width,
          height = jitter_height,
          size = point_size,
          alpha = point_alpha,
          stroke = 0
        ) +
        scale_color_viridis(discrete = TRUE, option = "A", drop = FALSE)
    } else {
      p <- p +
        geom_jitter(
          aes(text = .tooltip),
          width = jitter_width,
          height = jitter_height,
          size = point_size,
          alpha = point_alpha,
          stroke = 0,
          color = "black"
        )
    }
  }
  
  p +
    scale_x_discrete(labels = group_labels, drop = FALSE) +
    coord_cartesian(ylim = c(0, y_max)) +
    theme(
      legend.position = "none",
      plot.title = element_text(size = 12),
      axis.text.x = element_text(lineheight = 0.9)
    ) +
    labs(
      title = title_text,
      x = "Group",
      y = metric_name
    )
}

# -------------------------
# UI
# -------------------------
ui <- fluidPage(
  titlePanel("Interactive Violin + Points (hover shows row info)"),
  
  # JS handler: download SVG of CURRENT plotly view (zoom/pan included)
  tags$script(HTML("
    Shiny.addCustomMessageHandler('export_svg', function(msg) {
      var id = msg.id;
      var filename = msg.filename || 'violin_plot';
      var gd = document.getElementById(id);
      if (!gd || typeof Plotly === 'undefined') {
        console.error('Plotly graph not found for id:', id);
        return;
      }
      Plotly.downloadImage(gd, {format: 'svg', filename: filename});
    });
  ")),
  
  sidebarLayout(
    sidebarPanel(
      fileInput("csv", "Load CSV (optional; if empty uses Figure3_phage_coat_and_encapsulin_fraction.csv)", accept = ".csv"),
      tags$hr(),
      
      selectInput(
        "metric",
        "Metric",
        choices = c(
          "Encapsulin (enc_perc_fraction)" = "enc_perc_fraction",
          "Clan phagecoat (percentage_phage_coat)" = "clan_phagecoat_percentage_adjusted"
        ),
        selected = "enc_perc_fraction"
      ),
      
      checkboxInput("positive_only", "Positive values only (> 0)", TRUE),
      sliderInput("y_max", "Y max", min = 1, max = 100, value = 100, step = 1),
      
      tags$hr(),
      
      checkboxInput("show_violin", "Show violin (density)", TRUE),
      
      sliderInput("violin_width", "Violin width (base)", min = 0.1, max = 1.0, value = 0.5, step = 0.05),
      sliderInput("violin_width_factor", "Max width factor (only when scale = width)",
                  min = 0.1, max = 1.00, value = 0.5, step = 0.01),
      sliderInput("violin_adjust", "Violin smoothness (adjust)", min = 0.2, max = 3.0, value = 1.0, step = 0.1),
      
      selectInput(
        "violin_scale", "Violin scale",
        choices = c("width (same max width)" = "width",
                    "area (same area)" = "area",
                    "count (scaled by n)" = "count"),
        selected = "width"
      ),
      
      checkboxInput("violin_trim", "Trim tails (trim=TRUE)", FALSE),
      sliderInput("violin_alpha", "Violin alpha", min = 0.05, max = 1.0, value = 0.6, step = 0.05),
      
      tags$hr(),
      
      checkboxInput("show_points", "Show points", TRUE),
      selectInput("point_color_mode", "Point color",
                  choices = c("Group colors", "Black"),
                  selected = "Group colors"),
      
      sliderInput("point_size", "Point size", min = 0.1, max = 5.0, value = 0.8, step = 0.1),
      sliderInput("point_alpha", "Point alpha", min = 0.05, max = 1.0, value = 0.35, step = 0.05),
      sliderInput("jitter_width", "Jitter width (x)", min = 0, max = 0.8, value = 0.15, step = 0.01),
      sliderInput("jitter_height", "Jitter height (y)", min = 0, max = 2.0, value = 0.0, step = 0.05),
      
      tags$hr(),
      checkboxInput("show_counts_in_x", "Show (n=) in x labels", TRUE),
      sliderInput("plot_height", "Plot height (px)", min = 400, max = 1100, value = 650, step = 50),
      
      tags$hr(),
      
      # NEW: Export SVG button (exports current view)
      actionButton("export_svg", "Export SVG"),
      
      tags$hr(),
      uiOutput("tooltip_cols_ui")
    ),
    
    mainPanel(
      uiOutput("plt_ui")
    )
  )
)

# -------------------------
# Server
# -------------------------
server <- function(input, output, session) {
  
  df_raw <- reactive({
    default_path <- "Figure3_phage_coat_and_encapsulin_fraction.csv"
    path <- if (!is.null(input$csv)) input$csv$datapath else default_path
    
    validate(need(file.exists(path),
                  "CSV not found. Upload it, or put Figure3_phage_coat_and_encapsulin_fraction.csv next to app.R."))
    
    read.csv(path, stringsAsFactors = FALSE)
  })
  
  df_processed <- reactive({
    df <- df_raw()
    
    validate(need("percentage_phage_coat" %in% names(df),
                  "Missing column: percentage_phage_coat"))
    
    df %>%
      mutate(
        clan_phagecoat_percentage = percentage_phage_coat,
        clan_phagecoat_percentage_adjusted = percentage_phage_coat
      )
  })
  
  output$tooltip_cols_ui <- renderUI({
    df <- df_processed()
    default_cols <- intersect(
      c("protein_id", "group", "enc_perc_fraction", "percentage_phage_coat", "clan_phagecoat_percentage_adjusted"),
      names(df)
    )
    
    selectizeInput(
      "tooltip_cols",
      "Columns to show on hover (points)",
      choices = names(df),
      selected = default_cols,
      multiple = TRUE,
      options = list(plugins = list("remove_button"), maxItems = 25)
    )
  })
  
  output$plt_ui <- renderUI({
    plotlyOutput("plt", height = paste0(input$plot_height, "px"))
  })
  
  output$plt <- renderPlotly({
    df <- df_processed()
    
    metric <- input$metric
    title_text <- if (metric == "enc_perc_fraction") {
      "Encapsulin Percentage by Group"
    } else {
      "Clan Phagecoat Percentage by Group"
    }
    
    tooltip_cols <- if (is.null(input$tooltip_cols)) character(0) else input$tooltip_cols
    
    p <- make_violin_gg(
      df = df,
      value_col = metric,
      title_text = title_text,
      group_order = default_group_order,
      positive_only = input$positive_only,
      show_points = input$show_points,
      show_violin = input$show_violin,
      violin_width = input$violin_width,
      violin_width_factor = input$violin_width_factor,
      violin_adjust = input$violin_adjust,
      violin_scale = input$violin_scale,
      violin_trim_ui = input$violin_trim,
      violin_alpha = input$violin_alpha,
      point_size = input$point_size,
      point_alpha = input$point_alpha,
      jitter_width = input$jitter_width,
      jitter_height = input$jitter_height,
      y_max = input$y_max,
      show_counts_in_x = input$show_counts_in_x,
      tooltip_cols = tooltip_cols,
      point_color_mode = input$point_color_mode
    )
    
    ggplotly(p, tooltip = "text") %>%
      layout(
        hovermode = "closest",
        margin = list(l = 70, r = 20, b = 90, t = 60)
      ) %>%
      config(
        displayModeBar = TRUE
      )
  })
  
  # NEW: Export SVG of current view (client-side)
  observeEvent(input$export_svg, {
    ts <- format(Sys.time(), "%Y%m%d_%H%M%S")
    fname <- paste0("violin_", input$metric, "_", ts)
    session$sendCustomMessage("export_svg", list(id = "plt", filename = fname))
  })
}

shinyApp(ui, server)

# END
