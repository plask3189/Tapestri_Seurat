# =============================================================================
# Seurat NGT Layer Explorer — Interactive Shiny App
# Visualizes AF, DP, and GQ layers from seurat_obj@assays[["NGT"]]
# =============================================================================
# install.packages(c("shiny", "bslib", "ggplot2", "plotly",
#                    "dplyr", "tidyr", "DT", "shinyWidgets", "viridis", "scales"))
# =============================================================================

library(shiny)
library(bslib)
library(ggplot2)
library(plotly)
library(dplyr)
library(tidyr)
library(DT)
library(shinyWidgets)
library(viridis)

# ── Helper: safely extract a layer ───────────────────────────────────────────
get_layer <- function(seurat_obj, layer_name) {
  tryCatch({
    mat <- seurat_obj@assays[["NGT"]]@layers[[layer_name]]
    if (is.null(mat)) stop("Layer not found")
    as.matrix(mat)
  }, error = function(e) NULL)
}

# ── DEMO DATA ─────────────────────────────────────────────────────────────────
make_demo_data <- function(n_cells = 300, n_variants = 50) {
  set.seed(42)
  list(
    AF = matrix(rbeta(n_cells * n_variants, 0.5, 2),
                nrow = n_cells, ncol = n_variants,
                dimnames = list(paste0("Cell_", seq_len(n_cells)),
                                paste0("Var_",  seq_len(n_variants)))),
    DP = matrix(rnbinom(n_cells * n_variants, mu = 30, size = 5),
                nrow = n_cells, ncol = n_variants,
                dimnames = list(paste0("Cell_", seq_len(n_cells)),
                                paste0("Var_",  seq_len(n_variants)))),
    GQ = matrix(pmin(99, rnorm(n_cells * n_variants, mean = 60, sd = 20)),
                nrow = n_cells, ncol = n_variants,
                dimnames = list(paste0("Cell_", seq_len(n_cells)),
                                paste0("Var_",  seq_len(n_variants))))
  )
}

# ── Load data ─────────────────────────────────────────────────────────────────
if (exists("seurat_obj")) {
  af_mat <- get_layer(seurat_obj, "AF")
  dp_mat <- get_layer(seurat_obj, "DP")
  gq_mat <- get_layer(seurat_obj, "GQ")
  using_demo <- FALSE
  if (is.null(af_mat) || is.null(dp_mat) || is.null(gq_mat)) {
    message("One or more layers missing — falling back to demo data.")
    demo <- make_demo_data(); af_mat <- demo$AF; dp_mat <- demo$DP; gq_mat <- demo$GQ
    using_demo <- TRUE
  }
} else {
  message("seurat_obj not found — using demo data.")
  demo <- make_demo_data(); af_mat <- demo$AF; dp_mat <- demo$DP; gq_mat <- demo$GQ
  using_demo <- TRUE
}

layer_colors <- c(AF = "#00b4d8", DP = "#ff6b35", GQ = "#a78bfa")
n_cells    <- nrow(af_mat)
n_variants <- ncol(af_mat)
cell_ids   <- rownames(af_mat) %||% paste0("Cell_", seq_len(n_cells))
var_ids    <- colnames(af_mat) %||% paste0("Var_",  seq_len(n_variants))

# ── Pre-compute all per-cell metrics ─────────────────────────────────────────
message("Computing per-cell metrics…")
q25 <- function(x) as.numeric(quantile(x, 0.25, na.rm = TRUE))
q75 <- function(x) as.numeric(quantile(x, 0.75, na.rm = TRUE))
cv  <- function(x) { m <- mean(x, na.rm=TRUE); if (is.na(m)||m==0) NA else sd(x,na.rm=TRUE)/m }
skw <- function(x) {
  x <- x[!is.na(x)]; n <- length(x)
  if (n < 3) return(NA)
  m <- mean(x); s <- sd(x); if (s == 0) return(NA)
  sum((x-m)^3) / ((n-1)*s^3)
}

cell_metrics_raw <- data.frame(
  Cell = cell_ids,
  # AF — Central tendency
  AF_Mean        = rowMeans(af_mat, na.rm=TRUE),
  AF_Median      = apply(af_mat, 1, median, na.rm=TRUE),
  AF_Q25         = apply(af_mat, 1, q25),
  AF_Q75         = apply(af_mat, 1, q75),
  # AF — Spread / shape
  AF_SD          = apply(af_mat, 1, sd,  na.rm=TRUE),
  AF_CV          = apply(af_mat, 1, cv),
  AF_IQR         = apply(af_mat, 1, IQR, na.rm=TRUE),
  AF_Range       = apply(af_mat, 1, function(x) diff(range(x, na.rm=TRUE))),
  AF_Skew        = apply(af_mat, 1, skw),
  AF_Min         = apply(af_mat, 1, min, na.rm=TRUE),
  AF_Max         = apply(af_mat, 1, max, na.rm=TRUE),
  # AF — Genotype fractions
  AF_Hom_ref_pct = 100 * rowMeans(af_mat == 0,              na.rm=TRUE),
  AF_Het_pct     = 100 * rowMeans(af_mat > 0.2 & af_mat < 0.8, na.rm=TRUE),
  AF_Hom_alt_pct = 100 * rowMeans(af_mat >= 0.8,            na.rm=TRUE),
  AF_n_callable  = rowSums(!is.na(af_mat)),
  AF_NA_pct      = 100 * rowMeans(is.na(af_mat)),
  # DP — Coverage
  DP_Mean        = rowMeans(dp_mat, na.rm=TRUE),
  DP_Median      = apply(dp_mat, 1, median, na.rm=TRUE),
  DP_SD          = apply(dp_mat, 1, sd,  na.rm=TRUE),
  DP_CV          = apply(dp_mat, 1, cv),
  DP_IQR         = apply(dp_mat, 1, IQR, na.rm=TRUE),
  DP_Q25         = apply(dp_mat, 1, q25),
  DP_Q75         = apply(dp_mat, 1, q75),
  DP_Min         = apply(dp_mat, 1, min, na.rm=TRUE),
  DP_Max         = apply(dp_mat, 1, max, na.rm=TRUE),
  DP_Total       = rowSums(dp_mat, na.rm=TRUE),
  DP_n_low       = rowSums(dp_mat < 10,  na.rm=TRUE),
  DP_n_high      = rowSums(dp_mat > 100, na.rm=TRUE),
  DP_NA_pct      = 100 * rowMeans(is.na(dp_mat)),
  # GQ — Genotype quality
  GQ_Mean        = rowMeans(gq_mat, na.rm=TRUE),
  GQ_Median      = apply(gq_mat, 1, median, na.rm=TRUE),
  GQ_SD          = apply(gq_mat, 1, sd,  na.rm=TRUE),
  GQ_IQR         = apply(gq_mat, 1, IQR, na.rm=TRUE),
  GQ_Q25         = apply(gq_mat, 1, q25),
  GQ_Q75         = apply(gq_mat, 1, q75),
  GQ_Min         = apply(gq_mat, 1, min, na.rm=TRUE),
  GQ_n_hq        = rowSums(gq_mat >= 30, na.rm=TRUE),
  GQ_n_lq        = rowSums(gq_mat <  20, na.rm=TRUE),
  GQ_pct_hq      = 100 * rowMeans(gq_mat >= 30, na.rm=TRUE),
  GQ_NA_pct      = 100 * rowMeans(is.na(gq_mat))
)

# Cross-layer per-cell correlations
cell_metrics_raw$Cross_AF_DP_cor <- sapply(seq_len(n_cells), function(i)
  tryCatch(cor(af_mat[i,], dp_mat[i,], use="complete.obs"), error=function(e) NA))
cell_metrics_raw$Cross_AF_GQ_cor <- sapply(seq_len(n_cells), function(i)
  tryCatch(cor(af_mat[i,], gq_mat[i,], use="complete.obs"), error=function(e) NA))
cell_metrics_raw$Cross_DP_GQ_cor <- sapply(seq_len(n_cells), function(i)
  tryCatch(cor(dp_mat[i,], gq_mat[i,], use="complete.obs"), error=function(e) NA))
cell_metrics_raw$Cross_n_callable_all <- rowSums(!is.na(af_mat) & !is.na(dp_mat) & !is.na(gq_mat))

message("Done.")

# ── Metric catalogue ─────────────────────────────────────────────────────────
metric_groups <- list(
  "AF — Central Tendency" = c(
    "AF Mean"="AF_Mean", "AF Median"="AF_Median",
    "AF Q25"="AF_Q25",   "AF Q75"="AF_Q75"),
  "AF — Spread / Shape" = c(
    "AF SD"="AF_SD", "AF CV"="AF_CV", "AF IQR"="AF_IQR",
    "AF Range"="AF_Range", "AF Skewness"="AF_Skew",
    "AF Min"="AF_Min",  "AF Max"="AF_Max"),
  "AF — Genotype Fractions" = c(
    "AF Hom-Ref %"="AF_Hom_ref_pct", "AF Het %"="AF_Het_pct",
    "AF Hom-Alt %"="AF_Hom_alt_pct",
    "AF Callable Sites"="AF_n_callable", "AF Missing %"="AF_NA_pct"),
  "DP — Coverage" = c(
    "DP Mean"="DP_Mean", "DP Median"="DP_Median",
    "DP SD"="DP_SD",     "DP CV"="DP_CV",  "DP IQR"="DP_IQR",
    "DP Q25"="DP_Q25",   "DP Q75"="DP_Q75",
    "DP Total"="DP_Total","DP Min"="DP_Min","DP Max"="DP_Max",
    "DP n<10x"="DP_n_low","DP n>100x"="DP_n_high","DP Missing %"="DP_NA_pct"),
  "GQ — Genotype Quality" = c(
    "GQ Mean"="GQ_Mean", "GQ Median"="GQ_Median",
    "GQ SD"="GQ_SD",     "GQ IQR"="GQ_IQR",
    "GQ Q25"="GQ_Q25",   "GQ Q75"="GQ_Q75", "GQ Min"="GQ_Min",
    "GQ n HQ(≥30)"="GQ_n_hq", "GQ n LQ(<20)"="GQ_n_lq",
    "GQ % HQ"="GQ_pct_hq",    "GQ Missing %"="GQ_NA_pct"),
  "Cross-Layer" = c(
    "AF–DP corr"="Cross_AF_DP_cor",
    "AF–GQ corr"="Cross_AF_GQ_cor",
    "DP–GQ corr"="Cross_DP_GQ_cor",
    "Callable (all 3)"="Cross_n_callable_all")
)

all_metric_choices <- unlist(metric_groups)
inv_labels <- setNames(names(all_metric_choices), all_metric_choices)  # col→label

profile_metrics <- c("AF_Mean","AF_SD","AF_Het_pct","AF_Hom_ref_pct","AF_Hom_alt_pct",
                     "DP_Mean","DP_CV","DP_n_low","GQ_Mean","GQ_pct_hq",
                     "GQ_n_lq","Cross_n_callable_all")

# ── UI ────────────────────────────────────────────────────────────────────────
ui <- page_navbar(
  title = tags$span(tags$b("Seurat NGT"), tags$span(" Layer Explorer", style="font-weight:300")),
  theme = bs_theme(
    bg="#0d1117", fg="#e2e8f0",
    primary="#00b4d8", secondary="#a78bfa",
    base_font    = font_google("JetBrains Mono"),
    heading_font = font_google("Fraunces"),
    bootswatch   = "darkly"
  ),
  fillable = TRUE,
  
  # ════ TAB 1: Overview ════════════════════════════════════════════════════
  nav_panel("📊 Overview",
            layout_sidebar(
              sidebar = sidebar(title="Controls", width=280,
                                if (using_demo) tags$div(class="alert alert-warning", style="font-size:0.75rem;",
                                                         "⚠️ Demo data — load seurat_obj to use real data."),
                                tags$hr(),
                                tags$b("Dataset"),
                                tags$p(style="font-size:0.8rem;color:#8b949e;",
                                       paste(n_cells,"cells ×",n_variants,"variants")),
                                tags$hr(),
                                selectInput("overview_layer","Layer", choices=c("AF","DP","GQ"), selected="AF"),
                                sliderInput("hist_bins","Histogram bins",10,150,60),
                                checkboxInput("log_scale","Log Y-axis",FALSE),
                                tags$hr(),
                                sliderInput("gq_filter","Min GQ filter",0,99,0,step=1),
                                sliderInput("dp_filter","Min DP filter",0,200,0,step=1),
                                tags$hr(),
                                actionButton("reset_filters","Reset Filters",class="btn-outline-secondary btn-sm w-100")
              ),
              layout_columns(
                col_widths=c(6,6,12),
                card(card_header("Distribution"),           plotlyOutput("hist_plot",   height="280px")),
                card(card_header("Violin — per-cell mean"), plotlyOutput("violin_plot", height="280px")),
                card(card_header("Summary Statistics"),     DTOutput("summary_table"))
              )
            )
  ),
  
  # ════ TAB 2: Per-Cell Metrics ════════════════════════════════════════════
  nav_panel("🔬 Per-Cell Metrics",
            layout_sidebar(
              sidebar = sidebar(title="Controls", width=310,
                                
                                radioGroupButtons("percell_view","View",
                                                  choices=c("📈 Distributions","📊 Scatter","🎯 Cell Profile","📋 Table"),
                                                  direction="vertical", size="sm", status="outline-info"),
                                tags$hr(),
                                
                                # Distributions controls
                                conditionalPanel("input.percell_view == '📈 Distributions'",
                                                 pickerInput("dist_metrics","Metrics",
                                                             choices=metric_groups, selected=c("AF_Mean","DP_Mean","GQ_Mean"),
                                                             multiple=TRUE,
                                                             options=list(`live-search`=TRUE,`actions-box`=TRUE,
                                                                          `selected-text-format`="count > 2")),
                                                 selectInput("dist_plot_type","Plot type",
                                                             choices=c("Violin","Histogram","Box","ECDF")),
                                                 checkboxInput("dist_free_scales","Free scales per metric",TRUE)
                                ),
                                
                                # Scatter controls
                                conditionalPanel("input.percell_view == '📊 Scatter'",
                                                 selectInput("mm_x","X axis",    choices=all_metric_choices, selected="DP_Mean"),
                                                 selectInput("mm_y","Y axis",    choices=all_metric_choices, selected="AF_Mean"),
                                                 selectInput("mm_col","Color by",choices=all_metric_choices, selected="GQ_Mean"),
                                                 selectInput("mm_size","Size by",
                                                             choices=c("Fixed"="Fixed", all_metric_choices), selected="Fixed"),
                                                 sliderInput("mm_alpha","Opacity",0.1,1,0.55,step=0.05),
                                                 checkboxInput("mm_smooth","Add OLS trend line",FALSE)
                                ),
                                
                                # Cell profile controls
                                conditionalPanel("input.percell_view == '🎯 Cell Profile'",
                                                 selectizeInput("profile_cell","Select cell", choices=cell_ids,
                                                                options=list(placeholder="Type to search…")),
                                                 tags$p(style="font-size:0.75rem;color:#8b949e;",
                                                        "Radar = z-scores vs dataset. Bar = raw value vs dataset mean.")
                                ),
                                
                                # Table controls
                                conditionalPanel("input.percell_view == '📋 Table'",
                                                 pickerInput("table_cols","Columns",
                                                             choices=metric_groups,
                                                             selected=c("AF_Mean","AF_SD","AF_Het_pct","DP_Mean","DP_Total","GQ_Mean","GQ_pct_hq"),
                                                             multiple=TRUE,
                                                             options=list(`live-search`=TRUE,`actions-box`=TRUE,
                                                                          `selected-text-format`="count > 3")),
                                                 downloadButton("dl_table","⬇ Download CSV",
                                                                class="btn-sm btn-outline-secondary w-100 mt-2")
                                ),
                                
                                tags$hr(),
                                tags$p(style="font-size:0.72rem;color:#8b949e;",
                                       paste0(length(all_metric_choices)," metrics across ",n_cells," cells."))
              ),
              
              # Main panel content
              conditionalPanel("input.percell_view == '📈 Distributions'",
                               card(card_header("Per-Cell Metric Distributions"),
                                    plotlyOutput("dist_plot", height="540px"))
              ),
              conditionalPanel("input.percell_view == '📊 Scatter'",
                               layout_columns(col_widths=c(8,4),
                                              card(card_header("Scatter"), plotlyOutput("mm_scatter",   height="480px")),
                                              card(card_header("Marginals"),plotlyOutput("mm_marginals", height="480px"))
                               )
              ),
              conditionalPanel("input.percell_view == '🎯 Cell Profile'",
                               layout_columns(col_widths=c(6,6),
                                              card(card_header("Radar — z-score vs dataset"),    plotlyOutput("radar_plot",   height="460px")),
                                              card(card_header("All Metrics — Cell vs Mean"),    plotlyOutput("profile_bar",  height="460px"))
                               )
              ),
              conditionalPanel("input.percell_view == '📋 Table'",
                               card(card_header("Per-Cell Metrics Table"),
                                    DTOutput("cell_table", height="540px"))
              )
            )
  ),
  
  # ════ TAB 3: Per-Variant ═════════════════════════════════════════════════
  nav_panel("🧬 Per-Variant",
            layout_sidebar(
              sidebar = sidebar(title="Controls", width=280,
                                sliderInput("top_n_vars","Top N variants (by AF variance)",
                                            5, min(100,n_variants), min(30,n_variants)),
                                selectInput("sort_by","Sort by", choices=c("Variance","Mean","Missing %")),
                                checkboxInput("cluster_vars","Cluster variants",TRUE),
                                selectInput("heatmap_layer","Heatmap layer", choices=c("AF","DP","GQ"))
              ),
              layout_columns(col_widths=c(6,6),
                             card(card_header("Variant Ranking (AF mean ± SD)"), plotlyOutput("variant_bar",  height="380px")),
                             card(card_header("Heatmap (cells × variants)"),     plotlyOutput("heatmap_plot", height="380px"))
              )
            )
  ),
  
  # ════ TAB 4: Correlation ═════════════════════════════════════════════════
  nav_panel("📈 Correlation",
            layout_columns(col_widths=c(5,7),
                           card(card_header("Layer Correlation Matrix"),   plotlyOutput("corr_matrix", height="400px")),
                           card(card_header("3-D Scatter (AF×DP×GQ)"),    plotlyOutput("scatter_3d",  height="400px"))
            )
  ),
  
  # ════ TAB 5: QC Flags ════════════════════════════════════════════════════
  nav_panel("🚦 QC Flags",
            layout_sidebar(
              sidebar = sidebar(title="Thresholds", width=280,
                                sliderInput("qc_min_gq","Min GQ",  0,99, 20,step=1),
                                sliderInput("qc_min_dp","Min DP",  0,200,10,step=1),
                                sliderInput("qc_max_af","Max AF",  0.05,1.0,0.95,step=0.05),
                                tags$hr(),
                                uiOutput("qc_summary_ui")
              ),
              layout_columns(col_widths=c(6,6),
                             card(card_header("Cells passing QC per variant"), plotlyOutput("qc_cells_plot", height="350px")),
                             card(card_header("Variants passing QC per cell"), plotlyOutput("qc_vars_plot",  height="350px"))
              )
            )
  )
)

# ── SERVER ────────────────────────────────────────────────────────────────────
server <- function(input, output, session) {
  
  filtered_mats <- reactive({
    af <- af_mat; dp <- dp_mat; gq <- gq_mat
    mask <- gq < input$gq_filter | dp < input$dp_filter
    af[mask] <- NA; dp[mask] <- NA; gq[mask] <- NA
    list(AF=af, DP=dp, GQ=gq)
  })
  
  observeEvent(input$reset_filters, {
    updateSliderInput(session,"gq_filter",value=0)
    updateSliderInput(session,"dp_filter",value=0)
  })
  
  lc <- function(l) layer_colors[[l]]
  
  dark_layout <- function(p, xtitle="", ytitle="") {
    p %>% layout(
      paper_bgcolor="rgba(0,0,0,0)", plot_bgcolor="rgba(0,0,0,0)",
      font=list(color="#e2e8f0", family="JetBrains Mono"),
      xaxis=list(title=xtitle, gridcolor="#1e2d45"),
      yaxis=list(title=ytitle, gridcolor="#1e2d45"))
  }
  
  # ════ Tab 1 ══════════════════════════════════════════════════════════════
  output$hist_plot <- renderPlotly({
    vals <- as.vector(filtered_mats()[[input$overview_layer]]); vals <- vals[!is.na(vals)]
    plot_ly(x=vals, type="histogram", nbinsx=input$hist_bins,
            marker=list(color=lc(input$overview_layer),
                        line=list(color="#0d1117",width=0.3))) %>%
      dark_layout(xtitle=input$overview_layer, ytitle="Count") %>%
      layout(yaxis=list(type=if(input$log_scale)"log" else "linear",
                        gridcolor="#1e2d45"), bargap=0.05)
  })
  
  output$violin_plot <- renderPlotly({
    mats <- filtered_mats()
    df <- data.frame(AF=rowMeans(mats$AF,na.rm=TRUE),
                     DP=rowMeans(mats$DP,na.rm=TRUE),
                     GQ=rowMeans(mats$GQ,na.rm=TRUE)) %>%
      pivot_longer(everything(), names_to="Layer", values_to="Value")
    plot_ly(df, x=~Layer, y=~Value, color=~Layer, type="violin",
            colors=unname(layer_colors),
            box=list(visible=TRUE), meanline=list(visible=TRUE)) %>%
      dark_layout(ytitle="Per-cell mean") %>%
      layout(showlegend=FALSE)
  })
  
  output$summary_table <- renderDT({
    mats <- filtered_mats()
    rows <- lapply(c("AF","DP","GQ"), function(l) {
      v <- as.vector(mats[[l]]); v <- v[!is.na(v)]
      data.frame(Layer=l, N=length(v), Mean=round(mean(v),3), Median=round(median(v),3),
                 SD=round(sd(v),3), Min=round(min(v),3), Max=round(max(v),3),
                 `NA%`=round(100*mean(is.na(as.vector(mats[[l]]))),1), check.names=FALSE)
    })
    do.call(rbind,rows)
  }, options=list(dom="t",pageLength=3), rownames=FALSE, style="bootstrap4")
  
  # ════ Tab 2 — Per-Cell ═══════════════════════════════════════════════════
  
  # Palette helper
  metric_pal <- function(n) colorRampPalette(
    c("#00b4d8","#a78bfa","#ff6b35","#34d399","#fbbf24","#f472b6"))(n)
  
  # ── Distributions ────────────────────────────────────────────────────────
  output$dist_plot <- renderPlotly({
    cols <- input$dist_metrics
    req(length(cols) > 0)
    pal  <- metric_pal(length(cols))
    lbls <- inv_labels[cols]
    
    if (input$dist_plot_type == "ECDF") {
      p <- plot_ly()
      for (i in seq_along(cols)) {
        v   <- sort(cell_metrics_raw[[cols[i]]]); v <- v[!is.na(v)]
        cdf <- seq_along(v)/length(v)
        p   <- add_trace(p, x=v, y=cdf, type="scatter", mode="lines",
                         name=lbls[i], line=list(color=pal[i],width=2))
      }
      return(p %>% dark_layout("Value","Cumulative proportion"))
    }
    
    if (length(cols) == 1 || !input$dist_free_scales) {
      # Overlay on same axes
      df_long <- pivot_longer(
        cell_metrics_raw[, c("Cell", cols), drop=FALSE],
        -Cell, names_to="Metric", values_to="Value")
      df_long$Label <- inv_labels[df_long$Metric]
      
      if (input$dist_plot_type == "Violin") {
        p <- plot_ly(df_long, x=~Label, y=~Value, color=~Label, type="violin",
                     colors=pal, box=list(visible=TRUE), meanline=list(visible=TRUE))
      } else if (input$dist_plot_type == "Box") {
        p <- plot_ly(df_long, x=~Label, y=~Value, color=~Label, type="box",
                     colors=pal)
      } else {
        p <- plot_ly(df_long, x=~Value, color=~Label, type="histogram",
                     colors=pal, opacity=0.7, nbinsx=40) %>%
          layout(barmode="overlay")
      }
      return(p %>% dark_layout("","Value") %>% layout(showlegend=TRUE))
    }
    
    # Faceted subplots (free scales)
    n  <- length(cols)
    nc <- min(3, n)
    nr <- ceiling(n/nc)
    figs <- lapply(seq_along(cols), function(i) {
      v   <- cell_metrics_raw[[cols[i]]]
      lbl <- lbls[i]
      if (input$dist_plot_type == "Histogram") {
        plot_ly(x=v, type="histogram", name=lbl, nbinsx=40,
                marker=list(color=pal[i]), showlegend=FALSE) %>%
          layout(xaxis=list(title=lbl,gridcolor="#1e2d45"),
                 yaxis=list(gridcolor="#1e2d45"))
      } else if (input$dist_plot_type == "Box") {
        plot_ly(y=v, type="box", name=lbl, fillcolor=pal[i],
                line=list(color=pal[i]), showlegend=FALSE) %>%
          layout(xaxis=list(title=lbl,gridcolor="#1e2d45"),
                 yaxis=list(gridcolor="#1e2d45"))
      } else {  # Violin
        plot_ly(y=v, type="violin", name=lbl, fillcolor=pal[i],
                line=list(color=pal[i]), box=list(visible=TRUE),
                meanline=list(visible=TRUE), showlegend=FALSE) %>%
          layout(xaxis=list(title=lbl,gridcolor="#1e2d45"),
                 yaxis=list(gridcolor="#1e2d45"))
      }
    })
    subplot(figs, nrows=nr, shareX=FALSE, shareY=FALSE,
            titleX=TRUE, titleY=TRUE, margin=0.06) %>%
      layout(paper_bgcolor="rgba(0,0,0,0)", plot_bgcolor="rgba(0,0,0,0)",
             font=list(color="#e2e8f0",family="JetBrains Mono"), showlegend=FALSE)
  })
  
  # ── Scatter ──────────────────────────────────────────────────────────────
  output$mm_scatter <- renderPlotly({
    req(input$mm_x, input$mm_y)
    df  <- cell_metrics_raw
    xv  <- df[[input$mm_x]]; yv <- df[[input$mm_y]]; cv <- df[[input$mm_col]]
    sz  <- if (input$mm_size == "Fixed") 5 else {
      s <- df[[input$mm_size]]; rng <- range(s,na.rm=TRUE)
      if (diff(rng)==0) rep(5,length(s)) else 3 + 11*(s-rng[1])/diff(rng) }
    p <- plot_ly(x=xv, y=yv, color=cv, text=df$Cell,
                 type="scatter", mode="markers",
                 marker=list(size=sz, opacity=input$mm_alpha,
                             colorscale="Viridis",
                             colorbar=list(title=inv_labels[input$mm_col])),
                 hovertemplate=paste0("<b>%{text}</b><br>",
                                      inv_labels[input$mm_x],": %{x:.3f}<br>",
                                      inv_labels[input$mm_y],": %{y:.3f}<extra></extra>"))
    if (input$mm_smooth) {
      fit <- tryCatch(lm(yv~xv,na.action=na.omit), error=function(e) NULL)
      if (!is.null(fit)) {
        xr <- range(xv,na.rm=TRUE); xp <- seq(xr[1],xr[2],length.out=100)
        yp <- predict(fit,newdata=data.frame(xv=xp))
        p  <- add_trace(p, x=xp, y=yp, type="scatter", mode="lines",
                        line=list(color="#ff6b35",width=2,dash="dot"),
                        name="OLS fit", showlegend=TRUE, inherit=FALSE)
      }
    }
    p %>% dark_layout(inv_labels[input$mm_x], inv_labels[input$mm_y])
  })
  
  output$mm_marginals <- renderPlotly({
    req(input$mm_x, input$mm_y)
    xv <- cell_metrics_raw[[input$mm_x]]
    yv <- cell_metrics_raw[[input$mm_y]]
    px <- plot_ly(x=xv, type="histogram", nbinsx=40, name=inv_labels[input$mm_x],
                  marker=list(color="#00b4d8"), showlegend=FALSE)
    py <- plot_ly(x=yv, type="histogram", nbinsx=40, name=inv_labels[input$mm_y],
                  marker=list(color="#a78bfa"), showlegend=FALSE)
    subplot(px, py, nrows=2, shareX=FALSE, titleX=TRUE, margin=0.08) %>%
      layout(paper_bgcolor="rgba(0,0,0,0)", plot_bgcolor="rgba(0,0,0,0)",
             font=list(color="#e2e8f0",family="JetBrains Mono"),
             xaxis=list(title=inv_labels[input$mm_x],gridcolor="#1e2d45"),
             yaxis=list(title="Count",gridcolor="#1e2d45"),
             xaxis2=list(title=inv_labels[input$mm_y],gridcolor="#1e2d45"),
             yaxis2=list(title="Count",gridcolor="#1e2d45"))
  })
  
  # ── Cell Profile ──────────────────────────────────────────────────────────
  output$radar_plot <- renderPlotly({
    req(input$profile_cell)
    row <- cell_metrics_raw[cell_metrics_raw$Cell == input$profile_cell, ]
    req(nrow(row) > 0)
    zscores <- sapply(profile_metrics, function(m) {
      v <- cell_metrics_raw[[m]]; s <- sd(v,na.rm=TRUE); mn <- mean(v,na.rm=TRUE)
      if (is.na(s)||s==0) return(0)
      (as.numeric(row[[m]]) - mn) / s
    })
    lbls  <- inv_labels[profile_metrics]
    theta <- c(lbls, lbls[1]); r_val <- c(zscores, zscores[1])
    plot_ly(type="scatterpolar", mode="lines+markers",
            r=r_val, theta=theta, fill="toself",
            fillcolor="rgba(0,180,216,0.15)",
            line=list(color="#00b4d8",width=2),
            marker=list(color="#00b4d8",size=6)) %>%
      layout(polar=list(radialaxis=list(visible=TRUE,color="#8b949e",gridcolor="#1e2d45"),
                        angularaxis=list(color="#8b949e")),
             paper_bgcolor="rgba(0,0,0,0)",
             font=list(color="#e2e8f0",family="JetBrains Mono"),
             showlegend=FALSE,
             title=list(text=paste0("<b>",input$profile_cell,"</b> z-score radar"),
                        font=list(size=12)))
  })
  
  output$profile_bar <- renderPlotly({
    req(input$profile_cell)
    row       <- cell_metrics_raw[cell_metrics_raw$Cell == input$profile_cell, ]
    req(nrow(row) > 0)
    pm        <- profile_metrics
    cell_vals <- as.numeric(row[, pm])
    mean_vals <- colMeans(cell_metrics_raw[, pm], na.rm=TRUE)
    lbls      <- inv_labels[pm]
    df <- data.frame(
      Metric = rep(lbls, 2),
      Value  = c(cell_vals, mean_vals),
      Group  = rep(c("This cell","Dataset mean"), each=length(pm))
    )
    plot_ly(df, x=~Value, y=~Metric, color=~Group, type="bar", orientation="h",
            colors=c("This cell"="#00b4d8","Dataset mean"="#ff6b35")) %>%
      layout(barmode="group",
             paper_bgcolor="rgba(0,0,0,0)", plot_bgcolor="rgba(0,0,0,0)",
             font=list(color="#e2e8f0",family="JetBrains Mono"),
             xaxis=list(title="Value",gridcolor="#1e2d45"),
             yaxis=list(title="",gridcolor="#1e2d45"),
             legend=list(orientation="h",y=-0.12))
  })
  
  # ── Table ─────────────────────────────────────────────────────────────────
  cell_table_data <- reactive({
    cols <- input$table_cols; req(length(cols)>0)
    df   <- cell_metrics_raw[, c("Cell",cols), drop=FALSE]
    names(df)[-1] <- inv_labels[cols]
    df
  })
  
  output$cell_table <- renderDT({
    df <- cell_table_data()
    datatable(df, rownames=FALSE, style="bootstrap4", filter="top",
              options=list(pageLength=15, scrollX=TRUE, dom="lftip")) %>%
      formatRound(names(df)[-1], digits=3)
  })
  
  output$dl_table <- downloadHandler(
    filename = function() paste0("cell_metrics_",Sys.Date(),".csv"),
    content  = function(file) write.csv(cell_table_data(), file, row.names=FALSE)
  )
  
  # ════ Tab 3 — Per-Variant ════════════════════════════════════════════════
  output$variant_bar <- renderPlotly({
    af  <- filtered_mats()$AF; n <- input$top_n_vars
    ord <- switch(input$sort_by,
                  "Variance"  = order(apply(af,2,var, na.rm=TRUE),decreasing=TRUE)[1:n],
                  "Mean"      = order(apply(af,2,mean,na.rm=TRUE),decreasing=TRUE)[1:n],
                  "Missing %" = order(colMeans(is.na(af)),        decreasing=TRUE)[1:n])
    sub   <- af[,ord,drop=FALSE]
    means <- colMeans(sub,na.rm=TRUE); sds <- apply(sub,2,sd,na.rm=TRUE)
    df    <- data.frame(Variant=var_ids[ord],Mean=means,SD=sds)
    df    <- df[order(df$Mean,decreasing=TRUE),]
    df$Variant <- factor(df$Variant,levels=df$Variant)
    plot_ly(df, x=~Variant, y=~Mean, type="bar",
            error_y=list(array=~SD,color="#ffffff44"),
            marker=list(color=layer_colors["AF"],line=list(color="#0d1117",width=0.5))) %>%
      dark_layout("Variant","AF Mean") %>%
      layout(xaxis=list(tickangle=-45,gridcolor="#1e2d45"))
  })
  
  output$heatmap_plot <- renderPlotly({
    mat <- filtered_mats()[[input$heatmap_layer]]
    n   <- input$top_n_vars
    top <- order(apply(af_mat,2,var,na.rm=TRUE),decreasing=TRUE)[1:n]
    sub <- mat[,top,drop=FALSE]
    if (input$cluster_vars && ncol(sub)>1)
      sub <- sub[, hclust(dist(t(sub),method="euclidean"))$order, drop=FALSE]
    pal <- switch(input$heatmap_layer, AF="Blues", DP="Oranges", GQ="Purples")
    plot_ly(z=sub, x=colnames(sub), y=cell_ids, type="heatmap", colorscale=pal,
            hovertemplate="Cell:%{y}<br>Variant:%{x}<br>Value:%{z:.3f}<extra></extra>") %>%
      dark_layout("Variant","Cell") %>%
      layout(xaxis=list(tickangle=-45),
             yaxis=list(showticklabels=n_cells<=50))
  })
  
  # ════ Tab 4 — Correlation ════════════════════════════════════════════════
  output$corr_matrix <- renderPlotly({
    mats <- filtered_mats()
    df_c <- data.frame(AF=rowMeans(mats$AF,na.rm=TRUE),
                       DP=rowMeans(mats$DP,na.rm=TRUE),
                       GQ=rowMeans(mats$GQ,na.rm=TRUE))
    cm   <- round(cor(df_c,use="pairwise.complete.obs"),3)
    plot_ly(z=cm, x=colnames(cm), y=rownames(cm), type="heatmap",
            colorscale="RdBu", zmin=-1, zmax=1,
            text=cm, texttemplate="%{text}",
            hovertemplate="%{x} vs %{y}: %{z}<extra></extra>") %>%
      dark_layout()
  })
  
  output$scatter_3d <- renderPlotly({
    mats <- filtered_mats()
    df <- data.frame(AF=rowMeans(mats$AF,na.rm=TRUE),
                     DP=rowMeans(mats$DP,na.rm=TRUE),
                     GQ=rowMeans(mats$GQ,na.rm=TRUE), cell=cell_ids)
    plot_ly(df, x=~DP, y=~AF, z=~GQ, color=~GQ,
            colors=viridis::viridis(50), type="scatter3d", mode="markers",
            marker=list(size=3,opacity=0.7), text=~cell,
            hovertemplate="<b>%{text}</b><br>DP:%{x:.1f} AF:%{y:.3f} GQ:%{z:.1f}<extra></extra>") %>%
      layout(paper_bgcolor="rgba(0,0,0,0)",
             font=list(color="#e2e8f0",family="JetBrains Mono"),
             scene=list(
               xaxis=list(title="DP",gridcolor="#1e2d45",backgroundcolor="#111827"),
               yaxis=list(title="AF",gridcolor="#1e2d45",backgroundcolor="#111827"),
               zaxis=list(title="GQ",gridcolor="#1e2d45",backgroundcolor="#111827")))
  })
  
  # ════ Tab 5 — QC ═════════════════════════════════════════════════════════
  qc_mask <- reactive({
    gq_mat >= input$qc_min_gq & dp_mat >= input$qc_min_dp & af_mat <= input$qc_max_af
  })
  
  output$qc_cells_plot <- renderPlotly({
    pct <- colMeans(qc_mask(),na.rm=TRUE)*100
    df  <- data.frame(Variant=var_ids,Pct=pct)
    df  <- df[order(df$Pct),]; df$Variant <- factor(df$Variant,levels=df$Variant)
    plot_ly(df, x=~Variant, y=~Pct, type="bar",
            marker=list(color=~Pct,colorscale="Viridis",
                        line=list(color="#0d1117",width=0.3))) %>%
      dark_layout("Variant","% Cells Passing QC") %>%
      layout(yaxis=list(range=c(0,100),gridcolor="#1e2d45"),
             xaxis=list(showticklabels=n_variants<=60,tickangle=-45,gridcolor="#1e2d45"))
  })
  
  output$qc_vars_plot <- renderPlotly({
    pct <- rowMeans(qc_mask(),na.rm=TRUE)*100
    plot_ly(x=pct, type="histogram", nbinsx=40,
            marker=list(color=layer_colors["GQ"],line=list(color="#0d1117",width=0.3))) %>%
      dark_layout("% Variants Passing QC","# Cells")
  })
  
  output$qc_summary_ui <- renderUI({
    mask   <- qc_mask()
    n_pass <- sum(rowSums(mask,na.rm=TRUE)>0)
    v_pass <- sum(colSums(mask,na.rm=TRUE)>0)
    tags$div(
      tags$hr(), tags$b("QC Summary"),
      tags$p(style="font-size:0.8rem;color:#8b949e;margin-top:0.5rem;",
             paste0("Cells: ",n_pass," / ",n_cells," (",round(100*n_pass/n_cells,1),"%)")),
      tags$p(style="font-size:0.8rem;color:#8b949e;",
             paste0("Variants: ",v_pass," / ",n_variants," (",round(100*v_pass/n_variants,1),"%)"))
    )
  })
}

shinyApp(ui=ui, server=server)



