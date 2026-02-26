## ============================================================
##  Comprehensive Variant × Cell Explorer  —  R Shiny App
##  Designed for Seurat objects with NGT assay
##  (counts / AF / DP / GQ layers, 486 variants × 7412 cells)
## ============================================================

library(shiny)
library(shinydashboard)
library(shinyWidgets)
library(DT)
library(ggplot2)
library(plotly)
library(dplyr)
library(tidyr)
library(Matrix)
library(viridis)
library(pheatmap)
library(RColorBrewer)

# ── LOAD DATA ────────────────────────────────────────────────
# Replace with your actual seurat_obj path / object name
# seurat_obj <- readRDS("seurat_obj.rds")

# ── HELPERS ──────────────────────────────────────────────────
get_layer <- function(obj, layer) {
  LayerData(obj[["NGT"]], layer = layer)
}

# ── UI ───────────────────────────────────────────────────────
ui <- dashboardPage(
  skin = "black",
  
  # ---------- HEADER ----------
  dashboardHeader(
    title = tags$span(
      style = "font-family:'Courier New',monospace; font-weight:700; letter-spacing:2px; font-size:15px;",
      "⬡ VARIANT × CELL EXPLORER"
    ),
    titleWidth = 280,
    tags$li(
      class = "dropdown",
      tags$style(HTML("
        .main-header .logo { background:#0d1117 !important; border-bottom:1px solid #30363d; }
        .main-header .navbar { background:#0d1117 !important; border-bottom:1px solid #30363d; }
        .main-sidebar { background:#161b22 !important; }
        .sidebar-menu > li > a { color:#8b949e !important; font-family:'Courier New',monospace; font-size:12px; letter-spacing:1px; }
        .sidebar-menu > li.active > a, .sidebar-menu > li:hover > a { color:#58a6ff !important; background:#1f2937 !important; }
        .content-wrapper, .right-side { background:#0d1117 !important; }
        .box { background:#161b22 !important; border:1px solid #30363d !important; border-radius:6px !important; }
        .box-header { background:#161b22 !important; color:#e6edf3 !important; border-bottom:1px solid #30363d !important; font-family:'Courier New',monospace; letter-spacing:1px; }
        .box-title { color:#e6edf3 !important; }
        .nav-tabs-custom { background:#161b22 !important; border:1px solid #30363d !important; }
        .nav-tabs-custom > .nav-tabs { background:#0d1117 !important; border-bottom:1px solid #30363d !important; }
        .nav-tabs-custom > .nav-tabs > li > a { color:#8b949e !important; font-family:'Courier New',monospace; font-size:11px; }
        .nav-tabs-custom > .nav-tabs > li.active > a { color:#58a6ff !important; background:#161b22 !important; border-top:2px solid #58a6ff !important; }
        .nav-tabs-custom > .tab-content { background:#161b22 !important; }
        .small-box { border-radius:6px !important; }
        .dataTables_wrapper { color:#e6edf3 !important; }
        table.dataTable { color:#e6edf3 !important; background:#161b22 !important; border:none !important; }
        table.dataTable thead th { background:#0d1117 !important; color:#58a6ff !important; border-bottom:1px solid #30363d !important; font-family:'Courier New',monospace; font-size:11px; }
        table.dataTable tbody tr { background:#161b22 !important; }
        table.dataTable tbody tr:hover { background:#1f2937 !important; }
        .dataTables_filter input, .dataTables_length select { background:#0d1117 !important; color:#e6edf3 !important; border:1px solid #30363d !important; border-radius:4px; }
        .dataTables_info, .dataTables_paginate { color:#8b949e !important; }
        .paginate_button { color:#8b949e !important; background:#0d1117 !important; border:1px solid #30363d !important; border-radius:4px !important; margin:2px !important; }
        .paginate_button.current { color:#58a6ff !important; border-color:#58a6ff !important; }
        select, input[type='text'], .form-control { background:#0d1117 !important; color:#e6edf3 !important; border:1px solid #30363d !important; }
        label { color:#8b949e !important; font-family:'Courier New',monospace; font-size:11px; letter-spacing:1px; }
        h4, h3 { color:#e6edf3 !important; font-family:'Courier New',monospace; }
        .info-box { background:#161b22 !important; border:1px solid #30363d !important; border-radius:6px !important; }
        .info-box-icon { border-radius:6px 0 0 6px !important; }
        .info-box-content { color:#e6edf3 !important; }
        .shiny-notification { background:#1f2937; color:#e6edf3; border:1px solid #58a6ff; border-radius:6px; font-family:'Courier New',monospace; }
        ::-webkit-scrollbar { width:6px; height:6px; }
        ::-webkit-scrollbar-track { background:#0d1117; }
        ::-webkit-scrollbar-thumb { background:#30363d; border-radius:3px; }
        ::-webkit-scrollbar-thumb:hover { background:#58a6ff; }
        .selectize-input { background:#0d1117 !important; color:#e6edf3 !important; border:1px solid #30363d !important; }
        .selectize-dropdown { background:#161b22 !important; border:1px solid #30363d !important; color:#e6edf3 !important; }
        .selectize-dropdown-content .option:hover { background:#1f2937 !important; }
        .irs--shiny .irs-bar { background:#58a6ff !important; border-top:1px solid #58a6ff !important; border-bottom:1px solid #58a6ff !important; }
        .irs--shiny .irs-handle { background:#58a6ff !important; border:2px solid #58a6ff !important; }
        .irs--shiny .irs-from, .irs--shiny .irs-to, .irs--shiny .irs-single { background:#58a6ff !important; font-family:'Courier New',monospace; font-size:10px; }
        .btn-primary { background:#58a6ff !important; border-color:#58a6ff !important; font-family:'Courier New',monospace; letter-spacing:1px; }
        .btn-default { background:#1f2937 !important; border-color:#30363d !important; color:#e6edf3 !important; font-family:'Courier New',monospace; letter-spacing:1px; }
      "))
    )
  ),
  
  # ---------- SIDEBAR ----------
  dashboardSidebar(
    width = 280,
    sidebarMenu(
      id = "sidebar",
      menuItem("▸ OVERVIEW",         tabName = "overview",    icon = icon("grip-horizontal")),
      menuItem("▸ VARIANT TABLE",    tabName = "var_table",   icon = icon("table")),
      menuItem("▸ CELL TABLE",       tabName = "cell_table",  icon = icon("microscope")),
      menuItem("▸ GENOTYPE HEATMAP", tabName = "heatmap",     icon = icon("th")),
      menuItem("▸ VAF EXPLORER",     tabName = "vaf",         icon = icon("chart-bar")),
      menuItem("▸ QUALITY METRICS",  tabName = "quality",     icon = icon("heartbeat")),
      menuItem("▸ CO-MUTATION",      tabName = "comut",       icon = icon("project-diagram")),
      menuItem("▸ GENE SUMMARY",     tabName = "gene_sum",    icon = icon("dna")),
      menuItem("▸ CELL FILTER",      tabName = "cell_filter", icon = icon("filter")),
      br(),
      tags$div(
        style = "padding:10px 15px; color:#58a6ff; font-family:'Courier New',monospace; font-size:10px; border-top:1px solid #30363d;",
        "LOAD OBJECT",
        br(),
        fileInput("rds_file", NULL,
                  accept = ".rds",
                  placeholder = "seurat_obj.rds",
                  buttonLabel = "Browse"),
        tags$small(style = "color:#8b949e;", "or set seurat_obj in server.R")
      )
    )
  ),
  
  # ---------- BODY ----------
  dashboardBody(
    tabItems(
      
      # ════════════════════════════════════════
      # 1. OVERVIEW
      # ════════════════════════════════════════
      tabItem(tabName = "overview",
              fluidRow(
                infoBoxOutput("box_cells",    width = 3),
                infoBoxOutput("box_variants", width = 3),
                infoBoxOutput("box_genes",    width = 3),
                infoBoxOutput("box_samples",  width = 3)
              ),
              fluidRow(
                box(title = "CONSEQUENCE DISTRIBUTION", width = 4, height = 380,
                    plotlyOutput("plot_consequence", height = 300)),
                box(title = "VAF DISTRIBUTION (all variants)", width = 4, height = 380,
                    plotlyOutput("plot_vaf_dist", height = 300)),
                box(title = "GENOTYPING RATE DISTRIBUTION", width = 4, height = 380,
                    plotlyOutput("plot_geno_rate", height = 300))
              ),
              fluidRow(
                box(title = "TOP MUTATED GENES", width = 6, height = 360,
                    plotlyOutput("plot_top_genes", height = 280)),
                box(title = "CELLS PER SAMPLE", width = 6, height = 360,
                    plotlyOutput("plot_cells_sample", height = 280))
              )
      ),
      
      # ════════════════════════════════════════
      # 2. VARIANT TABLE
      # ════════════════════════════════════════
      tabItem(tabName = "var_table",
              fluidRow(
                box(title = "VARIANT METADATA", width = 12,
                    fluidRow(
                      column(3, selectInput("vt_consequence", "CONSEQUENCE",
                                            choices = c("All"), selected = "All", multiple = TRUE)),
                      column(3, selectInput("vt_gene_filter", "GENE",
                                            choices = c("All"), selected = "All", multiple = TRUE)),
                      column(3, sliderInput("vt_vaf", "MIN VAF (%)", 0, 100, 0, step = 1)),
                      column(3, sliderInput("vt_geno", "MIN GENOTYPING RATE (%)", 0, 100, 0, step = 1))
                    ),
                    DTOutput("variant_table")
                )
              )
      ),
      
      # ════════════════════════════════════════
      # 3. CELL TABLE
      # ════════════════════════════════════════
      tabItem(tabName = "cell_table",
              fluidRow(
                box(title = "CELL METADATA", width = 12,
                    fluidRow(
                      column(4, selectInput("ct_sample", "SAMPLE",
                                            choices = c("All"), selected = "All", multiple = TRUE)),
                      column(4, selectInput("ct_gene", "TOP MUTATED GENE",
                                            choices = c("All"), selected = "All", multiple = TRUE)),
                      column(4, sliderInput("ct_nfeat", "nFeature_NGT RANGE", 0, 500, c(0, 500), step = 1))
                    ),
                    DTOutput("cell_table")
                )
              )
      ),
      
      # ════════════════════════════════════════
      # 4. GENOTYPE HEATMAP
      # ════════════════════════════════════════
      tabItem(tabName = "heatmap",
              fluidRow(
                box(title = "CONTROLS", width = 3,
                    selectInput("hm_layer", "LAYER",
                                choices = c("counts", "AF", "DP", "GQ"), selected = "counts"),
                    selectInput("hm_gene_sel", "FILTER BY GENE (variants)",
                                choices = c("All"), selected = "All", multiple = TRUE),
                    sliderInput("hm_max_vars", "MAX VARIANTS TO SHOW", 10, 200, 50, step = 10),
                    sliderInput("hm_max_cells", "MAX CELLS TO SHOW", 100, 2000, 500, step = 100),
                    selectInput("hm_sort_cells", "SORT CELLS BY",
                                choices = c("nFeature_NGT", "sample_name", "Gene")),
                    actionButton("hm_run", "RENDER HEATMAP", class = "btn-primary",
                                 style = "width:100%; margin-top:10px;"),
                    br(), br(),
                    tags$small(style = "color:#8b949e; font-family:'Courier New',monospace;",
                               "NGT codes: 0=WT  1=HET  2=HOM  3=Missing")
                ),
                box(title = "GENOTYPE HEATMAP", width = 9, height = 720,
                    plotOutput("heatmap_plot", height = 680))
              )
      ),
      
      # ════════════════════════════════════════
      # 5. VAF EXPLORER
      # ════════════════════════════════════════
      tabItem(tabName = "vaf",
              fluidRow(
                box(title = "VAF PER VARIANT", width = 12,
                    fluidRow(
                      column(4, selectInput("vaf_gene", "GENE",
                                            choices = c("All"), selected = "All", multiple = TRUE)),
                      column(4, selectInput("vaf_consequence", "CONSEQUENCE",
                                            choices = c("All"), selected = "All", multiple = TRUE)),
                      column(4, selectInput("vaf_plot_type", "PLOT TYPE",
                                            choices = c("Box", "Violin", "Jitter", "Histogram")))
                    ),
                    plotlyOutput("vaf_plot", height = 500)
                )
              ),
              fluidRow(
                box(title = "AF LAYER — SINGLE-CELL VAF FOR SELECTED VARIANT", width = 12,
                    fluidRow(
                      column(6, selectInput("vaf_variant_sel", "SELECT VARIANT", choices = NULL)),
                      column(3, selectInput("vaf_group_by", "GROUP CELLS BY",
                                            choices = c("sample_name", "Gene"))),
                      column(3, sliderInput("vaf_dp_min", "MIN DEPTH (DP)", 0, 200, 5, step = 1))
                    ),
                    plotlyOutput("vaf_sc_plot", height = 400)
                )
              )
      ),
      
      # ════════════════════════════════════════
      # 6. QUALITY METRICS
      # ════════════════════════════════════════
      tabItem(tabName = "quality",
              fluidRow(
                box(title = "PER-CELL QUALITY", width = 6, height = 420,
                    selectInput("qc_x", "X-AXIS",
                                choices = c("nCount_NGT", "nFeature_NGT"), selected = "nCount_NGT"),
                    selectInput("qc_col", "COLOR BY",
                                choices = c("sample_name", "Gene"), selected = "sample_name"),
                    plotlyOutput("qc_cell_plot", height = 300)
                ),
                box(title = "PER-VARIANT QUALITY", width = 6, height = 420,
                    selectInput("qc_var_x", "X-AXIS",
                                choices = c("VAF", "genotyping_rate"), selected = "VAF"),
                    selectInput("qc_var_col", "COLOR BY",
                                choices = c("CONSEQUENCE", "Class"), selected = "CONSEQUENCE"),
                    plotlyOutput("qc_var_plot", height = 300)
                )
              ),
              fluidRow(
                box(title = "DEPTH (DP) ACROSS CELLS — SELECTED VARIANT", width = 6, height = 380,
                    selectInput("qc_var_dp", "SELECT VARIANT", choices = NULL),
                    plotlyOutput("qc_dp_plot", height = 280)
                ),
                box(title = "GQ DISTRIBUTION — SELECTED VARIANT", width = 6, height = 380,
                    selectInput("qc_var_gq", "SELECT VARIANT", choices = NULL),
                    plotlyOutput("qc_gq_plot", height = 280)
                )
              )
      ),
      
      # ════════════════════════════════════════
      # 7. CO-MUTATION
      # ════════════════════════════════════════
      tabItem(tabName = "comut",
              fluidRow(
                box(title = "CONTROLS", width = 3,
                    selectInput("cm_gene_sel", "GENES TO INCLUDE",
                                choices = c("All"), selected = "All", multiple = TRUE),
                    sliderInput("cm_max_vars", "MAX VARIANTS", 10, 100, 30, step = 5),
                    sliderInput("cm_min_geno", "MIN GENOTYPING RATE (%)", 0, 100, 50, step = 5),
                    selectInput("cm_metric", "CO-MUTATION METRIC",
                                choices = c("Jaccard", "Pearson", "Spearman"),
                                selected = "Jaccard"),
                    actionButton("cm_run", "COMPUTE", class = "btn-primary",
                                 style = "width:100%; margin-top:10px;")
                ),
                box(title = "CO-MUTATION MATRIX", width = 9, height = 620,
                    plotOutput("comut_plot", height = 580))
              )
      ),
      
      # ════════════════════════════════════════
      # 8. GENE SUMMARY
      # ════════════════════════════════════════
      tabItem(tabName = "gene_sum",
              fluidRow(
                box(title = "GENE-LEVEL SUMMARY TABLE", width = 7,
                    DTOutput("gene_summary_table")
                ),
                box(title = "GENE VARIANT COUNT", width = 5,
                    plotlyOutput("gene_var_count", height = 400))
              ),
              fluidRow(
                box(title = "VAF BY GENE", width = 12,
                    plotlyOutput("gene_vaf_box", height = 380))
              )
      ),
      
      # ════════════════════════════════════════
      # 9. CELL FILTER & EXPORT
      # ════════════════════════════════════════
      tabItem(tabName = "cell_filter",
              fluidRow(
                box(title = "FILTER CELLS", width = 4,
                    tags$p(style = "color:#8b949e; font-family:'Courier New',monospace; font-size:11px;",
                           "Define a cell population by mutation status"),
                    selectInput("cf_variant", "REQUIRE MUTATION IN VARIANT",
                                choices = NULL, multiple = TRUE),
                    selectInput("cf_wt_variant", "REQUIRE WT IN VARIANT",
                                choices = NULL, multiple = TRUE),
                    selectInput("cf_sample", "SAMPLE",
                                choices = c("All"), selected = "All", multiple = TRUE),
                    sliderInput("cf_min_dp", "MIN DEPTH AT ALL SELECTED VARIANTS", 0, 200, 5),
                    actionButton("cf_apply", "APPLY FILTER", class = "btn-primary",
                                 style = "width:100%; margin-top:10px;"),
                    br(), br(),
                    downloadButton("cf_download", "DOWNLOAD CELL BARCODES",
                                   style = "width:100%; background:#238636; border-color:#238636;
                                      font-family:'Courier New',monospace;")
                ),
                box(title = "FILTERED CELL RESULTS", width = 8,
                    infoBoxOutput("cf_n_cells", width = 6),
                    infoBoxOutput("cf_pct_cells", width = 6),
                    br(),
                    DTOutput("cf_cell_table"),
                    br(),
                    box(title = "GENOTYPE BREAKDOWN IN FILTERED CELLS", width = 12,
                        plotlyOutput("cf_geno_bar", height = 300))
                )
              )
      )
      
    ) # end tabItems
  ) # end dashboardBody
) # end dashboardPage


# ── SERVER ────────────────────────────────────────────────────
server <- function(input, output, session) {
  
  # ── reactive: load object ──
  seurat_r <- reactive({
    req(input$rds_file)
    withProgress(message = "Loading Seurat object...", {
      readRDS(input$rds_file$datapath)
    })
  })
  
  # Fallback: use globally loaded seurat_obj if present
  obj <- reactive({
    if (!is.null(input$rds_file)) return(seurat_r())
    if (exists("seurat_obj", envir = .GlobalEnv)) return(get("seurat_obj", envir = .GlobalEnv))
    NULL
  })
  
  # ── reactive: core data extracts ──
  var_meta <- reactive({
    req(obj())
    as.data.frame(obj()[["NGT"]]@meta.data)
  })
  
  cell_meta <- reactive({
    req(obj())
    obj()@meta.data
  })
  
  af_mat <- reactive({
    req(obj())
    as.matrix(get_layer(obj(), "AF"))
  })
  
  dp_mat <- reactive({
    req(obj())
    as.matrix(get_layer(obj(), "DP"))
  })
  
  gq_mat <- reactive({
    req(obj())
    as.matrix(get_layer(obj(), "GQ"))
  })
  
  ngt_mat <- reactive({
    req(obj())
    as.matrix(get_layer(obj(), "counts"))
  })
  
  # ── populate select inputs ──
  observe({
    req(var_meta(), cell_meta())
    vm  <- var_meta()
    cm  <- cell_meta()
    vars <- rownames(vm)
    if (is.null(vars)) vars <- as.character(seq_len(nrow(vm)))
    
    genes_v <- if ("AA_change" %in% colnames(vm)) {
      sort(unique(sub("\\..*", "", vm$AA_change)))
    } else character(0)
    
    conseqs <- if ("CONSEQUENCE" %in% colnames(vm)) sort(unique(as.character(vm$CONSEQUENCE))) else character(0)
    samples <- if ("sample_name" %in% colnames(cm)) sort(unique(cm$sample_name)) else character(0)
    genes_c <- if ("Gene" %in% colnames(cm)) sort(unique(cm$Gene)) else character(0)
    
    updateSelectInput(session, "vt_consequence",  choices = c("All", conseqs), selected = "All")
    updateSelectInput(session, "vt_gene_filter",   choices = c("All", genes_v), selected = "All")
    updateSelectInput(session, "vaf_gene",          choices = c("All", genes_v), selected = "All")
    updateSelectInput(session, "vaf_consequence",   choices = c("All", conseqs), selected = "All")
    updateSelectInput(session, "hm_gene_sel",       choices = c("All", genes_v), selected = "All")
    updateSelectInput(session, "cm_gene_sel",       choices = c("All", genes_v), selected = "All")
    updateSelectInput(session, "ct_sample",         choices = c("All", samples), selected = "All")
    updateSelectInput(session, "ct_gene",           choices = c("All", genes_c), selected = "All")
    updateSelectInput(session, "cf_sample",         choices = c("All", samples), selected = "All")
    updateSelectInput(session, "vaf_variant_sel",   choices = vars)
    updateSelectInput(session, "qc_var_dp",         choices = vars)
    updateSelectInput(session, "qc_var_gq",         choices = vars)
    updateSelectInput(session, "cf_variant",        choices = vars)
    updateSelectInput(session, "cf_wt_variant",     choices = vars)
    updateSliderInput(session, "ct_nfeat",
                      max   = max(cm$nFeature_NGT, na.rm = TRUE),
                      value = c(0, max(cm$nFeature_NGT, na.rm = TRUE)))
  })
  
  # ── OVERVIEW ──────────────────────────────────────────────
  
  output$box_cells    <- renderInfoBox({
    req(cell_meta())
    infoBox("CELLS", nrow(cell_meta()), icon = icon("circle"),
            color = "blue", fill = TRUE)
  })
  output$box_variants <- renderInfoBox({
    req(var_meta())
    infoBox("VARIANTS", nrow(var_meta()), icon = icon("dna"),
            color = "teal", fill = TRUE)
  })
  output$box_genes    <- renderInfoBox({
    req(var_meta())
    vm <- var_meta()
    n <- if ("AA_change" %in% colnames(vm)) length(unique(sub("\\..*", "", vm$AA_change))) else NA
    infoBox("GENES", n, icon = icon("code-branch"), color = "green", fill = TRUE)
  })
  output$box_samples  <- renderInfoBox({
    req(cell_meta())
    cm <- cell_meta()
    n <- if ("sample_name" %in% colnames(cm)) length(unique(cm$sample_name)) else NA
    infoBox("SAMPLES", n, icon = icon("vial"), color = "yellow", fill = TRUE)
  })
  
  output$plot_consequence <- renderPlotly({
    req(var_meta())
    vm <- var_meta()
    if (!"CONSEQUENCE" %in% colnames(vm)) return(NULL)
    df <- as.data.frame(table(vm$CONSEQUENCE)) %>% arrange(desc(Freq))
    plot_ly(df, x = ~Var1, y = ~Freq, type = "bar",
            marker = list(color = viridis(nrow(df)))) %>%
      layout(paper_bgcolor = "#161b22", plot_bgcolor = "#161b22",
             font = list(color = "#e6edf3", family = "Courier New"),
             xaxis = list(title = "", tickfont = list(size = 9), gridcolor = "#30363d"),
             yaxis = list(title = "Count", gridcolor = "#30363d"),
             margin = list(b = 60))
  })
  
  output$plot_vaf_dist <- renderPlotly({
    req(var_meta())
    vm <- var_meta()
    if (!"VAF" %in% colnames(vm)) return(NULL)
    plot_ly(x = vm$VAF, type = "histogram", nbinsx = 40,
            marker = list(color = "#58a6ff", line = list(color = "#0d1117", width = 0.5))) %>%
      layout(paper_bgcolor = "#161b22", plot_bgcolor = "#161b22",
             font = list(color = "#e6edf3", family = "Courier New"),
             xaxis = list(title = "VAF (%)", gridcolor = "#30363d"),
             yaxis = list(title = "Count", gridcolor = "#30363d"))
  })
  
  output$plot_geno_rate <- renderPlotly({
    req(var_meta())
    vm <- var_meta()
    if (!"genotyping_rate" %in% colnames(vm)) return(NULL)
    plot_ly(x = vm$genotyping_rate, type = "histogram", nbinsx = 40,
            marker = list(color = "#3fb950", line = list(color = "#0d1117", width = 0.5))) %>%
      layout(paper_bgcolor = "#161b22", plot_bgcolor = "#161b22",
             font = list(color = "#e6edf3", family = "Courier New"),
             xaxis = list(title = "Genotyping Rate (%)", gridcolor = "#30363d"),
             yaxis = list(title = "Count", gridcolor = "#30363d"))
  })
  
  output$plot_top_genes <- renderPlotly({
    req(var_meta())
    vm <- var_meta()
    if (!"AA_change" %in% colnames(vm)) return(NULL)
    vm$gene <- sub("\\..*", "", vm$AA_change)
    df <- sort(table(vm$gene), decreasing = TRUE)
    df <- data.frame(gene = names(df), n = as.numeric(df))[1:min(20, length(df)), ]
    plot_ly(df, x = ~reorder(gene, n), y = ~n, type = "bar",
            marker = list(color = viridis(nrow(df))),
            orientation = "v") %>%
      layout(paper_bgcolor = "#161b22", plot_bgcolor = "#161b22",
             font = list(color = "#e6edf3", family = "Courier New"),
             xaxis = list(title = "", tickangle = -45, tickfont = list(size = 9), gridcolor = "#30363d"),
             yaxis = list(title = "# Variants", gridcolor = "#30363d"))
  })
  
  output$plot_cells_sample <- renderPlotly({
    req(cell_meta())
    cm <- cell_meta()
    if (!"sample_name" %in% colnames(cm)) return(NULL)
    df <- as.data.frame(table(cm$sample_name)) %>% arrange(desc(Freq))
    plot_ly(df, x = ~Var1, y = ~Freq, type = "bar",
            marker = list(color = "#f78166")) %>%
      layout(paper_bgcolor = "#161b22", plot_bgcolor = "#161b22",
             font = list(color = "#e6edf3", family = "Courier New"),
             xaxis = list(title = "", tickangle = -45, tickfont = list(size = 9), gridcolor = "#30363d"),
             yaxis = list(title = "# Cells", gridcolor = "#30363d"))
  })
  
  # ── VARIANT TABLE ─────────────────────────────────────────
  
  filtered_variants <- reactive({
    req(var_meta())
    vm <- var_meta()
    if (!("All" %in% input$vt_consequence) && !is.null(input$vt_consequence))
      vm <- vm[as.character(vm$CONSEQUENCE) %in% input$vt_consequence, , drop = FALSE]
    if (!("All" %in% input$vt_gene_filter) && !is.null(input$vt_gene_filter)) {
      vm$gene <- sub("\\..*", "", vm$AA_change)
      vm <- vm[vm$gene %in% input$vt_gene_filter, , drop = FALSE]
    }
    if ("VAF" %in% colnames(vm))
      vm <- vm[vm$VAF >= input$vt_vaf, , drop = FALSE]
    if ("genotyping_rate" %in% colnames(vm))
      vm <- vm[vm$genotyping_rate >= input$vt_geno, , drop = FALSE]
    vm
  })
  
  output$variant_table <- renderDT({
    req(filtered_variants())
    datatable(filtered_variants(),
              options = list(pageLength = 15, scrollX = TRUE,
                             dom = "Bfrtip",
                             columnDefs = list(list(className = "dt-center", targets = "_all"))),
              rownames = TRUE,
              class = "cell-border stripe hover") %>%
      formatStyle(columns = names(filtered_variants()),
                  color = "#e6edf3", backgroundColor = "#161b22") %>%
      formatStyle("VAF", background = styleColorBar(c(0, 100), "#58a6ff"),
                  backgroundSize = "100% 70%", backgroundRepeat = "no-repeat",
                  backgroundPosition = "center")
  })
  
  # ── CELL TABLE ────────────────────────────────────────────
  
  filtered_cells <- reactive({
    req(cell_meta())
    cm <- cell_meta()
    if (!("All" %in% input$ct_sample) && !is.null(input$ct_sample))
      cm <- cm[cm$sample_name %in% input$ct_sample, , drop = FALSE]
    if (!("All" %in% input$ct_gene) && !is.null(input$ct_gene))
      cm <- cm[cm$Gene %in% input$ct_gene, , drop = FALSE]
    cm <- cm[cm$nFeature_NGT >= input$ct_nfeat[1] & cm$nFeature_NGT <= input$ct_nfeat[2], , drop = FALSE]
    cm
  })
  
  output$cell_table <- renderDT({
    req(filtered_cells())
    datatable(filtered_cells(),
              options = list(pageLength = 15, scrollX = TRUE,
                             columnDefs = list(list(className = "dt-center", targets = "_all"))),
              rownames = TRUE,
              class = "cell-border stripe hover") %>%
      formatStyle(columns = names(filtered_cells()),
                  color = "#e6edf3", backgroundColor = "#161b22")
  })
  
  # ── HEATMAP ───────────────────────────────────────────────
  
  heatmap_data <- eventReactive(input$hm_run, {
    req(obj(), var_meta(), cell_meta())
    mat <- get_layer(obj(), input$hm_layer)
    vm  <- var_meta()
    cm  <- cell_meta()
    
    # variant filter
    if (!("All" %in% input$hm_gene_sel) && !is.null(input$hm_gene_sel)) {
      vm$gene <- sub("\\..*", "", vm$AA_change)
      keep_vars <- which(vm$gene %in% input$hm_gene_sel)
      mat <- mat[keep_vars, , drop = FALSE]
      vm  <- vm[keep_vars, , drop = FALSE]
    }
    # subsample
    n_vars  <- min(input$hm_max_vars, nrow(mat))
    n_cells <- min(input$hm_max_cells, ncol(mat))
    
    # sort variants by VAF
    if ("VAF" %in% colnames(vm)) {
      vi <- order(vm$VAF, decreasing = TRUE)[1:n_vars]
    } else {
      vi <- seq_len(n_vars)
    }
    mat <- mat[vi, , drop = FALSE]
    vm  <- vm[vi, , drop = FALSE]
    var_names <- if ("AA_change" %in% colnames(vm)) vm$AA_change else rownames(vm)
    
    # sort cells
    sort_col <- input$hm_sort_cells
    if (sort_col %in% colnames(cm)) {
      ci <- order(cm[, sort_col])[1:n_cells]
    } else {
      ci <- seq_len(n_cells)
    }
    mat   <- mat[, ci, drop = FALSE]
    cm_s  <- cm[ci, , drop = FALSE]
    cell_names <- rownames(cm_s)
    
    rownames(mat) <- var_names
    colnames(mat) <- cell_names
    list(mat = mat, vm = vm, cm = cm_s)
  })
  
  output$heatmap_plot <- renderPlot({
    req(heatmap_data())
    d   <- heatmap_data()
    mat <- d$mat
    cm  <- d$cm
    
    # annotation
    ann_col <- data.frame(row.names = colnames(mat))
    if ("sample_name" %in% colnames(cm)) ann_col$Sample <- cm$sample_name
    if (nrow(ann_col) == 0) ann_col <- NA
    
    col_pal <- if (input$hm_layer == "counts") {
      c("0" = "#1f2937", "1" = "#3fb950", "2" = "#f85149", "3" = "#8b949e")
    } else {
      colorRampPalette(c("#0d1117", "#58a6ff", "#f85149"))(100)
    }
    
    breaks <- if (input$hm_layer == "counts") {
      c(0, 1, 2, 3)
    } else {
      seq(min(mat, na.rm = TRUE), max(mat, na.rm = TRUE), length.out = 100)
    }
    
    pheatmap(mat,
             color         = if (input$hm_layer == "counts") col_pal else col_pal,
             breaks        = if (input$hm_layer == "counts") NULL else breaks,
             cluster_rows  = TRUE,
             cluster_cols  = FALSE,
             show_colnames = FALSE,
             fontsize_row  = max(5, 9 - floor(nrow(mat) / 20)),
             annotation_col= if (is.data.frame(ann_col)) ann_col else NULL,
             border_color  = NA,
             main          = paste("Layer:", input$hm_layer),
             na_col        = "#30363d",
             annotation_colors = NULL)
  }, bg = "#161b22")
  
  # ── VAF EXPLORER ─────────────────────────────────────────
  
  vaf_filtered <- reactive({
    req(var_meta())
    vm <- var_meta()
    if (!"VAF" %in% colnames(vm)) return(vm)
    vm$gene <- sub("\\..*", "", vm$AA_change)
    if (!("All" %in% input$vaf_gene) && !is.null(input$vaf_gene))
      vm <- vm[vm$gene %in% input$vaf_gene, ]
    if (!("All" %in% input$vaf_consequence) && !is.null(input$vaf_consequence))
      vm <- vm[as.character(vm$CONSEQUENCE) %in% input$vaf_consequence, ]
    vm
  })
  
  output$vaf_plot <- renderPlotly({
    req(vaf_filtered())
    vm <- vaf_filtered()
    if (!"VAF" %in% colnames(vm)) return(NULL)
    if ("AA_change" %in% colnames(vm)) vm$label <- vm$AA_change
    pal <- viridis(max(length(unique(vm$gene)), 1))
    pt <- input$vaf_plot_type
    
    if (pt == "Histogram") {
      p <- plot_ly(vm, x = ~VAF, type = "histogram", nbinsx = 40,
                   marker = list(color = "#58a6ff")) %>%
        layout(xaxis = list(title = "VAF (%)"))
    } else if (pt == "Box") {
      p <- plot_ly(vm, x = ~gene, y = ~VAF, type = "box",
                   color = ~gene, colors = pal) %>%
        layout(xaxis = list(tickangle = -45))
    } else if (pt == "Violin") {
      p <- plot_ly(vm, x = ~gene, y = ~VAF, type = "violin",
                   color = ~gene, colors = pal) %>%
        layout(xaxis = list(tickangle = -45))
    } else {
      p <- plot_ly(vm, x = ~gene, y = ~VAF, type = "scatter", mode = "markers",
                   color = ~gene, colors = pal,
                   marker = list(size = 6, opacity = 0.7),
                   text = ~label) %>%
        layout(xaxis = list(tickangle = -45))
    }
    p %>% layout(paper_bgcolor = "#161b22", plot_bgcolor = "#161b22",
                 font = list(color = "#e6edf3", family = "Courier New"),
                 xaxis = list(gridcolor = "#30363d"),
                 yaxis = list(title = "VAF (%)", gridcolor = "#30363d"),
                 showlegend = FALSE)
  })
  
  output$vaf_sc_plot <- renderPlotly({
    req(af_mat(), dp_mat(), cell_meta(), input$vaf_variant_sel)
    af  <- af_mat()
    dp  <- dp_mat()
    cm  <- cell_meta()
    var <- input$vaf_variant_sel
    if (!var %in% rownames(af)) return(NULL)
    af_v <- af[var, ]
    dp_v <- dp[var, ]
    keep <- dp_v >= input$vaf_dp_min
    df <- data.frame(AF = af_v[keep],
                     DP = dp_v[keep],
                     group = cm[keep, input$vaf_group_by, drop = TRUE],
                     cell = names(af_v)[keep])
    pal <- viridis(length(unique(df$group)))
    plot_ly(df, x = ~AF, color = ~group, colors = pal,
            type = "histogram", nbinsx = 50,
            text = ~cell) %>%
      layout(barmode = "overlay",
             paper_bgcolor = "#161b22", plot_bgcolor = "#161b22",
             font = list(color = "#e6edf3", family = "Courier New"),
             xaxis = list(title = paste("AF (%) —", var), gridcolor = "#30363d"),
             yaxis = list(title = "Cells", gridcolor = "#30363d"),
             legend = list(font = list(size = 10)))
  })
  
  # ── QUALITY METRICS ───────────────────────────────────────
  
  output$qc_cell_plot <- renderPlotly({
    req(cell_meta())
    cm <- cell_meta()
    x_col <- input$qc_x
    c_col <- input$qc_col
    if (!x_col %in% colnames(cm)) return(NULL)
    pal <- viridis(length(unique(cm[[c_col]])))
    plot_ly(cm, x = ~get(x_col), color = ~get(c_col), colors = pal,
            type = "histogram", nbinsx = 50) %>%
      layout(barmode = "overlay",
             paper_bgcolor = "#161b22", plot_bgcolor = "#161b22",
             font = list(color = "#e6edf3", family = "Courier New"),
             xaxis = list(title = x_col, gridcolor = "#30363d"),
             yaxis = list(title = "Cells", gridcolor = "#30363d"))
  })
  
  output$qc_var_plot <- renderPlotly({
    req(var_meta())
    vm <- var_meta()
    x_col <- input$qc_var_x
    c_col <- input$qc_var_col
    if (!x_col %in% colnames(vm) || !c_col %in% colnames(vm)) return(NULL)
    pal <- viridis(length(unique(vm[[c_col]])))
    plot_ly(vm, x = ~get(x_col), color = ~as.character(get(c_col)), colors = pal,
            type = "histogram", nbinsx = 40) %>%
      layout(barmode = "stack",
             paper_bgcolor = "#161b22", plot_bgcolor = "#161b22",
             font = list(color = "#e6edf3", family = "Courier New"),
             xaxis = list(title = x_col, gridcolor = "#30363d"),
             yaxis = list(title = "Variants", gridcolor = "#30363d"))
  })
  
  output$qc_dp_plot <- renderPlotly({
    req(dp_mat(), input$qc_var_dp)
    dp  <- dp_mat()
    var <- input$qc_var_dp
    if (!var %in% rownames(dp)) return(NULL)
    dp_v <- dp[var, ]
    plot_ly(x = dp_v, type = "histogram", nbinsx = 50,
            marker = list(color = "#58a6ff")) %>%
      layout(paper_bgcolor = "#161b22", plot_bgcolor = "#161b22",
             font = list(color = "#e6edf3", family = "Courier New"),
             xaxis = list(title = paste("Depth —", var), gridcolor = "#30363d"),
             yaxis = list(title = "Cells", gridcolor = "#30363d"))
  })
  
  output$qc_gq_plot <- renderPlotly({
    req(gq_mat(), input$qc_var_gq)
    gq  <- gq_mat()
    var <- input$qc_var_gq
    if (!var %in% rownames(gq)) return(NULL)
    gq_v <- gq[var, ]
    plot_ly(x = gq_v, type = "histogram", nbinsx = 50,
            marker = list(color = "#d2a8ff")) %>%
      layout(paper_bgcolor = "#161b22", plot_bgcolor = "#161b22",
             font = list(color = "#e6edf3", family = "Courier New"),
             xaxis = list(title = paste("GQ —", var), gridcolor = "#30363d"),
             yaxis = list(title = "Cells", gridcolor = "#30363d"))
  })
  
  # ── CO-MUTATION ───────────────────────────────────────────
  
  comut_data <- eventReactive(input$cm_run, {
    req(ngt_mat(), var_meta())
    mat <- ngt_mat()  # 0=WT 1=HET 2=HOM 3=Missing
    vm  <- var_meta()
    
    # filter genotyping rate
    if ("genotyping_rate" %in% colnames(vm))
      keep <- vm$genotyping_rate >= input$cm_min_geno
    else
      keep <- rep(TRUE, nrow(vm))
    
    # gene filter
    if (!("All" %in% input$cm_gene_sel) && !is.null(input$cm_gene_sel)) {
      vm$gene <- sub("\\..*", "", vm$AA_change)
      keep <- keep & (vm$gene %in% input$cm_gene_sel)
    }
    mat <- mat[keep, , drop = FALSE]
    vm  <- vm[keep, , drop = FALSE]
    
    n <- min(input$cm_max_vars, nrow(mat))
    # rank by VAF
    if ("VAF" %in% colnames(vm)) idx <- order(vm$VAF, decreasing = TRUE)[1:n] else idx <- 1:n
    mat <- mat[idx, , drop = FALSE]
    vm  <- vm[idx, , drop = FALSE]
    
    # binary: mutant = 1/2, WT = 0, NA for 3
    bin <- mat
    bin[bin == 3] <- NA
    bin[bin > 0]  <- 1
    
    var_labels <- if ("AA_change" %in% colnames(vm)) vm$AA_change else paste0("V", seq_len(nrow(vm)))
    rownames(bin) <- var_labels
    
    metric <- input$cm_metric
    if (metric == "Pearson") {
      co <- cor(t(bin), use = "pairwise.complete.obs", method = "pearson")
    } else if (metric == "Spearman") {
      co <- cor(t(bin), use = "pairwise.complete.obs", method = "spearman")
    } else {
      # Jaccard
      n_vars2 <- nrow(bin)
      co <- matrix(NA, n_vars2, n_vars2, dimnames = list(var_labels, var_labels))
      for (i in seq_len(n_vars2)) {
        for (j in seq_len(n_vars2)) {
          a <- bin[i, ]; b <- bin[j, ]
          ok <- !is.na(a) & !is.na(b)
          co[i, j] <- sum(a[ok] == 1 & b[ok] == 1, na.rm = TRUE) /
            max(sum(a[ok] == 1 | b[ok] == 1, na.rm = TRUE), 1)
        }
      }
    }
    co
  })
  
  output$comut_plot <- renderPlot({
    req(comut_data())
    co <- comut_data()
    pal <- colorRampPalette(c("#0d1117", "#58a6ff", "#f0e68c", "#f85149"))(100)
    pheatmap(co,
             color        = pal,
             cluster_rows = TRUE, cluster_cols = TRUE,
             fontsize_row = max(5, 9 - floor(nrow(co) / 10)),
             fontsize_col = max(5, 9 - floor(ncol(co) / 10)),
             border_color = NA,
             main         = paste(input$cm_metric, "co-mutation matrix"),
             na_col       = "#30363d")
  }, bg = "#161b22")
  
  # ── GENE SUMMARY ──────────────────────────────────────────
  
  gene_summary <- reactive({
    req(var_meta())
    vm <- var_meta()
    if (!"AA_change" %in% colnames(vm)) return(NULL)
    vm$gene <- sub("\\..*", "", vm$AA_change)
    vm %>%
      group_by(gene) %>%
      summarise(
        n_variants    = n(),
        mean_VAF      = round(mean(VAF, na.rm = TRUE), 2),
        median_VAF    = round(median(VAF, na.rm = TRUE), 2),
        max_VAF       = round(max(VAF, na.rm = TRUE), 2),
        mean_geno_rate= round(mean(genotyping_rate, na.rm = TRUE), 2),
        consequences  = paste(sort(unique(as.character(CONSEQUENCE))), collapse = ", ")
      ) %>%
      arrange(desc(n_variants))
  })
  
  output$gene_summary_table <- renderDT({
    req(gene_summary())
    datatable(gene_summary(),
              options = list(pageLength = 15, scrollX = TRUE),
              rownames = FALSE,
              class = "cell-border stripe hover") %>%
      formatStyle(columns = names(gene_summary()),
                  color = "#e6edf3", backgroundColor = "#161b22") %>%
      formatStyle("mean_VAF", background = styleColorBar(c(0, 100), "#58a6ff"),
                  backgroundSize = "100% 70%", backgroundRepeat = "no-repeat",
                  backgroundPosition = "center")
  })
  
  output$gene_var_count <- renderPlotly({
    req(gene_summary())
    gs <- gene_summary()
    gs <- gs[order(gs$n_variants, decreasing = TRUE), ][1:min(20, nrow(gs)), ]
    plot_ly(gs, x = ~reorder(gene, n_variants), y = ~n_variants, type = "bar",
            marker = list(color = viridis(nrow(gs)))) %>%
      layout(paper_bgcolor = "#161b22", plot_bgcolor = "#161b22",
             font = list(color = "#e6edf3", family = "Courier New"),
             xaxis = list(title = "", tickangle = -45, tickfont = list(size = 9), gridcolor = "#30363d"),
             yaxis = list(title = "# Variants", gridcolor = "#30363d"))
  })
  
  output$gene_vaf_box <- renderPlotly({
    req(var_meta())
    vm <- var_meta()
    if (!"VAF" %in% colnames(vm)) return(NULL)
    vm$gene <- sub("\\..*", "", vm$AA_change)
    top_genes <- names(sort(table(vm$gene), decreasing = TRUE))[1:min(20, length(unique(vm$gene)))]
    vm <- vm[vm$gene %in% top_genes, ]
    pal <- viridis(length(top_genes))
    plot_ly(vm, x = ~gene, y = ~VAF, type = "box",
            color = ~gene, colors = pal) %>%
      layout(paper_bgcolor = "#161b22", plot_bgcolor = "#161b22",
             font = list(color = "#e6edf3", family = "Courier New"),
             xaxis = list(title = "", tickangle = -45, tickfont = list(size = 9), gridcolor = "#30363d"),
             yaxis = list(title = "VAF (%)", gridcolor = "#30363d"),
             showlegend = FALSE)
  })
  
  # ── CELL FILTER ───────────────────────────────────────────
  
  cf_result <- eventReactive(input$cf_apply, {
    req(ngt_mat(), dp_mat(), cell_meta())
    ngt <- ngt_mat()
    dp  <- dp_mat()
    cm  <- cell_meta()
    cells <- colnames(ngt)
    
    keep <- rep(TRUE, length(cells))
    
    # sample filter
    if (!("All" %in% input$cf_sample) && !is.null(input$cf_sample))
      keep <- keep & (cm$sample_name %in% input$cf_sample)
    
    # mutant variants
    if (!is.null(input$cf_variant) && length(input$cf_variant) > 0) {
      for (v in input$cf_variant) {
        if (!v %in% rownames(ngt)) next
        ngt_v <- ngt[v, ]
        dp_v  <- dp[v, ]
        keep  <- keep & (ngt_v %in% c(1, 2)) & (dp_v >= input$cf_min_dp)
      }
    }
    
    # WT variants
    if (!is.null(input$cf_wt_variant) && length(input$cf_wt_variant) > 0) {
      for (v in input$cf_wt_variant) {
        if (!v %in% rownames(ngt)) next
        ngt_v <- ngt[v, ]
        dp_v  <- dp[v, ]
        keep  <- keep & (ngt_v == 0) & (dp_v >= input$cf_min_dp)
      }
    }
    
    list(cells = cells[keep], cm = cm[keep, , drop = FALSE],
         ngt = ngt[, keep, drop = FALSE])
  })
  
  output$cf_n_cells <- renderInfoBox({
    infoBox("FILTERED CELLS", if (!is.null(cf_result())) length(cf_result()$cells) else 0,
            icon = icon("filter"), color = "blue", fill = TRUE)
  })
  
  output$cf_pct_cells <- renderInfoBox({
    req(cell_meta())
    total <- nrow(cell_meta())
    n <- if (!is.null(cf_result())) length(cf_result()$cells) else 0
    pct <- round(n / max(total, 1) * 100, 1)
    infoBox("% OF ALL CELLS", paste0(pct, "%"),
            icon = icon("percent"), color = "green", fill = TRUE)
  })
  
  output$cf_cell_table <- renderDT({
    req(cf_result())
    datatable(cf_result()$cm,
              options = list(pageLength = 10, scrollX = TRUE),
              rownames = TRUE,
              class = "cell-border stripe hover") %>%
      formatStyle(columns = names(cf_result()$cm),
                  color = "#e6edf3", backgroundColor = "#161b22")
  })
  
  output$cf_geno_bar <- renderPlotly({
    req(cf_result())
    ngt <- cf_result()$ngt
    df <- data.frame(
      code  = as.vector(ngt),
      variant = rep(rownames(ngt), ncol(ngt))
    )
    df$genotype <- factor(df$code,
                          levels = 0:3,
                          labels = c("WT", "HET", "HOM", "Missing"))
    df2 <- df %>% group_by(variant, genotype) %>% summarise(n = n(), .groups = "drop")
    pal <- c("WT" = "#3fb950", "HET" = "#58a6ff", "HOM" = "#f85149", "Missing" = "#8b949e")
    plot_ly(df2, x = ~variant, y = ~n, color = ~genotype, colors = pal,
            type = "bar") %>%
      layout(barmode = "stack",
             paper_bgcolor = "#161b22", plot_bgcolor = "#161b22",
             font = list(color = "#e6edf3", family = "Courier New"),
             xaxis = list(title = "", tickangle = -45, tickfont = list(size = 8), gridcolor = "#30363d"),
             yaxis = list(title = "# Cells", gridcolor = "#30363d"))
  })
  
  output$cf_download <- downloadHandler(
    filename = function() paste0("filtered_cells_", Sys.Date(), ".csv"),
    content  = function(file) {
      req(cf_result())
      df <- data.frame(barcode = cf_result()$cells)
      write.csv(df, file, row.names = FALSE)
    }
  )
  
} # end server

shinyApp(ui, server)