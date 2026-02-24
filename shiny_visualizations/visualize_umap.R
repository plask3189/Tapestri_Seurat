## ─────────────────────────────────────────────────────────────
##  Multi-Modal Explorer  ·  DSB_norm × NGT
##  Tabs: 1) Side-by-side 3D UMAP  2) Co-mutation  3) Protein×NGT
##  Requires: shiny, plotly, Seurat, SeuratObject, dplyr, tidyr, viridis
## ─────────────────────────────────────────────────────────────

library(shiny)
library(plotly)
library(Seurat)
library(SeuratObject)
library(dplyr)
library(tidyr)
library(viridis)

# ═══════════════════════════════════════════════════════════════
#  HELPERS
# ═══════════════════════════════════════════════════════════════

`%||%` <- function(a, b) if (!is.null(a)) a else b

NGT_COLORS <- c(
  "WT (0)"      = "#2166AC",
  "Het (1)"     = "#F4A582",
  "Hom (2)"     = "#D6604D",
  "Missing (3)" = "#AAAAAA"
)

CLUSTER_COLORS <- c(
  "0"="#E63946","1"="#2A9D8F","2"="#E9C46A",
  "3"="#F4A261","4"="#457B9D","5"="#A8DADC","6"="#6A4C93"
)

get_obj <- function() {
  if (exists("seurat_obj", envir=.GlobalEnv))
    get("seurat_obj", envir=.GlobalEnv)
  else {
    p <- "seurat_obj.rds"
    if (file.exists(p)) readRDS(p) else stop("seurat_obj not found")
  }
}

build_base_df <- function(obj) {
  umap_2d <- Embeddings(obj, "umap")
  pca_emb <- Embeddings(obj, "pca")
  has_3d  <- "umap_3" %in% colnames(umap_2d)
  data.frame(
    cell    = rownames(umap_2d),
    umap_1  = umap_2d[,1],
    umap_2  = umap_2d[,2],
    umap_3  = if (has_3d) umap_2d[,3] else pca_emb[,"PC_3"],
    cluster = as.character(obj$seurat_clusters),
    sample  = obj$sample_name,
    stringsAsFactors=FALSE
  )
}

add_dsb_col <- function(df, mat, feat) {
  df$color_val   <- as.numeric(mat[feat, df$cell])
  df$color_type  <- "continuous"
  df$color_label <- feat
  df
}

add_ngt_col <- function(df, mat, feat) {
  vals <- as.integer(mat[feat, df$cell])
  lv   <- factor(vals, levels=0:3,
                 labels=c("WT (0)","Het (1)","Hom (2)","Missing (3)"))
  df$color_val   <- as.character(lv)
  df$color_type  <- "categorical"
  df$color_label <- feat
  df
}

make_axis_style <- function(ax_vis, gd_vis, title_str) {
  list(
    showline=ax_vis, showticklabels=ax_vis,
    title=list(text=title_str, font=list(color="#8b949e",size=11)),
    tickfont=list(color="#8b949e",size=9),
    gridcolor=if(gd_vis)"#30363d" else "transparent",
    zerolinecolor=if(gd_vis)"#30363d" else "transparent",
    backgroundcolor="#0d1117", showbackground=TRUE
  )
}

scatter3d_continuous <- function(df, sz, alpha) {
  plot_ly(
    data=df, x=~umap_1, y=~umap_2, z=~umap_3,
    type="scatter3d", mode="markers",
    marker=list(
      size=sz, opacity=alpha, color=~color_val,
      colorscale=list(c(0,"#0d1117"),c(0.15,"#1a4a7a"),c(0.4,"#2166AC"),
                      c(0.65,"#4dac26"),c(0.85,"#f4a40e"),c(1,"#d73027")),
      colorbar=list(
        title=list(text=df$color_label[1], font=list(color="#e6edf3",size=12)),
        tickfont=list(color="#8b949e",size=10),
        bgcolor="#161b22",bordercolor="#30363d",borderwidth=1,len=0.5
      ),
      showscale=TRUE
    ),
    text=~paste0("<b>Cell:</b> ",cell,"<br><b>Cluster:</b> ",cluster,
                 "<br><b>",color_label,":</b> ",round(color_val,2)),
    hoverinfo="text"
  )
}

scatter3d_categorical <- function(df, sz, alpha, pal) {
  cats <- unique(df$color_val)
  p    <- plot_ly(type="scatter3d", mode="markers")
  for (cat in cats) {
    sub <- df[df$color_val==cat,]
    p   <- add_trace(p,
                     x=sub$umap_1, y=sub$umap_2, z=sub$umap_3, name=cat,
                     marker=list(color=unname(pal[cat]%||%"#aaaaaa"), size=sz, opacity=alpha),
                     text=paste0("<b>Cell:</b> ",sub$cell,"<br><b>Cluster:</b> ",sub$cluster,
                                 "<br><b>",sub$color_label,":</b> ",cat),
                     hoverinfo="text"
    )
  }
  p
}

apply_3d_layout <- function(p, ax_vis, gd_vis, z_label="PC 3") {
  p %>%
    layout(
      paper_bgcolor="#0d1117", plot_bgcolor="#0d1117",
      font=list(color="#e6edf3",family="IBM Plex Mono"),
      legend=list(bgcolor="#161b22",bordercolor="#30363d",borderwidth=1,
                  font=list(color="#e6edf3",size=10)),
      scene=list(
        bgcolor="#0d1117",
        xaxis=make_axis_style(ax_vis,gd_vis,"UMAP 1"),
        yaxis=make_axis_style(ax_vis,gd_vis,"UMAP 2"),
        zaxis=make_axis_style(ax_vis,gd_vis,z_label),
        camera=list(eye=list(x=1.4,y=1.4,z=0.9))
      ),
      margin=list(l=0,r=0,t=30,b=0)
    ) %>%
    config(displaylogo=FALSE,
           modeBarButtons=list(list("zoom3d","pan3d","orbitRotation",
                                    "resetCameraDefault3d","toImage")))
}

# ═══════════════════════════════════════════════════════════════
#  CSS
# ═══════════════════════════════════════════════════════════════

APP_CSS <- "
:root {
  --bg:#0d1117; --surf:#161b22; --border:#30363d;
  --acc:#58a6ff; --acc2:#3fb950; --text:#e6edf3; --muted:#8b949e;
}
*{box-sizing:border-box;margin:0;padding:0;}
html,body{height:100%;background:var(--bg);color:var(--text);
  font-family:'IBM Plex Sans',sans-serif;font-size:13px;overflow:hidden;}

/* topbar */
.topbar{display:flex;align-items:center;gap:14px;padding:0 20px;
  background:var(--surf);border-bottom:1px solid var(--border);height:50px;}
.logo{font-family:'IBM Plex Mono',monospace;font-size:14px;font-weight:600;
  color:var(--acc);letter-spacing:.04em;}
.subtitle{color:var(--muted);font-size:11px;font-family:'IBM Plex Mono',monospace;}
.topbar-spacer{flex:1;}
.pill{background:var(--bg);border:1px solid var(--border);border-radius:20px;
  padding:2px 10px;font-family:'IBM Plex Mono',monospace;font-size:11px;color:var(--muted);}
.pill span{color:var(--acc);}

/* nav */
.nav-tabs{border-bottom:1px solid var(--border)!important;
  background:var(--surf);padding:0 20px;}
.nav-tabs>li>a{font-family:'IBM Plex Mono',monospace!important;font-size:11px!important;
  letter-spacing:.06em!important;text-transform:uppercase!important;
  color:var(--muted)!important;border:none!important;
  border-bottom:2px solid transparent!important;background:transparent!important;
  padding:10px 16px!important;margin:0!important;border-radius:0!important;}
.nav-tabs>li.active>a,.nav-tabs>li>a:hover{
  color:var(--acc)!important;border-bottom-color:var(--acc)!important;
  background:transparent!important;}
.tab-content{height:calc(100vh - 96px);overflow:hidden;}
.tab-pane{height:100%;}

/* dual layout */
.dual-layout{display:flex;height:100%;overflow:hidden;}
.dual-sidebar{width:240px;min-width:240px;background:var(--surf);
  border-right:1px solid var(--border);padding:14px;overflow-y:auto;
  display:flex;flex-direction:column;gap:12px;}
.dual-plots{flex:1;display:flex;overflow:hidden;}
.plot-half{flex:1;border-right:1px solid var(--border);position:relative;}
.plot-half:last-child{border-right:none;}
.plot-title{position:absolute;top:8px;left:50%;transform:translateX(-50%);
  font-family:'IBM Plex Mono',monospace;font-size:11px;letter-spacing:.06em;
  text-transform:uppercase;color:var(--muted);z-index:10;pointer-events:none;}

/* analytics layout */
.analytics-layout{display:flex;height:100%;overflow:hidden;}
.analytics-sidebar{width:240px;min-width:240px;background:var(--surf);
  border-right:1px solid var(--border);padding:14px;overflow-y:auto;
  display:flex;flex-direction:column;gap:12px;}
.analytics-main{flex:1;overflow-y:auto;padding:14px;
  display:flex;flex-direction:column;gap:14px;}
.analytics-row{display:flex;gap:14px;min-height:420px;}
.analytics-card{flex:1;background:var(--surf);border:1px solid var(--border);
  border-radius:8px;padding:12px;display:flex;flex-direction:column;min-width:0;}
.card-title{font-family:'IBM Plex Mono',monospace;font-size:10px;
  letter-spacing:.1em;text-transform:uppercase;color:var(--muted);margin-bottom:8px;}

/* controls */
.cs{border:1px solid var(--border);border-radius:6px;padding:10px;}
.cs-title{font-family:'IBM Plex Mono',monospace;font-size:10px;
  letter-spacing:.1em;text-transform:uppercase;color:var(--muted);margin-bottom:8px;}
.shiny-input-container{margin-bottom:7px!important;}
.shiny-input-container label{color:var(--muted)!important;font-size:11px!important;
  font-family:'IBM Plex Mono',monospace!important;margin-bottom:3px!important;}
.form-control,select.form-control{background:var(--bg)!important;
  border:1px solid var(--border)!important;color:var(--text)!important;
  border-radius:4px!important;font-size:12px!important;
  font-family:'IBM Plex Mono',monospace!important;
  padding:5px 8px!important;height:auto!important;}
.form-control:focus{border-color:var(--acc)!important;outline:none!important;box-shadow:none!important;}
input[type=range]{accent-color:var(--acc);width:100%;}
.radio label,.checkbox label{color:var(--text)!important;font-size:12px!important;}

/* ngt legend */
.ngt-legend{display:flex;gap:8px;flex-wrap:wrap;margin-top:6px;}
.ngt-sw{display:flex;align-items:center;gap:4px;
  font-family:'IBM Plex Mono',monospace;font-size:10px;color:var(--text);}
.sw-dot{width:9px;height:9px;border-radius:50%;flex-shrink:0;}

/* stats */
.shiny-text-output{background:var(--bg)!important;color:var(--muted)!important;
  font-family:'IBM Plex Mono',monospace!important;font-size:10px!important;
  border:1px solid var(--border)!important;border-radius:4px!important;
  padding:7px!important;white-space:pre-wrap!important;}

::-webkit-scrollbar{width:4px;}
::-webkit-scrollbar-track{background:var(--bg);}
::-webkit-scrollbar-thumb{background:var(--border);border-radius:2px;}
"

# ═══════════════════════════════════════════════════════════════
#  UI
# ═══════════════════════════════════════════════════════════════

ui <- fluidPage(
  tags$head(
    tags$title("Multi-Modal Explorer"),
    tags$link(rel="stylesheet",
              href="https://fonts.googleapis.com/css2?family=IBM+Plex+Mono:wght@400;600&family=IBM+Plex+Sans:wght@300;400;600&display=swap"),
    tags$style(HTML(APP_CSS))
  ),
  
  ## topbar
  div(class="topbar",
      div(class="logo","Multi-Modal Explorer"),
      div(class="subtitle","DSB\u2009\u00d7\u2009NGT"),
      div(class="topbar-spacer"),
      div(class="pill","cells: ",tags$span("7,412")),
      div(class="pill","proteins: ",tags$span("46")),
      div(class="pill","variants: ",tags$span("126"))
  ),
  
  tabsetPanel(id="main_tabs",
              
              # ── TAB 1: Side-by-side UMAP ──────────────────────────────
              tabPanel("Side-by-Side UMAP",
                       div(class="dual-layout",
                           div(class="dual-sidebar",
                               div(class="cs",
                                   div(class="cs-title","Left — DSB Protein"),
                                   selectInput("left_feat","Marker",choices=character(0))
                               ),
                               div(class="cs",
                                   div(class="cs-title","Right — NGT Variant"),
                                   selectInput("right_feat","Variant",choices=character(0)),
                                   div(class="ngt-legend",
                                       div(class="ngt-sw",div(class="sw-dot",style="background:#2166AC;"),"WT (0)"),
                                       div(class="ngt-sw",div(class="sw-dot",style="background:#F4A582;"),"Het (1)"),
                                       div(class="ngt-sw",div(class="sw-dot",style="background:#D6604D;"),"Hom (2)"),
                                       div(class="ngt-sw",div(class="sw-dot",style="background:#AAAAAA;"),"Miss (3)")
                                   )
                               ),
                               div(class="cs",
                                   div(class="cs-title","Display"),
                                   sliderInput("sb_sz","Point Size",min=1,max=8,value=3,step=.5),
                                   sliderInput("sb_alpha","Opacity",min=.05,max=1,value=.65,step=.05),
                                   checkboxInput("sb_axes","Show Axes",TRUE),
                                   checkboxInput("sb_grid","Show Grid",FALSE)
                               ),
                               div(class="cs",
                                   div(class="cs-title","Filter"),
                                   uiOutput("sb_sample_ui"),
                                   uiOutput("sb_cluster_ui")
                               ),
                               div(class="cs",
                                   div(class="cs-title","NGT Stats"),
                                   verbatimTextOutput("sb_stats")
                               )
                           ),
                           div(class="dual-plots",
                               div(class="plot-half",
                                   div(class="plot-title","DSB Protein"),
                                   plotlyOutput("left_plot",width="100%",height="100%")
                               ),
                               div(class="plot-half",
                                   div(class="plot-title","NGT Variant"),
                                   plotlyOutput("right_plot",width="100%",height="100%")
                               )
                           )
                       )
              ),
              
              # ── TAB 2: Co-Mutation ────────────────────────────────────
              tabPanel("Co-Mutation",
                       div(class="analytics-layout",
                           div(class="analytics-sidebar",
                               div(class="cs",
                                   div(class="cs-title","Metric"),
                                   selectInput("comut_metric","Metric",
                                               choices=c("Co-occurring cells"="co_occur",
                                                         "Jaccard similarity"="jaccard",
                                                         "Obs - Expected (mutex)"="mutex"),
                                               selected="co_occur")
                               ),
                               div(class="cs",
                                   div(class="cs-title","Variant Filter"),
                                   numericInput("comut_vaf_min","Min VAF %",value=0,min=0,max=100,step=1),
                                   numericInput("comut_n_top","Top N variants",value=30,min=5,max=126,step=5),
                                   checkboxInput("comut_excl_missing","Missing = not mutated",value=TRUE)
                               ),
                               div(class="cs",
                                   div(class="cs-title","Labels"),
                                   checkboxInput("comut_showvals","Show values in heatmap",value=FALSE),
                                   sliderInput("comut_fontsize","Label size",min=6,max=14,value=9,step=1)
                               )
                           ),
                           div(class="analytics-main",
                               div(class="analytics-row",
                                   div(class="analytics-card",style="flex:2;",
                                       div(class="card-title","Co-Mutation Heatmap"),
                                       plotlyOutput("comut_heatmap",width="100%",height="390px")
                                   ),
                                   div(class="analytics-card",style="flex:1;",
                                       div(class="card-title","Mutation Frequency"),
                                       plotlyOutput("comut_freq",width="100%",height="390px")
                                   )
                               ),
                               div(class="analytics-row",style="min-height:320px;",
                                   div(class="analytics-card",
                                       div(class="card-title","% Mutated Cells per Cluster (Top 10 Variants)"),
                                       plotlyOutput("comut_cluster_bar",width="100%",height="280px")
                                   )
                               )
                           )
                       )
              ),
              
              # ── TAB 3: Protein x NGT ─────────────────────────────────
              tabPanel("Protein \u00d7 NGT",
                       div(class="analytics-layout",
                           div(class="analytics-sidebar",
                               div(class="cs",
                                   div(class="cs-title","Primary Features"),
                                   selectInput("pxn_protein","Protein",choices=character(0)),
                                   selectInput("pxn_variant","Variant",choices=character(0))
                               ),
                               div(class="cs",
                                   div(class="cs-title","Violin Options"),
                                   checkboxInput("pxn_jitter","Show jitter points",value=FALSE),
                                   checkboxInput("pxn_excl_miss","Exclude Missing (3)",value=TRUE)
                               ),
                               div(class="cs",
                                   div(class="cs-title","Scatter Options"),
                                   selectInput("pxn_prot2","Protein X",choices=character(0)),
                                   selectInput("pxn_prot3","Protein Y",choices=character(0)),
                                   selectInput("pxn_color_by","Color by",
                                               choices=c("NGT Genotype"="ngt","Cluster"="cluster"))
                               ),
                               div(class="cs",
                                   div(class="cs-title","Stats"),
                                   verbatimTextOutput("pxn_stats")
                               )
                           ),
                           div(class="analytics-main",
                               div(class="analytics-row",
                                   div(class="analytics-card",style="flex:1;",
                                       div(class="card-title","Protein Expression by Genotype"),
                                       plotlyOutput("pxn_violin",width="100%",height="390px")
                                   ),
                                   div(class="analytics-card",style="flex:1;",
                                       div(class="card-title","Median Protein per Genotype \u00d7 Top 20 Variants"),
                                       plotlyOutput("pxn_heatmap",width="100%",height="390px")
                                   )
                               ),
                               div(class="analytics-row",style="min-height:340px;",
                                   div(class="analytics-card",style="flex:1;",
                                       div(class="card-title","Protein Scatter — Coloured by NGT / Cluster"),
                                       plotlyOutput("pxn_scatter",width="100%",height="300px")
                                   ),
                                   div(class="analytics-card",style="flex:1;",
                                       div(class="card-title","Mutation Load vs Protein Expression"),
                                       plotlyOutput("pxn_mutload",width="100%",height="300px")
                                   )
                               )
                           )
                       )
              )
  )
)

# ═══════════════════════════════════════════════════════════════
#  SERVER
# ═══════════════════════════════════════════════════════════════

server <- function(input, output, session) {
  
  # ── Core data ────────────────────────────────────────────────
  obj <- reactive({ get_obj() })
  
  base_df <- reactive({ build_base_df(obj()) })
  
  dsb_mat <- reactive({
    as.matrix(LayerData(obj(), assay="DSB_norm", layer="counts"))
  })
  
  ngt_mat <- reactive({
    as.matrix(LayerData(obj(), assay="NGT", layer="counts"))
  })
  
  # binary mutation matrix (variants × cells); TRUE = Het or Hom
  mut_bin <- reactive({
    mat <- ngt_mat()
    mat == 1 | mat == 2
  })
  
  # ── Populate selectors ───────────────────────────────────────
  observe({
    dsb <- sort(rownames(dsb_mat()))
    ngt <- sort(rownames(ngt_mat()))
    updateSelectInput(session,"left_feat",   choices=dsb, selected=dsb[1])
    updateSelectInput(session,"right_feat",  choices=ngt, selected=ngt[1])
    updateSelectInput(session,"pxn_protein", choices=dsb, selected=dsb[1])
    updateSelectInput(session,"pxn_variant", choices=ngt, selected=ngt[1])
    updateSelectInput(session,"pxn_prot2",   choices=dsb, selected=dsb[1])
    updateSelectInput(session,"pxn_prot3",   choices=dsb,
                      selected=dsb[min(2,length(dsb))])
  })
  
  # ── TAB 1 FILTERS ────────────────────────────────────────────
  output$sb_sample_ui <- renderUI({
    sams <- sort(unique(base_df()$sample))
    checkboxGroupInput("sb_samples","Sample",choices=sams,selected=sams)
  })
  output$sb_cluster_ui <- renderUI({
    cls <- sort(unique(base_df()$cluster))
    checkboxGroupInput("sb_clusters","Cluster",choices=cls,selected=cls)
  })
  
  filt_df <- reactive({
    df <- base_df()
    if (!is.null(input$sb_samples))  df <- df[df$sample  %in% input$sb_samples,]
    if (!is.null(input$sb_clusters)) df <- df[df$cluster %in% input$sb_clusters,]
    df
  })
  
  z_label <- reactive({
    if ("umap_3" %in% colnames(Embeddings(obj(),"umap"))) "UMAP 3" else "PC 3"
  })
  
  # ── LEFT PLOT ────────────────────────────────────────────────
  output$left_plot <- renderPlotly({
    df   <- filt_df()
    feat <- req(input$left_feat)
    df   <- add_dsb_col(df, dsb_mat(), feat)
    scatter3d_continuous(df, input$sb_sz%||%3, input$sb_alpha%||%.65) %>%
      apply_3d_layout(isTRUE(input$sb_axes), isTRUE(input$sb_grid), z_label())
  })
  
  # ── RIGHT PLOT ───────────────────────────────────────────────
  output$right_plot <- renderPlotly({
    df   <- filt_df()
    feat <- req(input$right_feat)
    df   <- add_ngt_col(df, ngt_mat(), feat)
    scatter3d_categorical(df, input$sb_sz%||%3, input$sb_alpha%||%.65, NGT_COLORS) %>%
      apply_3d_layout(isTRUE(input$sb_axes), isTRUE(input$sb_grid), z_label())
  })
  
  # ── SIDE-BY-SIDE STATS ───────────────────────────────────────
  output$sb_stats <- renderText({
    df   <- filt_df()
    feat <- input$right_feat
    if (is.null(feat) || nchar(feat)==0)
      return(sprintf("Cells shown: %d", nrow(df)))
    ngt <- as.integer(ngt_mat()[feat, df$cell])
    gf  <- factor(ngt, levels=0:3,
                  labels=c("WT(0)","Het(1)","Hom(2)","Miss(3)"))
    tbl <- table(gf)
    pct <- round(100*tbl/sum(tbl),1)
    ln  <- mapply(function(nm,ct,pc) sprintf("  %-8s %4d  %5.1f%%",nm,ct,pc),
                  names(tbl),as.integer(tbl),pct,SIMPLIFY=FALSE)
    paste0("Cells: ",nrow(df),"\n\nVariant: ",feat,"\n",
           paste(ln,collapse="\n"))
  })
  
  # ═══════════════════════════════════════════════════════════
  #  TAB 2: CO-MUTATION
  # ═══════════════════════════════════════════════════════════
  
  top_vars <- reactive({
    mat    <- ngt_mat()
    vmeta  <- obj()@assays$NGT@meta.data
    vaf_mn <- input$comut_vaf_min %||% 0
    n_top  <- input$comut_n_top   %||% 30
    freq   <- rowMeans(mat==1|mat==2, na.rm=TRUE)
    if ("VAF" %in% colnames(vmeta)) {
      keep <- rownames(mat) %in% rownames(vmeta)[vmeta$VAF >= vaf_mn]
      freq <- freq[keep]
    }
    names(sort(freq, decreasing=TRUE))[seq_len(min(n_top,length(freq)))]
  })
  
  comut_mat_data <- reactive({
    bin  <- mut_bin()
    vars <- top_vars()
    sub  <- bin[vars,,drop=FALSE]
    n    <- nrow(sub)
    nc   <- ncol(sub)
    met  <- input$comut_metric %||% "co_occur"
    m    <- matrix(0L, n, n, dimnames=list(vars,vars))
    for (i in seq_len(n)) for (j in seq_len(n)) {
      a <- sub[i,]; b <- sub[j,]
      m[i,j] <- switch(met,
                       co_occur = sum(a & b),
                       jaccard  = { u<-sum(a|b); if(u==0) 0 else round(sum(a&b)/u,3) },
                       mutex    = round(sum(a&b) - sum(a)*sum(b)/nc, 2)
      )
    }
    m
  })
  
  output$comut_heatmap <- renderPlotly({
    cm   <- comut_mat_data()
    vars <- rownames(cm)
    show <- isTRUE(input$comut_showvals)
    fs   <- input$comut_fontsize%||%9
    met  <- input$comut_metric%||%"co_occur"
    lbl  <- c(co_occur="Co-occurring cells",
              jaccard="Jaccard similarity",
              mutex="Obs - Expected")[met]
    cs   <- if (met=="mutex")
      list(c(0,"#D6604D"),c(0.5,"#0d1117"),c(1,"#2166AC"))
    else
      list(c(0,"#0d1117"),c(0.3,"#1a4a7a"),c(0.6,"#2166AC"),c(1,"#d73027"))
    
    plot_ly(x=vars, y=vars, z=cm, type="heatmap",
            colorscale=cs,
            text=matrix(if(show) as.character(cm) else "",nrow(cm)),
            texttemplate=if(show)"%{text}" else "",
            textfont=list(size=fs,color="#e6edf3"),
            hovertemplate=paste0("<b>%{x}</b> \u00d7 <b>%{y}</b><br>",lbl,": %{z}<extra></extra>"),
            colorbar=list(
              title=list(text=lbl,font=list(color="#8b949e",size=10)),
              tickfont=list(color="#8b949e",size=9),
              bgcolor="#161b22",bordercolor="#30363d",borderwidth=1)
    ) %>% layout(
      paper_bgcolor="#0d1117",plot_bgcolor="#0d1117",
      font=list(color="#e6edf3",family="IBM Plex Mono"),
      xaxis=list(tickfont=list(size=8,color="#8b949e"),tickangle=-45,
                 showgrid=FALSE,zeroline=FALSE),
      yaxis=list(tickfont=list(size=8,color="#8b949e"),
                 showgrid=FALSE,zeroline=FALSE,autorange="reversed"),
      margin=list(l=130,r=20,t=20,b=130)
    ) %>% config(displaylogo=FALSE)
  })
  
  output$comut_freq <- renderPlotly({
    bin  <- mut_bin()
    vars <- top_vars()
    freq <- rowSums(bin[vars,,drop=FALSE]) / ncol(bin) * 100
    ord  <- order(freq,decreasing=TRUE)
    df   <- data.frame(
      variant=factor(vars[ord],levels=vars[ord]),
      pct=freq[ord])
    plot_ly(df, x=~pct, y=~variant, type="bar", orientation="h",
            marker=list(color="#2166AC",line=list(color="#30363d",width=.4)),
            hovertemplate="<b>%{y}</b><br>%{x:.1f}%<extra></extra>"
    ) %>% layout(
      paper_bgcolor="#0d1117",plot_bgcolor="#0d1117",
      font=list(color="#e6edf3",family="IBM Plex Mono"),
      xaxis=list(title="Mutation %",tickfont=list(size=9,color="#8b949e"),
                 gridcolor="#1e2530",zeroline=FALSE),
      yaxis=list(tickfont=list(size=8,color="#8b949e"),
                 autorange="reversed",showgrid=FALSE),
      margin=list(l=140,r=20,t=20,b=50)
    ) %>% config(displaylogo=FALSE)
  })
  
  output$comut_cluster_bar <- renderPlotly({
    mat  <- ngt_mat()
    vars <- top_vars()[seq_len(min(10,length(top_vars())))]
    bdf  <- base_df()
    
    rows <- lapply(vars, function(v) {
      geno <- (mat[v, bdf$cell]==1 | mat[v, bdf$cell]==2)
      data.frame(variant=v, cluster=bdf$cluster, mutated=geno,
                 stringsAsFactors=FALSE)
    })
    df <- do.call(rbind, rows)
    agg <- df %>%
      group_by(variant,cluster) %>%
      summarise(pct=mean(mutated)*100, .groups="drop")
    
    cl_ord <- sort(unique(agg$cluster))
    p <- plot_ly(type="bar")
    for (cl in cl_ord) {
      sub <- agg[agg$cluster==cl,]
      p   <- add_trace(p,
                       x=sub$variant, y=sub$pct, name=paste0("Cl ",cl),
                       marker=list(color=CLUSTER_COLORS[cl]%||%"#aaa"),
                       hovertemplate=paste0("<b>Cluster ",cl,"</b><br>%{x}<br>%{y:.1f}%<extra></extra>"))
    }
    p %>% layout(barmode="group",
                 paper_bgcolor="#0d1117",plot_bgcolor="#0d1117",
                 font=list(color="#e6edf3",family="IBM Plex Mono"),
                 xaxis=list(tickfont=list(size=9,color="#8b949e"),tickangle=-30,
                            showgrid=FALSE,zeroline=FALSE),
                 yaxis=list(title="% mutated",tickfont=list(size=9,color="#8b949e"),
                            gridcolor="#1e2530",zeroline=FALSE),
                 legend=list(bgcolor="#161b22",bordercolor="#30363d",borderwidth=1,
                             font=list(size=10)),
                 margin=list(l=60,r=20,t=10,b=80)
    ) %>% config(displaylogo=FALSE)
  })
  
  # ═══════════════════════════════════════════════════════════
  #  TAB 3: PROTEIN x NGT
  # ═══════════════════════════════════════════════════════════
  
  pxn_df <- reactive({
    bdf  <- base_df()
    prot <- req(input$pxn_protein)
    vari <- req(input$pxn_variant)
    geno <- factor(as.integer(ngt_mat()[vari, bdf$cell]), levels=0:3,
                   labels=c("WT (0)","Het (1)","Hom (2)","Missing (3)"))
    data.frame(cell=bdf$cell, cluster=bdf$cluster,
               protein=as.numeric(dsb_mat()[prot, bdf$cell]),
               genotype=as.character(geno),
               stringsAsFactors=FALSE)
  })
  
  pxn_filt <- reactive({
    df <- pxn_df()
    if (isTRUE(input$pxn_excl_miss)) df <- df[df$genotype!="Missing (3)",]
    df
  })
  
  output$pxn_violin <- renderPlotly({
    df   <- pxn_filt()
    prot <- input$pxn_protein
    vari <- input$pxn_variant
    genos<- intersect(c("WT (0)","Het (1)","Hom (2)"), unique(df$genotype))
    p    <- plot_ly(type="violin")
    for (g in genos) {
      sub <- df[df$genotype==g,]
      col <- NGT_COLORS[g]%||%"#aaa"
      p   <- add_trace(p, y=sub$protein, name=g, type="violin",
                       box=list(visible=TRUE), meanline=list(visible=TRUE),
                       fillcolor=paste0(col,"55"), line=list(color=col),
                       points=if(isTRUE(input$pxn_jitter))"all" else FALSE,
                       marker=list(size=2,opacity=.4,color=col),
                       hovertemplate=paste0("<b>",g,"</b><br>",prot,": %{y:.2f}<extra></extra>"))
    }
    p %>% layout(
      paper_bgcolor="#0d1117",plot_bgcolor="#0d1117",
      font=list(color="#e6edf3",family="IBM Plex Mono"),
      xaxis=list(tickfont=list(size=10,color="#8b949e"),showgrid=FALSE),
      yaxis=list(title=prot,tickfont=list(size=9,color="#8b949e"),
                 gridcolor="#1e2530",zeroline=FALSE),
      legend=list(bgcolor="#161b22",bordercolor="#30363d",borderwidth=1,
                  font=list(size=10)),
      violinmode="group",
      margin=list(l=60,r=20,t=20,b=50)
    ) %>% config(displaylogo=FALSE)
  })
  
  output$pxn_heatmap <- renderPlotly({
    bdf  <- base_df()
    prot <- req(input$pxn_protein)
    freq <- rowMeans(mut_bin(), na.rm=TRUE)
    top20<- names(sort(freq,decreasing=TRUE))[1:min(20,length(freq))]
    pv   <- as.numeric(dsb_mat()[prot, bdf$cell])
    
    rows <- lapply(top20, function(v) {
      geno <- factor(as.integer(ngt_mat()[v, bdf$cell]), levels=0:3,
                     labels=c("WT (0)","Het (1)","Hom (2)","Missing (3)"))
      data.frame(variant=v, genotype=as.character(geno),
                 protein=pv, stringsAsFactors=FALSE)
    })
    df2  <- do.call(rbind, rows)
    df2  <- df2[df2$genotype!="Missing (3)",]
    agg  <- df2 %>%
      group_by(variant,genotype) %>%
      summarise(med=median(protein,na.rm=TRUE),.groups="drop")
    
    gl   <- c("WT (0)","Het (1)","Hom (2)")
    zm   <- matrix(NA,length(top20),3,dimnames=list(top20,gl))
    for (i in seq_len(nrow(agg)))
      zm[agg$variant[i], agg$genotype[i]] <- agg$med[i]
    
    plot_ly(x=gl, y=top20, z=zm, type="heatmap",
            colorscale=list(c(0,"#0d1117"),c(0.5,"#1a4a7a"),c(1,"#d73027")),
            hovertemplate="<b>%{y}</b><br>%{x}<br>Median: %{z:.2f}<extra></extra>",
            colorbar=list(
              title=list(text=paste("Median",prot),font=list(color="#8b949e",size=10)),
              tickfont=list(color="#8b949e",size=9),
              bgcolor="#161b22",bordercolor="#30363d",borderwidth=1)
    ) %>% layout(
      paper_bgcolor="#0d1117",plot_bgcolor="#0d1117",
      font=list(color="#e6edf3",family="IBM Plex Mono"),
      xaxis=list(tickfont=list(size=10,color="#8b949e"),showgrid=FALSE),
      yaxis=list(tickfont=list(size=8,color="#8b949e"),
                 autorange="reversed",showgrid=FALSE),
      margin=list(l=140,r=20,t=20,b=60)
    ) %>% config(displaylogo=FALSE)
  })
  
  output$pxn_scatter <- renderPlotly({
    bdf    <- base_df()
    prot2  <- req(input$pxn_prot2)
    prot3  <- req(input$pxn_prot3)
    vari   <- input$pxn_variant
    by     <- input$pxn_color_by%||%"ngt"
    xv     <- as.numeric(dsb_mat()[prot2, bdf$cell])
    yv     <- as.numeric(dsb_mat()[prot3, bdf$cell])
    if (by=="ngt" && !is.null(vari) && nchar(vari)>0) {
      geno <- factor(as.integer(ngt_mat()[vari, bdf$cell]), levels=0:3,
                     labels=c("WT (0)","Het (1)","Hom (2)","Missing (3)"))
      cv <- as.character(geno); pal <- NGT_COLORS; ttl <- vari
    } else {
      cv <- bdf$cluster; pal <- CLUSTER_COLORS; ttl <- "Cluster"
    }
    p <- plot_ly(type="scatter",mode="markers")
    for (cat in unique(cv)) {
      idx <- cv==cat
      p   <- add_trace(p, x=xv[idx], y=yv[idx], name=cat,
                       marker=list(color=pal[cat]%||%"#aaa",size=3,opacity=.6),
                       hovertemplate=paste0("<b>",prot2,":</b> %{x:.2f}<br>",
                                            "<b>",prot3,":</b> %{y:.2f}<extra></extra>"))
    }
    p %>% layout(
      paper_bgcolor="#0d1117",plot_bgcolor="#0d1117",
      font=list(color="#e6edf3",family="IBM Plex Mono"),
      xaxis=list(title=prot2,tickfont=list(size=9,color="#8b949e"),
                 gridcolor="#1e2530",zeroline=FALSE),
      yaxis=list(title=prot3,tickfont=list(size=9,color="#8b949e"),
                 gridcolor="#1e2530",zeroline=FALSE),
      legend=list(title=list(text=ttl),bgcolor="#161b22",bordercolor="#30363d",
                  borderwidth=1,font=list(size=10)),
      margin=list(l=60,r=20,t=20,b=50)
    ) %>% config(displaylogo=FALSE)
  })
  
  output$pxn_mutload <- renderPlotly({
    bdf  <- base_df()
    prot <- req(input$pxn_protein)
    pv   <- as.numeric(dsb_mat()[prot, bdf$cell])
    load <- colSums(mut_bin()[, bdf$cell, drop=FALSE])
    df   <- data.frame(load_bin=factor(pmin(load,8)), protein=pv)
    plot_ly(df, x=~load_bin, y=~protein, type="box",
            color=~load_bin, colors=viridis(9),
            hovertemplate="Load %{x}<br>%{y:.2f}<extra></extra>",
            showlegend=FALSE
    ) %>% layout(
      paper_bgcolor="#0d1117",plot_bgcolor="#0d1117",
      font=list(color="#e6edf3",family="IBM Plex Mono"),
      xaxis=list(title="Mutation Load (# variants)",
                 tickfont=list(size=9,color="#8b949e"),showgrid=FALSE,zeroline=FALSE),
      yaxis=list(title=prot,tickfont=list(size=9,color="#8b949e"),
                 gridcolor="#1e2530",zeroline=FALSE),
      margin=list(l=60,r=20,t=20,b=50)
    ) %>% config(displaylogo=FALSE)
  })
  
  output$pxn_stats <- renderText({
    df   <- pxn_filt()
    prot <- input$pxn_protein
    vari <- input$pxn_variant
    if (nrow(df)==0) return("No data")
    tbl  <- table(df$genotype)
    ln   <- mapply(function(nm,ct) sprintf("  %-10s %d",nm,ct),
                   names(tbl),as.integer(tbl),SIMPLIFY=FALSE)
    wm   <- median(df$protein[df$genotype=="WT (0)"],  na.rm=TRUE)
    hm   <- median(df$protein[df$genotype=="Het (1)"], na.rm=TRUE)
    hom  <- median(df$protein[df$genotype=="Hom (2)"], na.rm=TRUE)
    paste0("Protein: ",prot,"\nVariant: ",vari,"\n\n",
           "Genotype N:\n",paste(ln,collapse="\n"),
           "\n\nMedian expression:\n",
           sprintf("  WT:  %.2f\n  Het: %.2f\n  Hom: %.2f",
                   wm%||%NA, hm%||%NA, hom%||%NA))
  })
}

shinyApp(ui, server)