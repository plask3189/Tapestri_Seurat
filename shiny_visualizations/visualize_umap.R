## ─────────────────────────────────────────────────────────────
##  Multi-Modal Explorer  ·  DSB_norm × NGT  (multi-sample)
##  Handles Seurat5 split layers (counts.1, counts.2 …)
##  Object name: seurat_merged  (or seurat_obj as fallback)
## ─────────────────────────────────────────────────────────────

library(shiny)
library(plotly)
library(Seurat)
library(SeuratObject)
library(dplyr)
library(tidyr)
library(viridis)
library(Matrix)

# ═══════════════════════════════════════════════════════════════
#  HELPERS
# ═══════════════════════════════════════════════════════════════

`%||%` <- function(a, b) if (!is.null(a)) a else b

NGT_COLORS <- c(
  "WT (0)"      = "#2166AC",
  "Het (1)"     = "#F4A582",
  "Hom (2)"     = "#D6604D",
  "Missing (3)" = "#888888",
  "N/A"         = "#444444"
)

CLUSTER_COLORS <- c(
  "0"="#E63946","1"="#2A9D8F","2"="#E9C46A","3"="#F4A261",
  "4"="#457B9D","5"="#A8DADC","6"="#6A4C93","7"="#FF9F1C",
  "8"="#2EC4B6","9"="#CBDF90"
)

get_obj <- function() {
  for (nm in c("seurat_merged","seurat_obj")) {
    if (exists(nm, envir=.GlobalEnv))
      return(get(nm, envir=.GlobalEnv))
  }
  p <- c("seurat_merged.rds","seurat_obj.rds")
  for (pp in p) if (file.exists(pp)) return(readRDS(pp))
  stop("No Seurat object found. Load seurat_merged or seurat_obj into the global env.")
}

join_layers <- function(obj, assay, layer_prefix = "counts") {
  av <- Assays(obj)
  if (!assay %in% av) stop("Assay '", assay, "' not found.")
  assay_obj <- obj@assays[[assay]]
  all_layers <- names(assay_obj@layers)
  target     <- grep(paste0("^", layer_prefix), all_layers, value=TRUE)
  if (length(target) == 0)
    stop("No layers matching '", layer_prefix, "' in assay '", assay, "'.")
  if (length(target) == 1) {
    mat <- as.matrix(LayerData(obj, assay=assay, layer=target[1]))
    if (is.null(rownames(mat))) rownames(mat) <- rownames(assay_obj@features@.Data)
    if (is.null(colnames(mat))) colnames(mat) <- rownames(assay_obj@cells@.Data)
    return(mat)
  }
  pieces <- lapply(target, function(lyr) {
    m <- as.matrix(LayerData(obj, assay=assay, layer=lyr))
    feat_logi  <- assay_obj@features@.Data[, lyr, drop=FALSE]
    cell_logi  <- assay_obj@cells@.Data[,    lyr, drop=FALSE]
    feat_names <- rownames(feat_logi)[feat_logi[, 1]]
    cell_names <- rownames(cell_logi)[cell_logi[, 1]]
    if (nrow(m) == length(feat_names)) rownames(m) <- feat_names
    if (ncol(m) == length(cell_names)) colnames(m) <- cell_names
    m
  })
  all_feats <- unique(unlist(lapply(pieces, rownames)))
  all_cells <- unique(unlist(lapply(pieces, colnames)))
  combined  <- matrix(NA_real_, nrow=length(all_feats), ncol=length(all_cells),
                      dimnames=list(all_feats, all_cells))
  for (m in pieces) combined[rownames(m), colnames(m)] <- m[rownames(m), colnames(m)]
  combined
}

build_base_df <- function(obj) {
  umap_2d <- Embeddings(obj, "umap")
  pca_emb <- Embeddings(obj, "pca")
  cells   <- rownames(umap_2d)
  meta    <- obj@meta.data[cells, , drop=FALSE]
  sample_col <- if ("sample_name" %in% colnames(meta)) meta$sample_name
  else sub("_[^_]+$", "", cells)
  data.frame(
    cell    = cells,
    umap_1  = umap_2d[, 1],
    umap_2  = umap_2d[, 2],
    umap_3  = if (ncol(umap_2d) >= 3) umap_2d[, 3] else pca_emb[, "PC_3"],
    cluster = as.character(meta$seurat_clusters),
    sample  = sample_col,
    stringsAsFactors = FALSE
  )
}

add_dsb_col <- function(df, mat, feat, clamp = NULL) {
  vals <- as.numeric(mat[feat, df$cell])
  if (!is.null(clamp) && length(clamp) == 2)
    vals <- pmax(pmin(vals, clamp[2]), clamp[1])
  df$color_val   <- vals
  df$color_type  <- "continuous"
  df$color_label <- feat
  df
}

add_ngt_col <- function(df, mat, feat) {
  if (feat %in% rownames(mat)) {
    raw <- mat[feat, df$cell]
    lv  <- ifelse(is.na(raw), "N/A",
                  ifelse(raw == 0, "WT (0)",
                         ifelse(raw == 1, "Het (1)",
                                ifelse(raw == 2, "Hom (2)", "Missing (3)"))))
  } else {
    lv <- rep("N/A", nrow(df))
  }
  df$color_val   <- lv
  df$color_type  <- "categorical"
  df$color_label <- feat
  df
}

make_axis_style <- function(ax_vis, gd_vis, title_str) {
  list(
    showline=ax_vis, showticklabels=ax_vis,
    title=list(text=title_str, font=list(color="#8b949e", size=11)),
    tickfont=list(color="#8b949e", size=9),
    gridcolor=if(gd_vis)"#30363d" else "transparent",
    zerolinecolor=if(gd_vis)"#30363d" else "transparent",
    backgroundcolor="#0d1117", showbackground=TRUE
  )
}

scatter3d_continuous <- function(df, sz, alpha, cmin=NULL, cmax=NULL) {
  plot_ly(data=df, x=~umap_1, y=~umap_2, z=~umap_3,
          type="scatter3d", mode="markers",
          marker=list(
            size=sz, opacity=alpha, color=~color_val,
            cmin=cmin, cmax=cmax,
            colorscale=list(c(0,"#0d1117"),c(0.15,"#1a4a7a"),c(0.4,"#2166AC"),
                            c(0.65,"#4dac26"),c(0.85,"#f4a40e"),c(1,"#d73027")),
            colorbar=list(
              title=list(text=df$color_label[1], font=list(color="#e6edf3",size=12)),
              tickfont=list(color="#8b949e",size=10),
              bgcolor="#161b22",bordercolor="#30363d",borderwidth=1,len=0.5),
            showscale=TRUE),
          text=~paste0("<b>Cell:</b> ",cell,"<br><b>Cluster:</b> ",cluster,
                       "<br><b>Sample:</b> ",sample,
                       "<br><b>",color_label,":</b> ",round(color_val,2)),
          hoverinfo="text")
}

scatter3d_categorical <- function(df, sz, alpha, pal) {
  lvl_order <- c("WT (0)","Het (1)","Hom (2)","Missing (3)","N/A")
  cats <- c(intersect(lvl_order, unique(df$color_val)),
            setdiff(unique(df$color_val), lvl_order))
  p    <- plot_ly(type="scatter3d", mode="markers")
  for (cat in cats) {
    sub <- df[df$color_val == cat, ]
    if (nrow(sub) == 0) next
    p   <- add_trace(p,
                     x=sub$umap_1, y=sub$umap_2, z=sub$umap_3, name=cat,
                     marker=list(color=unname(pal[cat] %||% "#aaaaaa"), size=sz, opacity=alpha),
                     text=paste0("<b>Cell:</b> ",sub$cell,
                                 "<br><b>Cluster:</b> ",sub$cluster,
                                 "<br><b>Sample:</b> ",sub$sample,
                                 "<br><b>",sub$color_label,":</b> ",cat),
                     hoverinfo="text")
  }
  p
}

apply_3d_layout <- function(p, ax_vis, gd_vis, z_label="PC 3", title=NULL) {
  p %>%
    layout(
      title=if(!is.null(title)) list(text=title,font=list(color="#8b949e",size=12)) else NULL,
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
      margin=list(l=0,r=0,t=if(!is.null(title))30 else 10,b=0)
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
  --acc:#58a6ff; --text:#e6edf3; --muted:#8b949e;
}
*{box-sizing:border-box;margin:0;padding:0;}
html,body{height:100%;background:var(--bg);color:var(--text);
  font-family:'IBM Plex Sans',sans-serif;font-size:13px;overflow:hidden;}

.topbar{display:flex;align-items:center;gap:14px;padding:0 20px;
  background:var(--surf);border-bottom:1px solid var(--border);height:50px;flex-shrink:0;}
.logo{font-family:'IBM Plex Mono',monospace;font-size:14px;font-weight:600;
  color:var(--acc);letter-spacing:.04em;}
.subtitle{color:var(--muted);font-size:11px;font-family:'IBM Plex Mono',monospace;}
.topbar-spacer{flex:1;}
.pill{background:var(--bg);border:1px solid var(--border);border-radius:20px;
  padding:2px 10px;font-family:'IBM Plex Mono',monospace;font-size:11px;color:var(--muted);}
.pill span{color:var(--acc);}

.nav-tabs{border-bottom:1px solid var(--border)!important;
  background:var(--surf);padding:0 20px;flex-shrink:0;}
.nav-tabs>li>a{font-family:'IBM Plex Mono',monospace!important;font-size:11px!important;
  letter-spacing:.06em!important;text-transform:uppercase!important;
  color:var(--muted)!important;border:none!important;
  border-bottom:2px solid transparent!important;background:transparent!important;
  padding:10px 16px!important;margin:0!important;border-radius:0!important;
  transition:color .15s,border-color .15s;}
.nav-tabs>li.active>a,.nav-tabs>li>a:hover{
  color:var(--acc)!important;border-bottom-color:var(--acc)!important;
  background:transparent!important;}
.tab-content{height:calc(100vh - 96px);overflow:hidden;}
.tab-pane{height:100%;}

.dual-layout{display:flex;height:100%;overflow:hidden;}
.dual-sidebar{width:260px;min-width:260px;background:var(--surf);
  border-right:1px solid var(--border);padding:14px;overflow-y:auto;
  display:flex;flex-direction:column;gap:12px;}
.dual-plots{flex:1;display:flex;overflow:hidden;}
.plot-half{flex:1;border-right:1px solid var(--border);position:relative;}
.plot-half:last-child{border-right:none;}
.plot-title{position:absolute;top:8px;left:50%;transform:translateX(-50%);
  font-family:'IBM Plex Mono',monospace;font-size:11px;letter-spacing:.06em;
  text-transform:uppercase;color:var(--muted);z-index:10;pointer-events:none;
  white-space:nowrap;}

.analytics-layout{display:flex;height:100%;overflow:hidden;}
.analytics-sidebar{width:260px;min-width:260px;background:var(--surf);
  border-right:1px solid var(--border);padding:14px;overflow-y:auto;
  display:flex;flex-direction:column;gap:12px;}
.analytics-main{flex:1;overflow-y:auto;padding:14px;
  display:flex;flex-direction:column;gap:14px;}
.analytics-row{display:flex;gap:14px;min-height:420px;}
.analytics-card{flex:1;background:var(--surf);border:1px solid var(--border);
  border-radius:8px;padding:12px;display:flex;flex-direction:column;min-width:0;}
.card-title{font-family:'IBM Plex Mono',monospace;font-size:10px;
  letter-spacing:.1em;text-transform:uppercase;color:var(--muted);margin-bottom:8px;}

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

/* color scale section */
.scale-row{display:flex;align-items:center;gap:6px;margin-bottom:4px;}
.scale-label{font-family:'IBM Plex Mono',monospace;font-size:10px;
  color:var(--muted);width:28px;flex-shrink:0;}
.scale-bar{flex:1;height:8px;border-radius:4px;
  background:linear-gradient(to right,#0d1117,#1a4a7a,#2166AC,#4dac26,#f4a40e,#d73027);
  border:1px solid var(--border);}

.ngt-legend{display:flex;gap:6px;flex-wrap:wrap;margin-top:6px;}
.ngt-sw{display:flex;align-items:center;gap:4px;
  font-family:'IBM Plex Mono',monospace;font-size:10px;color:var(--text);}
.sw-dot{width:9px;height:9px;border-radius:50%;flex-shrink:0;}

.var-tag{display:inline-block;font-family:'IBM Plex Mono',monospace;
  font-size:9px;padding:1px 6px;border-radius:3px;margin-left:4px;
  background:#1e3a5f;color:#58a6ff;vertical-align:middle;}
.var-tag.partial{background:#3a2a00;color:#d29922;}
.var-tag.na{background:#2a1a1a;color:#f85149;}

.shiny-text-output{background:var(--bg)!important;color:var(--muted)!important;
  font-family:'IBM Plex Mono',monospace!important;font-size:10px!important;
  border:1px solid var(--border)!important;border-radius:4px!important;
  padding:7px!important;white-space:pre-wrap!important;max-height:200px;overflow-y:auto;}

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
  
  div(class="topbar",
      div(class="logo","Multi-Modal Explorer"),
      div(class="subtitle","DSB\u2009\u00d7\u2009NGT"),
      div(class="topbar-spacer"),
      uiOutput("topbar_pills")
  ),
  
  tabsetPanel(id="main_tabs",
              
              # ── TAB 1: Side-by-Side UMAP ───────────────────────────────
              tabPanel("Side-by-Side UMAP",
                       div(class="dual-layout",
                           div(class="dual-sidebar",
                               
                               div(class="cs",
                                   div(class="cs-title","Left — DSB Protein"),
                                   selectInput("left_feat","Marker",choices=character(0))
                               ),
                               
                               # ── DSB Color Scale ──────────────────────────────────
                               div(class="cs",
                                   div(class="cs-title","DSB Color Scale"),
                                   div(class="scale-row",
                                       div(class="scale-label","low"),
                                       div(class="scale-bar"),
                                       div(class="scale-label",style="text-align:right;","high")
                                   ),
                                   uiOutput("dsb_range_ui"),
                                   actionButton("dsb_reset","Reset to data range",
                                                style="width:100%;margin-top:4px;background:#161b22;
                                  border:1px solid #30363d;color:#8b949e;
                                  font-family:'IBM Plex Mono',monospace;font-size:10px;
                                  padding:3px 0;border-radius:4px;")
                               ),
                               
                               div(class="cs",
                                   div(class="cs-title","Right — NGT Variant"),
                                   selectInput("right_feat","Variant",choices=character(0)),
                                   uiOutput("right_feat_availability"),
                                   div(class="ngt-legend",
                                       div(class="ngt-sw",div(class="sw-dot",style="background:#2166AC;"),"WT (0)"),
                                       div(class="ngt-sw",div(class="sw-dot",style="background:#F4A582;"),"Het (1)"),
                                       div(class="ngt-sw",div(class="sw-dot",style="background:#D6604D;"),"Hom (2)"),
                                       div(class="ngt-sw",div(class="sw-dot",style="background:#888888;"),"Miss (3)"),
                                       div(class="ngt-sw",div(class="sw-dot",style="background:#444444;"),"N/A")
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
                                   div(class="cs-title","Filter by Sample"),
                                   uiOutput("sb_sample_ui")
                               ),
                               div(class="cs",
                                   div(class="cs-title","Filter by Cluster"),
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
              
              # ── TAB 2: Co-Mutation ─────────────────────────────────────
              tabPanel("Co-Mutation",
                       div(class="analytics-layout",
                           div(class="analytics-sidebar",
                               div(class="cs",
                                   div(class="cs-title","Metric"),
                                   selectInput("comut_metric","",
                                               choices=c("Co-occurring cells"="co_occur",
                                                         "Jaccard similarity"="jaccard",
                                                         "Obs - Expected (mutex)"="mutex"),
                                               selected="co_occur")
                               ),
                               div(class="cs",
                                   div(class="cs-title","Sample Scope"),
                                   uiOutput("comut_sample_ui"),
                                   helpText(style="color:#8b949e;font-size:10px;font-family:'IBM Plex Mono',monospace;",
                                            "Variants not in a sample's panel are treated as N/A.")
                               ),
                               div(class="cs",
                                   div(class="cs-title","Variant Filter"),
                                   numericInput("comut_vaf_min","Min VAF %",value=0,min=0,max=100,step=1),
                                   numericInput("comut_n_top","Top N variants",value=30,min=5,max=200,step=5),
                                   checkboxInput("comut_excl_missing","Exclude Missing (3) from mut",value=TRUE)
                               ),
                               div(class="cs",
                                   div(class="cs-title","Labels"),
                                   checkboxInput("comut_showvals","Show values",value=FALSE),
                                   sliderInput("comut_fontsize","Label size",min=6,max=14,value=8,step=1)
                               )
                           ),
                           div(class="analytics-main",
                               div(class="analytics-row",
                                   div(class="analytics-card",style="flex:2;",
                                       div(class="card-title","Co-Mutation Heatmap"),
                                       plotlyOutput("comut_heatmap",width="100%",height="390px")
                                   ),
                                   div(class="analytics-card",style="flex:1;",
                                       div(class="card-title","Mutation Frequency by Sample"),
                                       plotlyOutput("comut_freq",width="100%",height="390px")
                                   )
                               ),
                               div(class="analytics-row",style="min-height:340px;",
                                   div(class="analytics-card",
                                       div(class="card-title","% Mutated per Cluster — Top 10 Variants"),
                                       plotlyOutput("comut_cluster_bar",width="100%",height="290px")
                                   )
                               )
                           )
                       )
              ),
              
              # ── TAB 3: Protein × NGT ───────────────────────────────────
              tabPanel("Protein \u00d7 NGT",
                       div(class="analytics-layout",
                           div(class="analytics-sidebar",
                               div(class="cs",
                                   div(class="cs-title","Primary Features"),
                                   selectInput("pxn_protein","Protein",choices=character(0)),
                                   selectInput("pxn_variant","Variant",choices=character(0)),
                                   uiOutput("pxn_variant_avail")
                               ),
                               div(class="cs",
                                   div(class="cs-title","Sample Filter"),
                                   uiOutput("pxn_sample_ui")
                               ),
                               div(class="cs",
                                   div(class="cs-title","Violin Options"),
                                   checkboxInput("pxn_jitter","Show jitter",value=FALSE),
                                   checkboxInput("pxn_excl_miss","Exclude Missing (3)",value=TRUE),
                                   checkboxInput("pxn_excl_na","Exclude N/A (panel mismatch)",value=TRUE)
                               ),
                               div(class="cs",
                                   div(class="cs-title","Scatter Options"),
                                   selectInput("pxn_prot2","Protein X",choices=character(0)),
                                   selectInput("pxn_prot3","Protein Y",choices=character(0)),
                                   selectInput("pxn_color_by","Color by",
                                               choices=c("NGT Genotype"="ngt","Cluster"="cluster","Sample"="sample"))
                               ),
                               div(class="cs",
                                   div(class="cs-title","Stats"),
                                   verbatimTextOutput("pxn_stats")
                               )
                           ),
                           div(class="analytics-main",
                               div(class="analytics-row",
                                   div(class="analytics-card",style="flex:1;",
                                       div(class="card-title","Protein by Genotype"),
                                       plotlyOutput("pxn_violin",width="100%",height="390px")
                                   ),
                                   div(class="analytics-card",style="flex:1;",
                                       div(class="card-title","Median Protein \u00d7 Top-20 Variants \u00d7 Genotype"),
                                       plotlyOutput("pxn_heatmap",width="100%",height="390px")
                                   )
                               ),
                               div(class="analytics-row",style="min-height:340px;",
                                   div(class="analytics-card",style="flex:1;",
                                       div(class="card-title","Protein Scatter — NGT / Cluster / Sample"),
                                       plotlyOutput("pxn_scatter",width="100%",height="300px")
                                   ),
                                   div(class="analytics-card",style="flex:1;",
                                       div(class="card-title","Mutation Load vs Protein"),
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
  
  # ── Core data ─────────────────────────────────────────────────
  obj <- reactive({ get_obj() })
  
  base_df <- reactive({ build_base_df(obj()) })
  
  dsb_mat <- reactive({
    withProgress(message="Joining DSB layers…", value=0.5,
                 join_layers(obj(), assay="DSB_norm", layer_prefix="counts"))
  })
  
  ngt_mat <- reactive({
    withProgress(message="Joining NGT layers…", value=0.5,
                 join_layers(obj(), assay="NGT", layer_prefix="counts"))
  })
  
  mut_bin <- reactive({
    m <- ngt_mat()
    m == 1 | m == 2
  })
  
  var_sample_map <- reactive({
    bdf   <- base_df(); mat <- ngt_mat()
    samps <- sort(unique(bdf$sample))
    lapply(setNames(samps, samps), function(s) {
      cells_s <- intersect(bdf$cell[bdf$sample == s], colnames(mat))
      apply(!is.na(mat[, cells_s, drop=FALSE]), 1, any)
    })
  })
  
  var_avail <- reactive({
    vsm    <- var_sample_map()
    n_samp <- length(vsm)
    vars   <- rownames(ngt_mat())
    sapply(vars, function(v) {
      n <- sum(sapply(vsm, function(x) isTRUE(x[v])))
      if (n == n_samp) "all" else if (n == 0) "none" else "partial"
    })
  })
  
  # ── DSB color scale ───────────────────────────────────────────
  # Compute per-feature 1st/99th pctile range
  dsb_feat_range <- reactive({
    feat <- req(input$left_feat)
    mat  <- dsb_mat()
    req(feat %in% rownames(mat))
    vals <- as.numeric(mat[feat, ])
    list(
      lo = floor(quantile(vals, 0.01, na.rm=TRUE) * 10) / 10,
      hi = ceiling(quantile(vals, 0.99, na.rm=TRUE) * 10) / 10
    )
  })
  
  # Track whether user has manually moved the slider
  dsb_user_range <- reactiveVal(NULL)
  
  # Reset user range when feature changes or reset button pressed
  observeEvent(input$left_feat, { dsb_user_range(NULL) })
  observeEvent(input$dsb_reset,  { dsb_user_range(NULL) })
  
  # Capture manual slider moves
  observeEvent(input$dsb_range, {
    rng <- dsb_feat_range()
    # Only flag as user-changed if different from computed defaults
    if (!is.null(input$dsb_range) &&
        (abs(input$dsb_range[1] - rng$lo) > 0.05 ||
         abs(input$dsb_range[2] - rng$hi) > 0.05)) {
      dsb_user_range(input$dsb_range)
    }
  }, ignoreInit=TRUE)
  
  output$dsb_range_ui <- renderUI({
    rng <- dsb_feat_range()
    # Use user-set value if available, else data-driven default
    cur <- dsb_user_range() %||% c(rng$lo, rng$hi)
    # Expand slider bounds slightly beyond data range
    bound_lo <- min(rng$lo, cur[1]) - 1
    bound_hi <- max(rng$hi, cur[2]) + 1
    sliderInput("dsb_range", NULL,
                min=bound_lo, max=bound_hi,
                value=cur, step=0.1,
                ticks=FALSE)
  })
  
  # Effective clamp values for plotting
  dsb_clamp <- reactive({
    input$dsb_range %||% {
      rng <- dsb_feat_range()
      c(rng$lo, rng$hi)
    }
  })
  
  # ── Topbar pills ──────────────────────────────────────────────
  output$topbar_pills <- renderUI({
    bdf <- base_df()
    tagList(
      div(class="pill","cells: ",    tags$span(format(nrow(bdf),big.mark=","))),
      div(class="pill","samples: ",  tags$span(length(unique(bdf$sample)))),
      div(class="pill","clusters: ", tags$span(length(unique(bdf$cluster)))),
      div(class="pill","proteins: ", tags$span(nrow(dsb_mat()))),
      div(class="pill","variants: ", tags$span(nrow(ngt_mat())))
    )
  })
  
  # ── Populate selectors ────────────────────────────────────────
  observe({
    dsb <- sort(rownames(dsb_mat()))
    ngt <- sort(rownames(ngt_mat()))
    updateSelectInput(session,"left_feat",   choices=dsb, selected=dsb[1])
    updateSelectInput(session,"right_feat",  choices=ngt, selected=ngt[1])
    updateSelectInput(session,"pxn_protein", choices=dsb, selected=dsb[1])
    updateSelectInput(session,"pxn_variant", choices=ngt, selected=ngt[1])
    updateSelectInput(session,"pxn_prot2",   choices=dsb, selected=dsb[1])
    updateSelectInput(session,"pxn_prot3",   choices=dsb, selected=dsb[min(2,length(dsb))])
  })
  
  # ── Availability badges ───────────────────────────────────────
  avail_badge <- function(feat, av) {
    if (is.null(feat) || !feat %in% names(av)) return(NULL)
    switch(av[feat],
           all     = tags$span(class="var-tag","all samples"),
           partial = tags$span(class="var-tag partial","partial coverage"),
           none    = tags$span(class="var-tag na","no coverage")
    )
  }
  
  output$right_feat_availability <- renderUI({
    avail_badge(input$right_feat, var_avail())
  })
  output$pxn_variant_avail <- renderUI({
    avail_badge(input$pxn_variant, var_avail())
  })
  
  # ── Filter UIs ────────────────────────────────────────────────
  make_sample_ui <- function(id) renderUI({
    sams <- sort(unique(base_df()$sample))
    checkboxGroupInput(id, NULL, choices=sams, selected=sams)
  })
  output$sb_sample_ui    <- make_sample_ui("sb_samples")
  output$comut_sample_ui <- make_sample_ui("comut_samples")
  output$pxn_sample_ui   <- make_sample_ui("pxn_samples")
  
  output$sb_cluster_ui <- renderUI({
    cls <- sort(unique(base_df()$cluster))
    checkboxGroupInput("sb_clusters", NULL, choices=cls, selected=cls)
  })
  
  filt_df <- reactive({
    df <- base_df()
    if (!is.null(input$sb_samples))  df <- df[df$sample  %in% input$sb_samples, ]
    if (!is.null(input$sb_clusters)) df <- df[df$cluster %in% input$sb_clusters, ]
    df
  })
  
  z_label <- reactive({
    if (ncol(Embeddings(obj(),"umap")) >= 3) "UMAP 3" else "PC 3"
  })
  
  # ── LEFT PLOT (DSB) ───────────────────────────────────────────
  output$left_plot <- renderPlotly({
    df   <- filt_df()
    feat <- req(input$left_feat)
    mat  <- dsb_mat()
    req(feat %in% rownames(mat))
    clamp <- dsb_clamp()
    df    <- add_dsb_col(df, mat, feat, clamp=clamp)
    scatter3d_continuous(df, input$sb_sz %||% 3, input$sb_alpha %||% .65,
                         cmin=clamp[1], cmax=clamp[2]) %>%
      apply_3d_layout(isTRUE(input$sb_axes), isTRUE(input$sb_grid), z_label())
  })
  
  # ── RIGHT PLOT (NGT) ──────────────────────────────────────────
  output$right_plot <- renderPlotly({
    df   <- filt_df()
    feat <- req(input$right_feat)
    df   <- add_ngt_col(df, ngt_mat(), feat)
    scatter3d_categorical(df, input$sb_sz %||% 3, input$sb_alpha %||% .65, NGT_COLORS) %>%
      apply_3d_layout(isTRUE(input$sb_axes), isTRUE(input$sb_grid), z_label())
  })
  
  # ── Side-by-side stats ────────────────────────────────────────
  output$sb_stats <- renderText({
    df   <- filt_df()
    feat <- input$right_feat %||% ""
    samp_lines <- paste(mapply(function(s,n) sprintf("  %-12s %d", s, n),
                               names(table(df$sample)),
                               as.integer(table(df$sample))), collapse="\n")
    if (nchar(feat) > 0 && feat %in% rownames(ngt_mat())) {
      raw  <- ngt_mat()[feat, df$cell]
      lv   <- ifelse(is.na(raw),"N/A",
                     ifelse(raw==0,"WT(0)",
                            ifelse(raw==1,"Het(1)",
                                   ifelse(raw==2,"Hom(2)","Miss(3)"))))
      tbl  <- table(factor(lv, levels=c("WT(0)","Het(1)","Hom(2)","Miss(3)","N/A")))
      pct  <- round(100*tbl/sum(tbl),1)
      geno_lines <- paste(mapply(function(nm,ct,pc)
        sprintf("  %-8s %4d  %5.1f%%", nm, ct, pc),
        names(tbl), as.integer(tbl), pct), collapse="\n")
      paste0("Cells: ",nrow(df),"\n\nBy sample:\n",samp_lines,
             "\n\nVariant: ",feat,"\n",geno_lines)
    } else {
      paste0("Cells: ",nrow(df),"\n\nBy sample:\n",samp_lines)
    }
  })
  
  # ═══════════════════════════════════════════════════════════
  #  TAB 2: CO-MUTATION
  # ═══════════════════════════════════════════════════════════
  
  comut_cells <- reactive({
    bdf  <- base_df()
    sels <- input$comut_samples %||% unique(bdf$sample)
    intersect(bdf$cell[bdf$sample %in% sels], colnames(ngt_mat()))
  })
  
  top_vars <- reactive({
    cells  <- comut_cells()
    mat    <- ngt_mat()
    vmeta  <- obj()@assays$NGT@meta.data
    vaf_mn <- input$comut_vaf_min %||% 0
    n_top  <- input$comut_n_top   %||% 30
    excl   <- isTRUE(input$comut_excl_missing)
    sub    <- mat[, cells, drop=FALSE]
    mut_fn <- if (excl) function(x) x==1|x==2 else function(x) x>=1&x<=3
    freq   <- rowMeans(mut_fn(sub), na.rm=TRUE)
    if ("VAF" %in% colnames(vmeta)) {
      vaf_vals <- setNames(vmeta$VAF, rownames(vmeta))
      freq <- freq[is.na(vaf_vals[names(freq)]) | vaf_vals[names(freq)] >= vaf_mn]
    }
    names(sort(freq, decreasing=TRUE))[seq_len(min(n_top, length(freq)))]
  })
  
  comut_mat_data <- reactive({
    cells <- comut_cells()
    mat   <- ngt_mat()
    vars  <- top_vars()
    excl  <- isTRUE(input$comut_excl_missing)
    met   <- input$comut_metric %||% "co_occur"
    sub   <- mat[vars, cells, drop=FALSE]
    bin   <- if (excl) sub==1|sub==2 else sub==1|sub==2|sub==3
    bin[is.na(bin)] <- FALSE
    nc <- length(cells); n <- length(vars)
    m  <- matrix(0, n, n, dimnames=list(vars,vars))
    for (i in seq_len(n)) for (j in seq_len(n)) {
      a <- bin[i,]; b <- bin[j,]
      m[i,j] <- switch(met,
                       co_occur = sum(a & b),
                       jaccard  = { u <- sum(a|b); if(u==0) 0 else round(sum(a&b)/u,3) },
                       mutex    = round(sum(a&b) - sum(a)*sum(b)/nc, 2)
      )
    }
    m
  })
  
  output$comut_heatmap <- renderPlotly({
    cm   <- comut_mat_data()
    vars <- rownames(cm)
    show <- isTRUE(input$comut_showvals)
    fs   <- input$comut_fontsize %||% 8
    met  <- input$comut_metric   %||% "co_occur"
    lbl  <- c(co_occur="Co-occurring cells",jaccard="Jaccard",mutex="Obs-Exp")[met]
    cs   <- if (met=="mutex")
      list(c(0,"#D6604D"),c(0.5,"#0d1117"),c(1,"#2166AC"))
    else
      list(c(0,"#0d1117"),c(0.3,"#1a4a7a"),c(0.6,"#2166AC"),c(1,"#d73027"))
    plot_ly(x=vars, y=vars, z=cm, type="heatmap",
            colorscale=cs,
            text=matrix(if(show) as.character(cm) else "", nrow(cm)),
            texttemplate=if(show) "%{text}" else "",
            textfont=list(size=fs,color="#e6edf3"),
            hovertemplate=paste0("<b>%{x}</b> × <b>%{y}</b><br>",lbl,": %{z}<extra></extra>"),
            colorbar=list(
              title=list(text=lbl,font=list(color="#8b949e",size=10)),
              tickfont=list(color="#8b949e",size=9),
              bgcolor="#161b22",bordercolor="#30363d",borderwidth=1)
    ) %>% layout(
      paper_bgcolor="#0d1117",plot_bgcolor="#0d1117",
      font=list(color="#e6edf3",family="IBM Plex Mono"),
      xaxis=list(tickfont=list(size=8,color="#8b949e"),tickangle=-45,showgrid=FALSE,zeroline=FALSE),
      yaxis=list(tickfont=list(size=8,color="#8b949e"),showgrid=FALSE,zeroline=FALSE,autorange="reversed"),
      margin=list(l=130,r=20,t=20,b=130)
    ) %>% config(displaylogo=FALSE)
  })
  
  output$comut_freq <- renderPlotly({
    bdf   <- base_df()
    sels  <- input$comut_samples %||% unique(bdf$sample)
    cells <- comut_cells()
    mat   <- ngt_mat()
    vars  <- top_vars()
    samps <- sort(unique(bdf$sample[bdf$sample %in% sels]))
    rows  <- lapply(samps, function(s) {
      sc <- intersect(bdf$cell[bdf$sample==s], cells)
      if (length(sc)==0) return(NULL)
      sub  <- mat[vars, sc, drop=FALSE]
      freq <- rowMeans(sub==1|sub==2, na.rm=TRUE)*100
      data.frame(variant=vars, pct=freq, sample=s, stringsAsFactors=FALSE)
    })
    df <- do.call(rbind, Filter(Negate(is.null), rows))
    df$variant <- factor(df$variant, levels=rev(top_vars()))
    samp_pal <- setNames(
      colorRampPalette(c("#2166AC","#4dac26","#d73027","#f4a40e","#6A4C93"))(length(samps)),
      samps)
    p <- plot_ly(type="bar")
    for (s in samps) {
      sub <- df[df$sample==s,]
      p   <- add_trace(p, x=sub$pct, y=sub$variant, name=s,
                       type="bar", orientation="h",
                       marker=list(color=samp_pal[s]),
                       hovertemplate=paste0("<b>",s,"</b><br>%{y}<br>%{x:.1f}%<extra></extra>"))
    }
    p %>% layout(barmode="stack",
                 paper_bgcolor="#0d1117",plot_bgcolor="#0d1117",
                 font=list(color="#e6edf3",family="IBM Plex Mono"),
                 xaxis=list(title="Mutation %",tickfont=list(size=9,color="#8b949e"),
                            gridcolor="#1e2530",zeroline=FALSE),
                 yaxis=list(tickfont=list(size=8,color="#8b949e"),showgrid=FALSE),
                 legend=list(bgcolor="#161b22",bordercolor="#30363d",borderwidth=1,font=list(size=10)),
                 margin=list(l=140,r=20,t=20,b=50)
    ) %>% config(displaylogo=FALSE)
  })
  
  output$comut_cluster_bar <- renderPlotly({
    bdf  <- base_df()
    sels <- input$comut_samples %||% unique(bdf$sample)
    bdf  <- bdf[bdf$sample %in% sels,]
    mat  <- ngt_mat()
    vars <- top_vars()[seq_len(min(10,length(top_vars())))]
    excl <- isTRUE(input$comut_excl_missing)
    rows <- lapply(vars, function(v) {
      raw <- mat[v, bdf$cell]
      mut <- if (excl) raw==1|raw==2 else raw==1|raw==2|raw==3
      mut[is.na(mut)] <- FALSE
      data.frame(variant=v, cluster=bdf$cluster, mutated=mut, stringsAsFactors=FALSE)
    })
    df  <- do.call(rbind, rows)
    agg <- df %>% group_by(variant,cluster) %>%
      summarise(pct=mean(mutated)*100, .groups="drop")
    cl_ord <- sort(unique(agg$cluster))
    p <- plot_ly(type="bar")
    for (cl in cl_ord) {
      sub <- agg[agg$cluster==cl,]
      p   <- add_trace(p, x=sub$variant, y=sub$pct, name=paste0("Cl ",cl),
                       marker=list(color=CLUSTER_COLORS[cl] %||% "#aaa"),
                       hovertemplate=paste0("<b>Cl ",cl,"</b><br>%{x}<br>%{y:.1f}%<extra></extra>"))
    }
    p %>% layout(barmode="group",
                 paper_bgcolor="#0d1117",plot_bgcolor="#0d1117",
                 font=list(color="#e6edf3",family="IBM Plex Mono"),
                 xaxis=list(tickfont=list(size=9,color="#8b949e"),tickangle=-30,
                            showgrid=FALSE,zeroline=FALSE),
                 yaxis=list(title="% mutated",tickfont=list(size=9,color="#8b949e"),
                            gridcolor="#1e2530",zeroline=FALSE),
                 legend=list(bgcolor="#161b22",bordercolor="#30363d",borderwidth=1,font=list(size=10)),
                 margin=list(l=60,r=20,t=10,b=80)
    ) %>% config(displaylogo=FALSE)
  })
  
  # ═══════════════════════════════════════════════════════════
  #  TAB 3: PROTEIN × NGT
  # ═══════════════════════════════════════════════════════════
  
  pxn_base <- reactive({
    bdf  <- base_df()
    sels <- input$pxn_samples %||% unique(bdf$sample)
    bdf[bdf$sample %in% sels,]
  })
  
  pxn_df <- reactive({
    bdf   <- pxn_base()
    prot  <- req(input$pxn_protein)
    vari  <- req(input$pxn_variant)
    mat_d <- dsb_mat(); mat_n <- ngt_mat()
    req(prot %in% rownames(mat_d))
    raw_ngt <- if (vari %in% rownames(mat_n)) mat_n[vari, bdf$cell]
    else rep(NA_real_, nrow(bdf))
    geno <- ifelse(is.na(raw_ngt), "N/A",
                   ifelse(raw_ngt==0, "WT (0)",
                          ifelse(raw_ngt==1, "Het (1)",
                                 ifelse(raw_ngt==2, "Hom (2)", "Missing (3)"))))
    data.frame(cell=bdf$cell, cluster=bdf$cluster, sample=bdf$sample,
               protein=as.numeric(mat_d[prot, bdf$cell]),
               genotype=geno, stringsAsFactors=FALSE)
  })
  
  pxn_filt <- reactive({
    df <- pxn_df()
    if (isTRUE(input$pxn_excl_miss)) df <- df[df$genotype != "Missing (3)", ]
    if (isTRUE(input$pxn_excl_na))   df <- df[df$genotype != "N/A", ]
    df
  })
  
  output$pxn_violin <- renderPlotly({
    df    <- pxn_filt()
    prot  <- input$pxn_protein
    genos <- intersect(c("WT (0)","Het (1)","Hom (2)"), unique(df$genotype))
    p     <- plot_ly(type="violin")
    for (g in genos) {
      sub <- df[df$genotype==g,]
      col <- NGT_COLORS[g] %||% "#aaa"
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
      legend=list(bgcolor="#161b22",bordercolor="#30363d",borderwidth=1,font=list(size=10)),
      violinmode="group",margin=list(l=60,r=20,t=20,b=50)
    ) %>% config(displaylogo=FALSE)
  })
  
  output$pxn_heatmap <- renderPlotly({
    bdf   <- pxn_base()
    prot  <- req(input$pxn_protein)
    mat_d <- dsb_mat(); mat_n <- ngt_mat()
    req(prot %in% rownames(mat_d))
    freq  <- rowMeans(mat_n[, bdf$cell, drop=FALSE]==1 |
                        mat_n[, bdf$cell, drop=FALSE]==2, na.rm=TRUE)
    top20 <- names(sort(freq, decreasing=TRUE))[1:min(20,length(freq))]
    pv    <- as.numeric(mat_d[prot, bdf$cell])
    rows  <- lapply(top20, function(v) {
      raw  <- if (v %in% rownames(mat_n)) mat_n[v, bdf$cell] else rep(NA, nrow(bdf))
      geno <- ifelse(is.na(raw),"N/A",
                     ifelse(raw==0,"WT (0)",
                            ifelse(raw==1,"Het (1)",
                                   ifelse(raw==2,"Hom (2)","Missing (3)"))))
      data.frame(variant=v,genotype=geno,protein=pv,stringsAsFactors=FALSE)
    })
    df2  <- do.call(rbind, rows)
    df2  <- df2[!df2$genotype %in% c("Missing (3)","N/A"),]
    agg  <- df2 %>% group_by(variant,genotype) %>%
      summarise(med=median(protein,na.rm=TRUE),.groups="drop")
    gl   <- c("WT (0)","Het (1)","Hom (2)")
    zm   <- matrix(NA, length(top20), 3, dimnames=list(top20,gl))
    for (i in seq_len(nrow(agg))) zm[agg$variant[i], agg$genotype[i]] <- agg$med[i]
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
      yaxis=list(tickfont=list(size=8,color="#8b949e"),autorange="reversed",showgrid=FALSE),
      margin=list(l=140,r=20,t=20,b=60)
    ) %>% config(displaylogo=FALSE)
  })
  
  output$pxn_scatter <- renderPlotly({
    bdf   <- pxn_base()
    prot2 <- req(input$pxn_prot2)
    prot3 <- req(input$pxn_prot3)
    vari  <- input$pxn_variant
    by    <- input$pxn_color_by %||% "ngt"
    mat_d <- dsb_mat(); mat_n <- ngt_mat()
    req(prot2 %in% rownames(mat_d), prot3 %in% rownames(mat_d))
    xv <- as.numeric(mat_d[prot2, bdf$cell])
    yv <- as.numeric(mat_d[prot3, bdf$cell])
    if (by=="ngt" && !is.null(vari) && nchar(vari)>0) {
      raw <- if (vari %in% rownames(mat_n)) mat_n[vari, bdf$cell] else rep(NA,nrow(bdf))
      cv  <- ifelse(is.na(raw),"N/A",
                    ifelse(raw==0,"WT (0)",
                           ifelse(raw==1,"Het (1)",
                                  ifelse(raw==2,"Hom (2)","Missing (3)"))))
      pal <- NGT_COLORS; ttl <- vari
    } else if (by=="sample") {
      cv   <- bdf$sample
      samps<- sort(unique(cv))
      pal  <- setNames(
        colorRampPalette(c("#E63946","#2A9D8F","#E9C46A","#457B9D","#6A4C93"))(length(samps)),
        samps)
      ttl  <- "Sample"
    } else {
      cv <- bdf$cluster; pal <- CLUSTER_COLORS; ttl <- "Cluster"
    }
    p <- plot_ly(type="scatter", mode="markers")
    for (cat in unique(cv)) {
      idx <- cv==cat
      p   <- add_trace(p, x=xv[idx], y=yv[idx], name=cat,
                       marker=list(color=pal[cat] %||% "#aaa", size=3, opacity=.6),
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
    bdf   <- pxn_base()
    prot  <- req(input$pxn_protein)
    mat_d <- dsb_mat(); mat_n <- ngt_mat()
    req(prot %in% rownames(mat_d))
    pv   <- as.numeric(mat_d[prot, bdf$cell])
    load <- colSums(mat_n[, bdf$cell, drop=FALSE]==1 |
                      mat_n[, bdf$cell, drop=FALSE]==2, na.rm=TRUE)
    df   <- data.frame(load_bin=factor(pmin(load,8L)), protein=pv, sample=bdf$sample)
    plot_ly(df, x=~load_bin, y=~protein, type="box",
            color=~load_bin, colors=viridis(9),
            hovertemplate="Load %{x}<br>%{y:.2f}<extra></extra>",
            showlegend=FALSE
    ) %>% layout(
      paper_bgcolor="#0d1117",plot_bgcolor="#0d1117",
      font=list(color="#e6edf3",family="IBM Plex Mono"),
      xaxis=list(title="Mutation Load",tickfont=list(size=9,color="#8b949e"),
                 showgrid=FALSE,zeroline=FALSE),
      yaxis=list(title=prot,tickfont=list(size=9,color="#8b949e"),
                 gridcolor="#1e2530",zeroline=FALSE),
      margin=list(l=60,r=20,t=20,b=50)
    ) %>% config(displaylogo=FALSE)
  })
  
  output$pxn_stats <- renderText({
    df   <- pxn_filt()
    prot <- input$pxn_protein %||% ""
    vari <- input$pxn_variant %||% ""
    if (nrow(df)==0) return("No data after filters")
    samp_n <- table(df$sample)
    samp_l <- paste(mapply(function(s,n) sprintf("  %-12s %d",s,n),
                           names(samp_n), as.integer(samp_n)), collapse="\n")
    tbl    <- table(df$genotype)
    geno_l <- paste(mapply(function(nm,ct) sprintf("  %-12s %d",nm,ct),
                           names(tbl), as.integer(tbl)), collapse="\n")
    wm  <- median(df$protein[df$genotype=="WT (0)"],  na.rm=TRUE)
    hm  <- median(df$protein[df$genotype=="Het (1)"], na.rm=TRUE)
    hom <- median(df$protein[df$genotype=="Hom (2)"], na.rm=TRUE)
    paste0("Protein: ",prot,"\nVariant: ",vari,
           "\n\nSamples:\n",samp_l,
           "\n\nGenotype N:\n",geno_l,
           "\n\nMedian ",prot,":\n",
           sprintf("  WT:  %.2f\n  Het: %.2f\n  Hom: %.2f",
                   wm %||% NA, hm %||% NA, hom %||% NA))
  })
}

shinyApp(ui, server)