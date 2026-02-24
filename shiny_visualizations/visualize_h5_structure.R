# HDF5 Explorer — R Shiny
# Install: install.packages(c("shiny", "DT")); BiocManager::install("rhdf5")
# Run: shiny::runApp("app.R")

library(shiny)
library(rhdf5)
library(DT)

FILE <- "/Users/kateplas/Library/Mobile Documents/com~apple~CloudDocs/Documents/Miles_Lab/BRAF/Pt_H.dna+protein.h5"

# ── Helpers ───────────────────────────────────────────────────────────────────

fmt_bytes <- function(b) {
  if (b < 1024)       paste0(b, " B")
  else if (b < 1048576)    paste0(round(b/1024, 1), " KB")
  else if (b < 1073741824) paste0(round(b/1048576, 1), " MB")
  else                     paste0(round(b/1073741824, 2), " GB")
}

build_tree_nodes <- function(filepath) {
  info  <- h5ls(filepath, all = TRUE)
  nodes <- vector("list", nrow(info))
  for (i in seq_len(nrow(info))) {
    row   <- info[i, ]
    path  <- if (row$group == "/") paste0("/", row$name)
              else paste0(row$group, "/", row$name)
    pid   <- if (row$group == "/") "#" else row$group
    is_ds <- row$otype == "H5I_DATASET"
    dim_s <- if (is_ds && nchar(trimws(row$dim)) > 0)
               paste0("  [", gsub(" x ", "\u00d7", trimws(row$dim)), "]") else ""
    nodes[[i]] <- list(
      id     = path,
      parent = pid,
      text   = paste0(row$name, dim_s),
      icon   = if (is_ds) "ds-icon" else "grp-icon",
      data   = list(is_ds = is_ds, path = path,
                    dtype = if (is_ds) row$dclass else "",
                    dim   = if (is_ds) trimws(row$dim) else "")
    )
  }
  nodes
}

get_preview <- function(filepath, path, max_rows = 12, max_cols = 10) {
  tryCatch({
    info <- h5ls(filepath)
    info$full_path <- ifelse(info$group == "/",
                             paste0("/", info$name),
                             paste0(info$group, "/", info$name))
    row <- info[info$full_path == path, ]
    if (nrow(row) == 0) return(list(kind = "error", msg = "Path not found"))

    is_ds <- row$otype[1] == "H5I_DATASET"
    if (!is_ds) {
      children <- info[info$group == path, ]
      return(list(kind = "group", path = path,
                  n_children = nrow(children),
                  children = children$name))
    }

    dclass <- row$dclass[1]
    dim_str <- trimws(row$dim[1])
    # rhdf5 h5ls dim strings vary: "7412 x 312", "( 7412, 312)", or "7412"
    # Parse robustly, then reverse (Fortran->C order) to match HDF5/Python convention
    dims <- if (nchar(dim_str) == 0) integer(0) else {
      # strip parens, split on ' x ', ',', or whitespace runs
      clean <- gsub("[()\\[\\]]", "", dim_str)
      parts <- trimws(unlist(strsplit(clean, "[, ]+|\\s+x\\s+", perl=TRUE)))
      parts <- parts[nchar(parts) > 0]
      vals  <- suppressWarnings(as.integer(parts))
      if (any(is.na(vals))) integer(0) else rev(vals)
    }
    nbytes <- if (length(dims) == 0) 8L else prod(dims) * 8L

    # Scalar
    if (length(dims) == 0 || all(dims == 0)) {
      val <- h5read(filepath, path)
      return(list(kind = "scalar", value = as.character(val),
                  dtype = dclass, bytes = nbytes, shape = dims))
    }

    # String
    if (dclass == "STRING") {
      n   <- min(60, dims[1])
      idx <- if (length(dims) == 1) list(seq_len(n)) else list(seq_len(n), seq_len(dims[2]))
      val <- tryCatch(h5read(filepath, path, index = list(seq_len(n))),
                      error = function(e) h5read(filepath, path)[seq_len(n)])
      return(list(kind = "strings", values = as.character(as.vector(val)),
                  dtype = dclass, bytes = nbytes, shape = dims))
    }

    # 1-D numeric
    if (length(dims) == 1) {
      n   <- min(200, dims[1])
      idx <- sort(sample(dims[1], n, replace = FALSE))
      val <- as.numeric(h5read(filepath, path, index = list(idx)))
      stats <- list(min = min(val, na.rm=T), max = max(val, na.rm=T),
                    mean = mean(val, na.rm=T), median = median(val, na.rm=T))
      return(list(kind = "1d", values = val, stats = stats,
                  dtype = dclass, bytes = nbytes, shape = dims))
    }

    # 2-D numeric
    # dims are now [nrows, ncols] in HDF5/Python convention.
    # h5read index list is [col_idx, row_idx] (Fortran order), so we swap.
    if (length(dims) == 2) {
      r  <- min(max_rows, dims[1])   # rows in HDF5 sense
      cc <- min(max_cols, dims[2])   # cols in HDF5 sense
      ri <- sort(sample(dims[1], r, replace = FALSE))
      # h5read index: first index = cols (dims[2]), second = rows (dims[1])
      mat  <- t(h5read(filepath, path, index = list(seq_len(cc), ri)))
      # wider sample for stats (read cols first, then rows)
      r2   <- min(300, dims[1])
      ri2  <- sort(sample(dims[1], r2, replace = FALSE))
      samp <- h5read(filepath, path, index = list(seq_len(min(30, dims[2])), ri2))
      sv   <- as.numeric(samp)
      stats <- list(min = min(sv, na.rm=T), max = max(sv, na.rm=T),
                    mean = mean(sv, na.rm=T), median = median(sv, na.rm=T))
      df <- as.data.frame(mat)
      colnames(df) <- paste0("col_", seq_len(ncol(df)))
      return(list(kind = "2d", data = df, stats = stats,
                  rows_total = dims[1], cols_total = dims[2],
                  rows_shown = r, cols_shown = cc,
                  dtype = dclass, bytes = nbytes, shape = dims))
    }

    # ND: 2-D slice of last two dims
    mid  <- ceiling(dims[seq_len(length(dims)-2)] / 2)
    idx_l <- c(as.list(mid),
               list(seq_len(min(max_rows, dims[length(dims)-1])),
                    seq_len(min(max_cols, dims[length(dims)]))))
    sl <- h5read(filepath, path, index = idx_l)
    df <- as.data.frame(sl)
    colnames(df) <- paste0("col_", seq_len(ncol(df)))
    return(list(kind = "nd_slice", data = df, slice_at = mid,
                dtype = dclass, bytes = nbytes, shape = dims))

  }, error = function(e) list(kind = "error", msg = conditionMessage(e)))
}

# ── CSS ───────────────────────────────────────────────────────────────────────

CSS <- "
@import url('https://fonts.googleapis.com/css2?family=JetBrains+Mono:wght@300;400;600&family=Syne:wght@700;800&display=swap');
* { box-sizing: border-box; margin: 0; padding: 0; }
html, body, .shiny-page, #app-root { height: 100%; background: #080c10 !important; }

.app-shell  { display:flex; flex-direction:column; height:100vh; overflow:hidden;
              font-family:'JetBrains Mono',monospace; color:#cdd9e5; background:#080c10; }

/* header */
.hdr { background:#0d1117; border-bottom:1px solid #1e2633; padding:11px 22px;
       display:flex; align-items:center; gap:14px; flex-shrink:0; }
.hdr-title { font-family:'Syne',sans-serif; font-weight:800; font-size:16px; color:#58a6ff; }
.hdr-file  { color:#545d68; font-size:11px; margin-top:2px; }
.hdr-badge { background:#182030; border:1px solid #1e2633; border-radius:4px;
             padding:2px 9px; font-size:10px; color:#545d68; }
.search-wrap { margin-left:auto; position:relative; }
#node_search { background:#080c10; border:1px solid #1e2633; border-radius:6px;
               color:#cdd9e5; font-family:'JetBrains Mono',monospace; font-size:12px;
               padding:5px 10px 5px 26px; width:210px; outline:none; }
#node_search:focus { border-color:#58a6ff; }
.search-wrap::before { content:'⌕'; position:absolute; left:8px; top:50%;
                        transform:translateY(-50%); color:#545d68; pointer-events:none; }

/* body split */
.body-split { display:flex; flex:1; overflow:hidden; }

/* tree */
.tree-pane { width:360px; min-width:200px; background:#0d1117;
             border-right:1px solid #1e2633; overflow-y:auto; padding:6px 0; flex-shrink:0; }

/* jstree overrides */
.jstree-default .jstree-anchor { color:#cdd9e5 !important; font-size:12px;
  font-family:'JetBrains Mono',monospace !important; }
.jstree-default .jstree-hovered { background:#131c26 !important; border-radius:4px !important; box-shadow:none !important; }
.jstree-default .jstree-clicked { background:#132040 !important; border-radius:4px !important; box-shadow:none !important; }
.jstree-default { background:#0d1117 !important; }
.grp-icon:before { content:'⬡'; font-style:normal; color:#4fa3e8; }
.ds-icon:before  { content:'▪'; font-style:normal; color:#3fb950; }
.jstree-icon { font-style:normal !important; }

/* preview */
.preview-pane { flex:1; overflow-y:auto; padding:26px 30px; background:#080c10; }
.ph { display:flex; flex-direction:column; align-items:center; justify-content:center;
      height:70vh; color:#545d68; gap:10px; font-size:13px; }
.ph .big { font-size:44px; }

.prev-title { font-family:'Syne',sans-serif; font-weight:800; font-size:15px; color:#58a6ff; margin-bottom:3px; }
.prev-path  { color:#545d68; font-size:11px; padding-bottom:14px; border-bottom:1px solid #1e2633; margin-bottom:18px; }

.meta-grid { display:grid; grid-template-columns:repeat(auto-fill,minmax(140px,1fr));
             gap:10px; margin-bottom:20px; }
.meta-card { background:#0d1117; border:1px solid #1e2633; border-radius:8px; padding:10px 14px; }
.meta-card .k { font-size:10px; color:#545d68; margin-bottom:3px; text-transform:uppercase; letter-spacing:.5px; }
.meta-card .v { font-size:13px; color:#cdd9e5; font-weight:600; }

.sec { font-family:'Syne',sans-serif; font-size:10px; font-weight:700; color:#545d68;
       text-transform:uppercase; letter-spacing:1px; margin-bottom:9px; margin-top:4px; }

.stats-row { display:flex; gap:9px; margin-bottom:18px; flex-wrap:wrap; }
.spill { background:#182030; border:1px solid #1e2633; border-radius:20px;
         padding:4px 13px; font-size:11px; }
.spill .l { color:#c9921a; margin-right:3px; }

.str-grid { display:flex; flex-wrap:wrap; gap:6px; margin-bottom:16px; }
.schip { background:#0d1117; border:1px solid #1e2633; border-radius:4px;
         padding:3px 10px; font-size:11px; color:#3fb950; }
.gchip { color:#4fa3e8; }

.spark { display:flex; align-items:flex-end; gap:2px; height:58px; margin-bottom:16px; flex-wrap:wrap; }
.sbar  { flex:1; min-width:3px; max-width:15px; background:#58a6ff; opacity:.65;
         border-radius:2px 2px 0 0; transition:opacity .1s; cursor:default; }
.sbar:hover { opacity:1; }

.trunc { font-size:10px; color:#545d68; margin-top:5px; }

/* DT */
table.dataTable { background:#0d1117 !important; color:#cdd9e5 !important;
  font-size:12px; font-family:'JetBrains Mono',monospace !important; border:none !important; }
table.dataTable thead th { background:#0d1520 !important; color:#c9921a !important;
  border-bottom:1px solid #1e2633 !important; font-size:10px; text-transform:uppercase; }
table.dataTable tbody tr { background:#0d1117 !important; }
table.dataTable tbody tr:hover td { background:#131c26 !important; }
table.dataTable tbody td { border-bottom:1px solid #111820 !important; padding:5px 10px !important; }
.dataTables_info, .dataTables_filter label, .dataTables_length label { color:#545d68 !important; font-size:11px; }
.dataTables_filter input { background:#0d1117 !important; border:1px solid #1e2633 !important;
  color:#cdd9e5 !important; border-radius:4px; padding:3px 7px; font-size:11px; }
.paginate_button { color:#545d68 !important; background:transparent !important; border:none !important; font-size:11px; }
.paginate_button.current { color:#58a6ff !important; font-weight:bold; }
.dataTables_wrapper { background:transparent !important; }

::-webkit-scrollbar { width:5px; height:5px; }
::-webkit-scrollbar-track { background:transparent; }
::-webkit-scrollbar-thumb { background:#2a3240; border-radius:3px; }
"

# ── UI ────────────────────────────────────────────────────────────────────────

ui <- fluidPage(
  tags$head(
    tags$style(HTML(CSS)),
    tags$link(rel="stylesheet",
      href="https://cdnjs.cloudflare.com/ajax/libs/jstree/3.3.16/themes/default/style.min.css"),
    tags$script(src="https://cdnjs.cloudflare.com/ajax/libs/jquery/3.7.1/jquery.min.js"),
    tags$script(src="https://cdnjs.cloudflare.com/ajax/libs/jstree/3.3.16/jstree.min.js")
  ),

  div(class="app-shell",
    # Header
    div(class="hdr",
      div(
        div(class="hdr-title", "HDF5 Explorer"),
        div(class="hdr-file",  basename(FILE))
      ),
      div(class="hdr-badge", id="node_count_lbl", "loading…"),
      div(class="search-wrap",
        tags$input(id="node_search", type="text", placeholder="Search nodes…",
                   oninput="searchTree(this.value)")
      )
    ),

    # Body
    div(class="body-split",
      div(class="tree-pane", div(id="jstree_div")),
      div(class="preview-pane", uiOutput("preview_ui"))
    )
  ),

  tags$script(HTML("
    var searchTimer;
    function searchTree(q) {
      clearTimeout(searchTimer);
      searchTimer = setTimeout(function() {
        $('#jstree_div').jstree(true).search(q);
      }, 200);
    }

    Shiny.addCustomMessageHandler('build_tree', function(nodes) {
      document.getElementById('node_count_lbl').textContent = nodes.length + ' nodes';
      $('#jstree_div').jstree('destroy').empty();
      $('#jstree_div').jstree({
        core: { data: nodes, themes: { name:'default', dots:true, icons:true } },
        plugins: ['search', 'wholerow'],
        search: { show_only_matches: true, show_only_matches_children: true }
      });
      $('#jstree_div').on('select_node.jstree', function(e, d) {
        Shiny.setInputValue('selected_path', d.node.data.path, {priority:'event'});
      });
    });
  "))
)

# ── Server ────────────────────────────────────────────────────────────────────

server <- function(input, output, session) {

  observe({
    nodes <- build_tree_nodes(FILE)
    jsnodes <- lapply(nodes, function(n) {
      list(id=n$id, parent=n$parent, text=n$text, icon=n$icon, data=n$data,
           state=list(opened=(n$parent=="#")))
    })
    session$sendCustomMessage("build_tree", jsnodes)
  })

  output$preview_ui <- renderUI({
    path <- req(input$selected_path)

    p <- get_preview(FILE, path)

    if (p$kind == "error")
      return(div(class="prev-title", style="color:#f78166;", "Error: ", p$msg))

    name <- basename(path)
    hdr  <- tagList(div(class="prev-title", name),
                    div(class="prev-path",  path))

    # Group
    if (p$kind == "group") {
      chips <- lapply(p$children, function(c)
        span(class="schip gchip", paste0("⬡ ", c)))
      return(tagList(hdr,
        div(class="meta-grid",
          div(class="meta-card", div(class="k","Type"),     div(class="v","Group")),
          div(class="meta-card", div(class="k","Children"), div(class="v", p$n_children))
        ),
        div(class="sec","Children"),
        do.call(div, c(list(class="str-grid"), chips))
      ))
    }

    shape_str <- if (length(p$shape)==0) "scalar" else paste(p$shape, collapse=" \u00d7 ")
    meta <- div(class="meta-grid",
      div(class="meta-card", div(class="k","Shape"),  div(class="v", shape_str)),
      div(class="meta-card", div(class="k","dtype"),  div(class="v", p$dtype)),
      div(class="meta-card", div(class="k","Size"),   div(class="v", fmt_bytes(p$bytes))),
      div(class="meta-card", div(class="k","Kind"),   div(class="v", p$kind))
    )

    # Scalar
    if (p$kind == "scalar")
      return(tagList(hdr, meta, div(class="sec","Value"),
        span(class="schip", style="color:#cdd9e5;font-size:14px;", p$value)))

    # Strings
    if (p$kind == "strings") {
      chips <- lapply(head(p$values, 60), function(v) span(class="schip", v))
      return(tagList(hdr, meta,
        div(class="sec", paste0("Values (first ", length(chips), ")")),
        do.call(div, c(list(class="str-grid"), chips))))
    }

    # 1-D
    if (p$kind == "1d") {
      vals <- p$values
      mn <- min(vals,na.rm=T); mx <- max(vals,na.rm=T); rng <- max(mx-mn, 1e-9)
      bars <- lapply(vals, function(v) {
        h <- max(4, round(((v-mn)/rng)*52))
        div(class="sbar", style=paste0("height:",h,"px"), title=signif(v,5))
      })
      sr <- div(class="stats-row",
        div(class="spill", span(class="l","min"),    signif(p$stats$min,5)),
        div(class="spill", span(class="l","max"),    signif(p$stats$max,5)),
        div(class="spill", span(class="l","mean"),   signif(p$stats$mean,5)),
        div(class="spill", span(class="l","median"), signif(p$stats$median,5))
      )
      chips <- lapply(head(vals,80), function(v)
        span(class="schip", style="color:#cdd9e5;", signif(v,5)))
      return(tagList(hdr, meta, sr,
        div(class="sec", paste0("Distribution (", length(vals), " samples)")),
        do.call(div, c(list(class="spark"), bars)),
        div(class="sec","Raw values"),
        do.call(div, c(list(class="str-grid"), chips))))
    }

    # 2-D / ND slice
    if (p$kind %in% c("2d","nd_slice")) {
      sr <- NULL
      if (!is.null(p$stats)) {
        sr <- div(class="stats-row",
          div(class="spill", span(class="l","min"),    signif(p$stats$min,5)),
          div(class="spill", span(class="l","max"),    signif(p$stats$max,5)),
          div(class="spill", span(class="l","mean"),   signif(p$stats$mean,5)),
          div(class="spill", span(class="l","median"), signif(p$stats$median,5))
        )
      }
      sec_lbl <- if (p$kind=="nd_slice")
        div(class="sec", paste0("2D slice at [", paste(p$slice_at,collapse=","), "]"))
      else
        div(class="sec", paste0(p$rows_shown," rows \u00d7 ",p$cols_shown," cols shown"))

      note <- if (!is.null(p$rows_shown))
        div(class="trunc", paste0("Showing ",p$rows_shown," of ",p$rows_total,
            " rows, ",p$cols_shown," of ",p$cols_total," cols — random sample"))
      else NULL

      return(tagList(hdr, meta, sr, sec_lbl, DTOutput("tbl"), note))
    }
  })

  output$tbl <- renderDT({
    path <- req(input$selected_path)
    p    <- get_preview(FILE, path)
    if (!p$kind %in% c("2d","nd_slice") || is.null(p$data)) return(NULL)
    datatable(p$data,
      options = list(pageLength=10, scrollX=TRUE, dom="ftp",
        initComplete = JS("function(s,j){
          $(this.api().table().header()).css({'background':'#0d1520','color':'#c9921a'});
        }")),
      rownames = TRUE, class = "compact hover"
    )
  })
}

shinyApp(ui, server)