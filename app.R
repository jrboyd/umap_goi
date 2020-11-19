#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(data.table)
library(ggplot2)
library(magrittr)
library(cowplot)
library(shinycssloaders)
library(tippy)

source("tippy.R")
source("functions_setup.R")

gene_name_gs2mm = function(x){
    x = tolower(x)
    substr(x, 1, 1) = toupper(substr(x, 1, 1))
    x = sub("(?<=[0-9])rik$", "Rik", x, perl = TRUE)
    is_rik = grepl("Rik$", x)
    substr(x[is_rik], 8, 8) = toupper(substr(x[is_rik], 8, 8))
    x
}

js <- '$("#selDataSource :input").each(function() {
    $(this).attr("id", "radio_" + $(this).val());
});'

js2 <- '$("#selDataSource span").each(function() {
    $(this).attr("id", "radioLabel_" + $(this).text());
});'
js2a <- '$("#selDataSource span").each(function() {
    $(this).parents()[1].id = "radioItem_" + $(this).text();
});'

js3 <- '$("#txtGene").parents()[1].id = "txtGene-container"'
js_selGeneLookup <- '$("#selGeneLookup").parents()[1].id = "selGeneLookup-container"'
# get_meta_dt = function(seurat_obj){
#     meta_dt = as.data.table(seurat_obj@meta.data)
#     meta_dt$id = seurat_obj@meta.data %>% rownames
#     
#     umap_dt = as.data.table(seurat_obj@reductions$umap@cell.embeddings, keep.rownames = TRUE) %>% setnames(., "rn", "id")
#     meta_dt =merge(meta_dt, umap_dt, by = "id")
#     meta_dt$orig.ident = factor(meta_dt$orig.ident, levels = c("wt", 'df4'))
#     meta_dt
# }
# 
# get_rna_dt = function(seurat_obj, sel_genes = NULL){
#     if(is.null(sel_genes)){
#         rna_dt = as.data.frame(seurat_obj@assays$RNA)
#     }else{
#         rna_dt = as.data.frame(seurat_obj@assays$RNA[sel_genes, ])
#         
#     }
#     rna_dt$gene_name = rownames(rna_dt)
#     rna_dt = melt(as.data.table(rna_dt), variable.name = "id", value.name = "expression", id.vars = "gene_name")
#     rna_dt
# }

#bulk rnaseq
cnt_dt = ssvRecipes::bfcif(BiocFileCache::BiocFileCache(), "umap_goi_bulk_v1", function(){
    cnt_res = setup_RNA_counts_with_proB()
    cnt_dt = cnt_res$tidy
    cnt_dt[, sum(count), .(file)]
    cnt_dt[, cpm := count / sum(count) * 1e6, .(file)]
    cnt_dt[, zscore := (cpm - mean(cpm)) / sd(cpm), .(gene_name)]
    ko = cnt_dt[is.na(zscore)]$gene_name %>% unique
    cnt_dt = cnt_dt[!gene_name %in% ko]
    cnt_dt = cnt_dt[order(rep)][order(treatment)][order(status)][order(stage)]
    cnt_dt$file = factor(cnt_dt$file, levels = unique(as.character(cnt_dt$file)))
    
    cnt_dt[, xgroup := paste(stage, treatment)]
    cnt_dt$xgroup = factor(cnt_dt$xgroup, levels = unique(as.character(cnt_dt$xgroup)))
    cnt_dt    
})

pbmc_sources_dir = "/home/user/Documents/R_workspace/Bcell_IKdf4"
if(!dir.exists(pbmc_sources_dir)){
    pbmc_sources_dir = "/slipstream/home/joeboyd/R/SF_AutoImmune_ssv/"  
}
if(!dir.exists(pbmc_sources_dir)) stop("bad dir")
# 
remove_toy = TRUE
pbmc_sources = list(
    "toy" = "toy.rds",
    "refix" = "Bcell_combined.refixed_120919.rds",
    "SF_original" = "Bcell_combined.rds", 
    "DEGs" = "Bcell_combined.DEGs_umap.rds", 
    
    "refix_Bcell" = "Bcell_combined.refixed_refined_012720.rds")


pbmc_loaded = lapply(pbmc_sources, function(x)NULL)
pbmc_sources = lapply(pbmc_sources, function(x)file.path(pbmc_sources_dir, x))

if(!file.exists(pbmc_sources$toy)){
    pbmc_small = readRDS(pbmc_sources$refix)
    DefaultAssay(pbmc_small) = "RNA"
    k = grepl("^Cd", rownames(pbmc_small))
    pbmc_small = pbmc_small[k,]
    k = sample(seq(ncol(pbmc_small), 500))
    pbmc_small = pbmc_small[k,]
    saveRDS(pbmc_small, pbmc_sources$toy)
}



k = sapply(pbmc_sources, file.exists)
if(!all(k)) warning("some data sources not found, will be ignored.")
pbmc_loaded = pbmc_loaded[k]
pbmc_sources = pbmc_sources[k]

if(remove_toy){
    pbmc_loaded$toy = NULL
    pbmc_sources$toy = NULL
}else{
    pbmc_loaded$toy = readRDS(pbmc_sources$toy)    
}


stopifnot(names(pbmc_loaded) == names(pbmc_sources))

coarse_file = file.path(pbmc_sources_dir, "coarse_scores.csv")
coarse_dt = fread(coarse_file)
fine_file = file.path(pbmc_sources_dir, "fine_scores.csv")
fine_dt = fread(fine_file)

module_genes_dt = openxlsx::read.xlsx(file.path(pbmc_sources_dir, "Haemopedia/Immgen_gene_assignment.xlsx")) %>% as.data.table
module_genes_dt[, gene_name := gene_name_gs2mm(Gene)]
coarse_genes_dt = module_genes_dt[, .(gene_name, module = paste0("coarse_", Coarse.module))]
fine_genes_dt = module_genes_dt[, .(gene_name, module = paste0("fine_", Fine.module))]

coarse_descriptions_dt = openxlsx::read.xlsx(file.path(pbmc_sources_dir, "Haemopedia/Immgen_coarse_module_descriptions.NIHMS455856-supplement-4.xlsx"), startRow = 2) %>% as.data.table
coarse_descriptions_dt[,module := paste0("coarse_", as.character(module.number))]

fine_descriptions_dt = openxlsx::read.xlsx(file.path(pbmc_sources_dir, "Haemopedia/Immgen_fine_module_descriptions.NIHMS455856-supplement-5.xlsx"), startRow = 2) %>% as.data.table
fine_descriptions_dt[,module := paste0("fine_", as.character(Fine.module))]

# Define UI for application that draws a histogram
ui <- fluidPage(
    shinyjs::useShinyjs(),
    tags$a(href="https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3690947/", "Immgen Modules Paper"),
    actionButton("btnSelectUMAP", "Select UMAP"),
    checkboxInput("checkSplitGenotype", "Split By Genotype", value = TRUE),
    tippy_this("checkSplitGenotype", "If checked, all plots will facet data based on genotype."),
        # Application title
    titlePanel("scRNAseq of wt and IKdF4 B-cells"),
    
    tabsetPanel(id = "tabset1", 
                # Sidebar with a slider input for number of bins
                tabPanel("Gene UMAP Query", id = "tag_goi",
                         sidebarLayout(
                             sidebarPanel(
                                 selectizeInput(inputId = "txtGene", label = "Select Genes", choices = NULL, multiple = FALSE),
                                 tags$script(js3),
                                 tippy_this(
                                     "txtGene-container",
                                     tooltip = "Delete and start typing gene names to see matching genes."
                                 )
                             ),
                             
                             # Show a plot of the generated distribution
                             mainPanel(
                                 withSpinner(plotOutput("plotGOI", width = "550px", height = "310px")),
                                 tippy_this(
                                     "plotGOI",
                                     tooltip = paste("scRNAseq UMAP facetted by mutation status of cells with",
                                                     "expression values of the selected gene mapped to color.", 
                                                     "Select a different UMAP or gene on the left.",
                                                     "Right-click and save or copy image to export plot.")
                                 ),
                                 withSpinner(
                                     plotOutput("plotGOI_perCluster", width = "550px", height = "380px")
                                 ),
                                 withSpinner(
                                     plotOutput("plotBulkGOI", width = "550px", height = "310px")
                                 )
                             )
                         )
                ),
                tabPanel("Explore Immgen Modules", id = "tab_immgen",
                         sidebarLayout(
                             sidebarPanel(
                                 radioButtons("selModuleResolution", label = "Module Resolution", choices = c("coarse", "fine")),
                                 tippy_this("selModuleResolution", "Immgen Modules exist at coarse (81 modules) and fine (334 modules) resolutions.  Fine modules are derived from coarse."),
                                 radioButtons("selClusterSet", 
                                              label = "scRNAseq Clustering Method", 
                                              choiceNames = c("seurat_clusters", 
                                                              "meta_cluster", 
                                                              "seurat_clusters+genotype",
                                                              "meta_cluster+genotype"), 
                                              choiceValues = c("seurat_clusters", 
                                                               "meta_cluster", 
                                                               "seurat_clustersAndGenotype",
                                                               "meta_clusterAndGenotype")),
                                 tippy_this("selClusterSet", paste("scRNAseq cell clustering are provided at 2 resolutions.", 
                                                                   "There are 17 Seurat clusters and 5 meta clusters created", 
                                                                   "by comibining the most similar Seurat clusters.",  
                                                                   "Variants of these two clusterings are provided that are further split by cell genotype.")),
                                 uiOutput("clustersAvailable"),
                                 tippy_this("clustersAvailable", paste("Based on your cell clustering method, select which clusters to calculate ranking stats on.")),
                                 withSpinner(plotOutput("plotMiniClusters", width = "210px", height = "200px")),
                                 tippy_this("plotMiniClusters", paste("UMAP projection highlighting selected cell clusters.")),
                                 radioButtons("selClusterMethod",
                                              label = "Ranking Method",
                                              choices = c("identity", "difference", "variance", "average", "max", "min")),
                                 tippy_this("selClusterMethod", paste("Ranking metrics to provide in module browsing table.")),
                                 selectInput("selGeneLookup", label = "Lookup By Gene", 
                                             choices = sort(module_genes_dt$gene_name), selected = "Cd19"),
                                 tags$script(js_selGeneLookup),
                                 tippy_this("selGeneLookup-container", paste("Type a gene name here to see what Immgen module it belongs to. Result is specific to current Module Resolution.")),
                                 withSpinner(uiOutput("moduleLookup")),
                                 textInput("txtModuleQuery", "Module Query", value = ""),
                                 tippy_this("txtModuleQuery", paste("Type an Immgen Module id here to manually specify. Result is specific to current Module Resolution.")),
                                 tags$span("invalid module ID", id = "msgBadModule", style = "color:red"),
                                 tags$span("valid module ID", id = "msgGoodModule", style = "color:green"),
                                 withSpinner(plotOutput("plotMiniModule", width = "210px", height = "200px")),
                                 tippy_this("plotMiniModule", paste("UMAP projection with selected odule scores mapped to color.")),
                                 actionButton("btnViewModuleDetail", label = "Module Details"),
                                 tippy_this("btnViewModuleDetail", paste("Launch a more detailed gene by gene browser for modules."))
                             ),
                             mainPanel(
                                 withSpinner(DT::dataTableOutput("coarseTable")),
                                 tippy_this("coarseTable", paste("Sort by ranking metrics or search for terms to explore Immgen Modules. Click on entries to select."))
                             )
                         )
                )
    ),
    tippy_this("tabset1", "Gene UMAP Query\nExplore Immgen Modules")
)

# Define server logic required to draw a histogram
server <- function(input, output, session) {
    shinyjs::hide(id = "msgGoodModule")
    shinyjs::disable(id = "btnViewModuleDetail")
    pbmc_data = reactiveVal(NULL)
    active_source = reactiveVal(names(pbmc_loaded)[1])
    
    observeEvent(input$btnSelectUMAP, {
        showModal(modalDialog(
            radioButtons(inputId = "selDataSource", label = "Select UMAP", choices = names(pbmc_loaded), selected = active_source()),
            tags$script(js),
            tags$script(js2),
            tags$script(js2a),
            tippy_datasets(),
            actionButton("btnLoadModal", "Load"),
            title = "Select a different UMAP",
            "These can take a minute or two to load.",
            easyClose = TRUE,
            footer = NULL
        ))
    })
    
    observeEvent({
        input$btnLoadModal
    }, {
        removeModal()
        sel = input$selDataSource
        stopifnot(sel %in% names(pbmc_loaded))
        active_source(sel)
        
    })
    
    observeEvent({
        active_source()
    }, {
        sel = active_source()
        if(is.null(pbmc_loaded[[sel]])){
            showNotification(paste0("loading ", sel, " UMAP..."), duration = NULL, id = "load_note")
            pbmc_loaded[[sel]] = readRDS(pbmc_sources[[sel]])
            Seurat::DefaultAssay(pbmc_loaded[[sel]]) = "RNA"
            removeNotification(id = "load_note")
        }
        
        pbmc_data(pbmc_loaded[[sel]])
        #invalidate derived data
        meta_dtR(NULL)
        module_dtR(NULL)
        
        all_genes = sort(rownames(pbmc_data()@assays$RNA))
        start = min(which(grepl("^Cd", all_genes)))
        
        all_genes = all_genes[c(seq(start, length(all_genes)), seq(1, start - 1))]
        def = "Cd84"
        if(def %in% all_genes){
            updateSelectizeInput(session, 'txtGene', choices = all_genes, selected = def, server = TRUE)    
        }else{
            updateSelectizeInput(session, 'txtGene', choices = all_genes, server = TRUE)    
        }
        
    })
    
    output$plotBulkGOI = renderPlot({
        goi = input$txtGene
        req(goi)
        ggplot(cnt_dt[gene_name == goi], aes(x = xgroup, y = cpm, color = status)) +
            geom_jitter(width = .05) +
            scale_color_manual(values = c("WT" = "black", "dF4" = "red")) +
            theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
            labs(title = paste(goi, "in bulk RNA-seq"), x = "group")
        
    })
    
    output$plotGOI <- renderPlot({
        goi = input$txtGene
        req(goi)
        
        x = pbmc_data()
        req(x)
        req(goi != '')
        meta_dt = get_meta_dt(x)
        meta_dt$orig.ident = factor(meta_dt$orig.ident, levels = c("wt", "df4"))
        rna_dt = get_rna_dt(x, goi)
        rna_dt = merge(meta_dt, rna_dt, by = "id")
        rna_dt = rna_dt[order(expression, decreasing = FALSE)]
        p = ggplot(rna_dt, aes(x = UMAP_1, y = UMAP_2, color = expression)) + 
            geom_point(size= .1) +
            theme_cowplot() +
            # facet_wrap(~orig.ident) +
            theme(panel.background = element_rect(fill = "gray70"), 
                  strip.background = element_blank()) +
            scale_color_viridis_c() +
            labs(title = goi, color = "expression")    
        if(input$checkSplitGenotype){
            p = p + facet_wrap(~orig.ident)    
        }
        p
    })
    
    output$plotGOI_perCluster = renderPlot({
        goi = input$txtGene
        req(goi)
        
        x = pbmc_data()
        req(x)
        req(goi != '')
        meta_dt = get_meta_dt(x)
        meta_dt$orig.ident = factor(meta_dt$orig.ident, levels = c("wt", "df4"))
        
        lab_dt = meta_dt[, .(UMAP_1 = median(UMAP_1), UMAP_2 = median(UMAP_2)), .(meta_cluster)]
        
        rna_dt = get_rna_dt(x, goi)
        rna_dt = merge(meta_dt, rna_dt, by = "id")
        rna_dt = rna_dt[order(expression, decreasing = FALSE)]
        p1 = ggplot(rna_dt, aes(x = UMAP_1, y = UMAP_2, color = meta_cluster)) +
            geom_point() +
            geom_label(data = lab_dt, aes(label = meta_cluster), show.legend = FALSE) +
            scale_color_manual(values = get_clusters_colors())
        p2 = ggplot(rna_dt, aes(x = meta_cluster, y = expression, 
                           #fill = meta_cluster, 
                           color = orig.ident)) +
            geom_boxplot() +
            # scale_fill_manual(values = get_clusters_colors()) +
            scale_color_manual(values = c("df4" = 'red', "wt" = "black")) + 
            labs(title = goi) +
            theme(axis.text.x = element_text(color = get_clusters_colors(), size = 14))
        cowplot::plot_grid(p1 + theme(legend.position = "bottom"), p2 + theme(legend.position = "bottom"))
    })
    
    loaded_coarse = reactiveVal()
    loaded_fine = reactiveVal()
    meta_dtR = reactiveVal()
    
    observe({
        x = req(pbmc_data())
        module_resolution = req(input$selModuleResolution)
        
        meta_dt = get_meta_dt(x)
        
        td = list(
            "main_continent" = c(7, 0, 2, 16, 2, 1, 3, 5, 10),
            "north_archipelago" = c(9, 8, 14, 12), 
            "west_island" = c(15, 6), 
            "west_projection" = c(13, 4), 
            "lonely_island" = c(11)
        )
        
        td_dt = lapply(td, function(x){
            data.table(seurat_clusters = factor(x, levels = levels(meta_dt$seurat_clusters)))
        }) %>% rbindlist(idcol = "meta_cluster")
        td_dt$meta_cluster = factor(td_dt$meta_cluster, levels = names(td))
        meta_dt = merge(meta_dt, td_dt, by = 'seurat_clusters', allow.cartesian = TRUE) %>% unique
        meta_dt$orig.ident = factor(meta_dt$orig.ident, levels = c("wt", "df4"))
        
        if(module_resolution == "coarse"){
            cn = colnames(coarse_dt)[grepl("coarse_", colnames(coarse_dt))]
            meta_dt = merge(meta_dt, coarse_dt[, c("id", cn), with = FALSE], by = "id")    
        }else if(module_resolution == "fine"){
            cn = colnames(fine_dt)[grepl("fine_", colnames(fine_dt))]
            meta_dt = merge(meta_dt, fine_dt[, c("id", cn), with = FALSE], by = "id")    
        }
        meta_dt[["seurat_clustersAndGenotype"]] = factor(paste(meta_dt[["seurat_clusters"]], meta_dt[["orig.ident"]]))
        meta_dt[["meta_clusterAndGenotype"]] = factor(paste(meta_dt[["meta_cluster"]], meta_dt[["orig.ident"]]))
        meta_dtR(meta_dt)
    })
    
    
    
    output$clustersAvailable = renderUI({
        cluster_variable = req(input$selClusterSet)
        meta_dt = req(meta_dtR())
        selectInput("selClusterIDs", label = "Select Clusters", choices = levels(meta_dt[[cluster_variable]]), 
                    selected = levels(meta_dt[[cluster_variable]])[1:3],
                    multiple = TRUE)
    })
    
    module_dtR = reactiveVal()
    
    output$coarseTable = DT::renderDT({
        # radioButtons("selClusterMethod",
        #              label = "Method",
        #              choices = c("identity", "variance", "average", "max", "min", "difference"))
        module_resolution = req(input$selModuleResolution)
        meta_dt = req(meta_dtR())
        sel_clust = req(input$selClusterIDs)
        cluster_variable = req(input$selClusterSet)
        sort_method = req(input$selClusterMethod)
        # browser()
        k = as.character(meta_dt[[cluster_variable]]) %in% sel_clust
        if(!any(k)) req(NULL)
        
        if(module_resolution == "coarse"){
            cn = colnames(coarse_dt)[grepl("coarse_", colnames(coarse_dt))]
        }else if(module_resolution == "fine"){
            cn = colnames(fine_dt)[grepl("fine_", colnames(fine_dt))]
        }
        
        suppressWarnings({
            melted = melt(meta_dt[k, ][, c(cluster_variable, cn), with = FALSE], id.vars = cluster_variable, variable.name = "module", measure.vars = cn)
            module_dt = melted[, .(mean_score = mean(value)), c(cluster_variable, "module")]    
        })
        
        if(sort_method == "identity"){
            lev = module_dt[[cluster_variable]] %>% unique %>% sort %>% as.character()
            module_dt = dcast(module_dt, paste0("module~", cluster_variable), value.var = "mean_score")
        }else if(sort_method == "variance"){
            lev = "variance"
            module_dt = module_dt[, .(variance = var(mean_score)), .(module)]
            
        }else if(sort_method == "average"){
            lev = "average"
            module_dt = module_dt[, .(average = mean(mean_score)), .(module)]
            
        }else if(sort_method == "max"){
            lev = "maximum"
            module_dt = module_dt[, .(maximum = max(mean_score)), .(module)]
            
        }else if(sort_method == "min"){
            lev = "minimum"
            module_dt = module_dt[, .(minimum = min(mean_score)), .(module)]
            
        }else if(sort_method == "difference"){
            lev = (module_dt[[cluster_variable]] %>% unique %>% as.character())
            if(length(lev) < 2) showNotification("select 2 clusters for difference", type = "warning")
            if(length(lev) > 2) showNotification("only first 2 clusters used for difference", type = "warning")
            lev = lev[1:2]
            module_dt = dcast(module_dt[get(cluster_variable) %in% lev], paste0("module~", cluster_variable), value.var = "mean_score")
            module_dt[, difference := get(lev[1]) - get(lev[2])]
            module_dt = module_dt[, c(1,4,2,3)]
            lev = c(lev, "difference")
            
        }else{
            stop("Unrecognized selClusterMethod: ", sort_method)
        }
        
        # add extra info
        if(module_resolution == "coarse"){
            module_dt = merge(module_dt, coarse_descriptions_dt, by = "module")        
        }else if(module_resolution == "fine"){
            module_dt = merge(module_dt, fine_descriptions_dt, by = "module")        
        }
        # sort
        module_dt = module_dt[order(module_dt[, 2][[1]], decreasing = TRUE),]
        module_dtR(module_dt)
        DT::datatable(module_dt, selection = "single") %>%
            DT::formatRound(lev, 3)
    })
    
    output$moduleLookup = renderUI({
        goi = req(input$selGeneLookup)
        sel = module_genes_dt[gene_name == goi,]
        stopifnot(nrow(sel) == 1)
        
        if(input$selModuleResolution == "coarse"){
            tags$span(paste0(goi, ": ", paste0("coarse_", sel$Coarse.module)),
                      style="color:blue; font-weight:bold")    
        }else{
            tags$span(paste0(goi, ": ", paste0("fine_", sel$Fine.module)),
                      style="color:blue; font-weight:bold")    
        }
    })
    
    observeEvent({
        input$coarseTable_rows_selected
    }, {
        sel = req(input$coarseTable_rows_selected)
        sel_dt = req(module_dtR())[sel,]
        stopifnot(nrow(sel_dt) == 1)
        updateTextInput(session = session, "txtModuleQuery", value = sel_dt$module)
        # showNotification(paste("selected", sel_dt$module))
    })
    
    isValidModule = reactiveVal(FALSE)
    
    observeEvent({
        input$txtModuleQuery
        meta_dtR()
    },{
        meta_dt = req(meta_dtR())
        module_query = req(input$txtModuleQuery)
        if(module_query %in% colnames(meta_dt)){
            shinyjs::hide(id = "msgBadModule")
            shinyjs::show(id = "msgGoodModule")
            shinyjs::enable(id = "btnViewModuleDetail")
            isValidModule(TRUE)
        }else{
            shinyjs::hide(id = "msgGoodModule")
            shinyjs::show(id = "msgBadModule")
            shinyjs::disable(id = "btnViewModuleDetail")
            isValidModule(FALSE)
        }
    })
    
    output$plotMiniClusters = renderPlot({
        meta_dt = req(meta_dtR())
        sel_set = req(input$selClusterSet)
        sel_id = req(input$selClusterIDs)
        
        k = as.character(meta_dt[[sel_set]]) %in% sel_id
        if(!any(k)) req(NULL)
        p = ggplot(data = NULL, aes_string(x = "UMAP_1", y = "UMAP_2", color = sel_set)) + 
            geom_point(data = meta_dt[!k,], color = "gray10", size = .1) +
            geom_point(data = meta_dt[k,], size = .2) +
            labs(color = '', x= "", y = "") +
            guides(color = guide_legend(override.aes = aes(size = 3), nrow = 3)) +
            theme_cowplot() +
            theme(legend.position = "bottom", 
                  axis.text = element_blank(), 
                  axis.ticks = element_blank(), 
                  axis.line = element_blank())
        if(input$checkSplitGenotype){
            p = p + facet_wrap(~orig.ident) 
        }
        p
    })
    
    output$plotMiniModule = renderPlot({
        meta_dt = req(meta_dtR())
        sel_set = req(input$selClusterSet)
        sel_id = req(input$selClusterIDs)
        mod_query = input$txtModuleQuery
        if(is.null(mod_query)) req(NULL)
        if(isValidModule()){
            p = ggplot(data = NULL, aes_string(x = "UMAP_1", y = "UMAP_2", color = mod_query)) + 
                geom_point(data = meta_dt, size = .2) +
                labs(color = '', x= "", y = "") +
                scale_color_viridis_c(option = "magma") +
                theme_cowplot() +
                theme(legend.position = "bottom", 
                      axis.text = element_blank(), 
                      axis.ticks = element_blank(), 
                      axis.line = element_blank(), 
                      legend.text = element_text(size = 6))
            if(input$checkSplitGenotype){
                p = p + facet_wrap(~orig.ident) 
            }
        }else{
            p = ggplot() + 
                annotate("text", x=0, y = 0, label = 'waiting for valid module') + 
                theme_void()
        }
        
        # ggsave("tmp.png", p)
        p
    })
    
    geneDetailDT = reactiveVal()
    
    output$dtGeneDetail = DT::renderDataTable({
        meta_dt = req(meta_dtR())
        sel_set = req(input$selClusterSet)
        sel_res = req(input$selModuleResolution)
        sel_mod = req(input$txtModuleQuery)
        if(sel_res == "fine"){
            gois = fine_genes_dt[module == sel_mod]$gene_name
        }else{
            gois = coarse_genes_dt[module == sel_mod]$gene_name
        }
        missed = setdiff(gois, rownames(pbmc_data()))
        rna_dt = get_rna_dt(pbmc_data(), intersect(gois, rownames(pbmc_data())))
        rna_dt = rna_dt[, .(mean_expression = mean(expression)), .(gene_name)]
        rna_dt = rna_dt[order(mean_expression, decreasing = TRUE)]
        rna_dt = rbind(rna_dt, data.table(gene_name = missed, mean_expression = NA))
        geneDetailDT(rna_dt)
        DT::datatable(rna_dt, selection = list(mode = 'single', selected = c(1), target = 'row')) %>%
            DT::formatRound("mean_expression", 3)
    })
    
    output$plotModuleDetail = renderPlot({
        meta_dt = req(meta_dtR())
        sel_set = req(input$selClusterSet)
        sel_id = req(input$selClusterIDs)
        mod_query = input$txtModuleQuery
        # k = as.character(meta_dt[[sel_set]]) %in% sel_id
        p = ggplot(data = NULL, aes_string(x = "UMAP_1", y = "UMAP_2", color = mod_query)) + 
            # geom_point(data = meta_dt[!k,], color = "gray10", size = .1) +
            # geom_point(data = meta_dt[k,], size = .2) +
            geom_point(data = meta_dt, size = .2) +
            labs(color = '', x= "", y = "") +
            scale_color_viridis_c(option = "magma") +
            # guides(color = guide_legend(override.aes = aes(size = 3), nrow = 3)) +
            theme_cowplot() +
            theme(legend.position = "bottom", 
                  axis.text = element_blank(), 
                  axis.ticks = element_blank(), 
                  axis.line = element_blank(),
                  legend.text = element_text(size = 6)) +
            labs(title = mod_query)
        if(input$checkSplitGenotype){
            p = p + facet_wrap(~orig.ident) 
        }
        p
    })
    
    output$plotGeneDetail = renderPlot({
        meta_dt = req(meta_dtR())
        sel_set = req(input$selClusterSet)
        sel_id = req(input$selClusterIDs)
        sel_dt = req(geneDetailDT())
        if(nrow(sel_dt) == 0) req(NULL)
        
        sel_gene = sel_dt$gene_name[req(input$dtGeneDetail_rows_selected)]
        if(sel_gene %in% rownames(pbmc_data())){
            rna_dt = get_rna_dt(pbmc_data(), sel_gene)
            rna_dt = merge(meta_dt[, .(id, UMAP_1, UMAP_2, orig.ident)], rna_dt, by = "id")
            p = ggplot(data = rna_dt, aes_string(x = "UMAP_1", y = "UMAP_2", color = "expression")) + 
                geom_point(size = .2) +
                labs(color = '', x= "", y = "") +
                scale_color_viridis_c() +
                # guides(color = guide_legend(override.aes = aes(size = 3), nrow = 3)) +
                theme_cowplot() +
                theme(legend.position = "bottom", 
                      axis.text = element_blank(), 
                      axis.ticks = element_blank(), 
                      axis.line = element_blank()) +
                labs(title = sel_gene)
            if(input$checkSplitGenotype){
                p = p + facet_wrap(~orig.ident) 
            }
        }else{
            p = ggplot() + annotate("text", x = .5, y = .5, label = "gene not found in scRNA data") +
                theme_void()
        }
        p
    })
    
    observeEvent({
        input$btnViewModuleDetail
    }, {
        showModal(modalDialog(
            fluidRow(
                column(width = 12,
                       withSpinner(uiOutput("infoModuleDetail"))
                )
            ),
            fluidRow(
                column(width = 12,
                       withSpinner(plotOutput("plotModuleDetail", width = "560px", height = "400px")),
                       tippy_this("plotModuleDetail", "Module score of selected Immgen Module in UMAP.")
                       
                )
            ),
            fluidRow(
                column(width = 12,
                       div(withSpinner(DT::dataTableOutput("dtGeneDetail")), style = "font-size: 80%; width: 80%"),
                       tippy_this("dtGeneDetail", "Genes in selected Immgen Module ranked by average scaled expression.")
                )
            ),
            fluidRow(
                column(width = 12,
                       withSpinner(plotOutput("plotGeneDetail", width = "560px", height = "400px")),
                       tippy_this("plotGeneDetail", "Scaled expression value of selected gene in UMAP.")
                       
                )
            ),
            actionButton("btnModalDetailDone", "Done"),
            # sidebarLayout(
            #     sidebarPanel = sidebarPanel(
            #         div(withSpinner(DT::dataTableOutput("dtGeneDetail")), style = "font-size: 50%; width: 50%"),
            #         actionButton("btnModalDetailDone", "Done")),
            #     mainPanel = mainPanel(
            #         withSpinner(plotOutput("plotModuleDetail")),
            #         withSpinner(plotOutput("plotGeneDetail"))
            #     )
            # ),
            title = "View module details",
            size = "l",
            easyClose = TRUE,
            footer = NULL
        ))
    })
    
    output$infoModuleDetail = renderUI({
        sel = req(input$coarseTable_rows_selected)
        sel_dt = req(module_dtR())[sel,]
        if(input$selModuleResolution == "coarse"){
            list(
                tags$h6(paste("ID:", sel_dt$module)),
                tags$h6(paste("Title:", sel_dt$title)),
                tags$h6(paste("Induction Pattern:", sel_dt$Induction.Pattern)),
                tags$h6(paste("Pattern Type:", sel_dt$Pattern.type)),
                tags$h6(paste("Population:", sel_dt$Population)),
                tags$h6(paste("Gene Families:", sel_dt$families))
            )    
        }else{
            list(
                tags$h6(paste("ID:", sel_dt$module)),
                tags$h6(paste("Parent:", paste0("coarse_", sel_dt$Coarse.module))),
                tags$h6(paste("Title:", sel_dt$Title)),
                tags$h6(paste("Induction Pattern:", sel_dt$Induction.Pattern)),
                tags$h6(paste("Gene Families:", sel_dt$families)),
                tags$h6(paste("MSIgDB:", sel_dt$MSIgDB)),
                tags$h6(paste("ExpressioninGeneAtlas:", sel_dt$ExpressioninGeneAtlas)),
                tags$h6(paste("TopRegulators:", sel_dt$TopRegulators)),
                tags$h6(paste("BS_binding:", sel_dt$BS_binding)),
                tags$h6(paste("Remarks:", sel_dt$Remarks))
            )  
        }
        
            
    })
    
    observeEvent({
        input$btnModalDetailDone
    }, {
        removeModal()
    })
}

# Run the application 
shinyApp(ui = ui, server = server)

