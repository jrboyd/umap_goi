library(data.table)
library(seqsetvis)
library(magrittr)
library(GenomicRanges)
library(BiocFileCache)
library(Seurat)


plot0 = function(width = 1, height = 1){
    fudge = 0.037037037037
    plot(c(0+fudge*width, width-fudge * width), c(0+fudge*height, height-fudge * height), type = 'n', xlab = '', ylab = '', axes = F, )
}

startRasterMode = function(width = 1, height = 1){
    png('tmp.png', units = 'in', width = width, height = height, res = 600)
    par(mai = rep(0,4))
}

stopRasterMode = function(mai = NA){
    dev.off()
    if(length(mai) > 1 && !is.na(mai)){
        par(mai = mai)
    }
    plot0()
    rasterImage(png::readPNG("tmp.png", native = FALSE), xleft = 0, xright = 1, ytop = 1, ybottom = 0,
                interpolate = FALSE)
}


get_cluster_rename = function(){

    rename_clust = c("main" = "1", 
                     "10" = "2", 
                     "11" = "3", 
                     "7" = "4", 
                     "4" = "5", 
                     "13" = "6", 
                     "6"= "7", 
                     "15" = "8",
                     "12" = "9", 
                     "9" = "10", 
                     "8" = "11", 
                     "14" = "12"
    )
    rename_clust
}

get_clusters_to_combine = function(){
    all_clust = 0:16
    main_group = c(0,1,2,3,5,16)
    other_groups = as.list(setdiff(all_clust, main_group))
    names(other_groups) = setdiff(all_clust, main_group)
    to_comb = c(list(main = main_group), other_groups)
    
    names(to_comb)
    names(to_comb) = get_cluster_rename()[names(to_comb)]
    to_comb
}

get_cluster_to_combine.old = function(){
    list("west_island" = c(15, 6), 
         "west_projection" = c(13, 4), 
         "lonely_island" = c(11), 
         "north_archipelago" = c(9, 8, 14, 12), 
         "main_continent" = c(7, 0, 2, 16, 2, 1, 3, 5, 10))
}

get_clusters_colors = function(){
    col.blue = safeBrew(7, "blues")[2:5+2]
    names(col.blue) = c("main", 11, 10, 7)
    col.purp = safeBrew(5, "BuPu")[c(3, 5)]
    names(col.purp) = c(4, 13)
    col.oranges = safeBrew(6, "Oranges")[c(3, 6)]
    names(col.oranges) = c(6, 15)
    col.reds = safeBrew(5, "Reds")[2:5]
    names(col.reds) = c(8, 9, 12, 14)
    col.all = c(col.blue, col.purp, col.oranges, col.reds)
    
    names(col.all) = get_cluster_rename()[names(col.all)]
    col.all
}

#' get_meta_dt
#'
#' @param seurat_obj an object of class Seurat to extract meta data from
#' @param to_combine 
#' @param cluster_vars cluster variables, will be made factors with levels ordered by decreasing size.
#' @param manual_factors named list, names specify variables and entries specify factor levels to apply.
#' @param reduction string to specify reduction, must be in seurat_obj@reductions
#'
#' @return
#' @export
#'
#' @examples
get_meta_dt = function(seurat_obj, 
                       combine_source = "seurat_clusters",
                       to_combine = get_clusters_to_combine(),
                       cluster_vars = c("seurat_clusters", "meta_cluster"),
                       manual_factors = list(orig.ident = c("wt", "df4")),
                       reduction = "umap"){
    meta_dt = as.data.table(seurat_obj@meta.data)
    meta_dt$id = seurat_obj@meta.data %>% rownames
    
    if(reduction %in% names(seurat_obj@reductions)){
        stopifnot(reduction %in% names(seurat_obj@reductions))
        umap_dt = as.data.table(seurat_obj@reductions[[reduction]]@cell.embeddings, keep.rownames = TRUE) %>% setnames(., "rn", "id")
        meta_dt = merge(meta_dt, umap_dt, by = "id")    
    }else{
        meta_dt$UMAP_1 = NA
        meta_dt$UMAP_2 = NA
    }
    
    if(!is.null(to_combine)){
        to_combine_dt = lapply(to_combine, function(x){
            data.table(indi_clust = factor(x, levels = levels(meta_dt[[combine_source]])))
        }) %>% rbindlist(idcol = "meta_cluster")
        setnames(to_combine_dt, "indi_clust", combine_source)
        meta_dt = merge(meta_dt, to_combine_dt, by = combine_source, allow.cartesian = TRUE) %>% unique    
    }
    
    for(i in seq_along(manual_factors)){
        meta_dt[[names(manual_factors)[i]]] = factor(meta_dt[[names(manual_factors)[i]]], levels = manual_factors[[i]])    
    }
    
    
    stopifnot(cluster_vars %in% colnames(meta_dt))
    for(cv in cluster_vars){
        lev = meta_dt[, .N, c(cv)][order(N, decreasing = TRUE)][[cv]] %>% as.character
        meta_dt[[cv]] = factor(meta_dt[[cv]], levels = lev)
    }
    
    if(!is.null(meta_dt$meta_cluster)){
        meta_dt$meta_cluster = factor(meta_dt$meta_cluster, levels = get_cluster_rename())    
    }
    
    
    meta_dt
}

get_meta_dt.raw = function(seurat_obj, manual_factors = list(orig.ident = c("wt", "df4"))){
    get_meta_dt(seurat_obj, to_combine = NULL, cluster_vars = character(), manual_factors = manual_factors)    
}


get_rna_dt = function(seurat_obj, sel_genes = NULL, raw_counts = FALSE, assay_name = "RNA"){
    if(is.null(sel_genes)){
        rna_dt = as.data.frame(seurat_obj@assays[[assay_name]]@counts)
    }else{
        len_input = length(sel_genes)
        sel_genes = intersect(sel_genes, seurat_obj@assays[[assay_name]] %>% rownames)
        if(length(sel_genes) != len_input){
            d = len_input - length(sel_genes)
            perc = round(100 * d / len_input, 2)
            warning(perc, "% (", d, " of ", len_input, ") genes discarded to match scRNAseq")
        }
        if(!raw_counts){
            rna_dt = as.data.frame(seurat_obj@assays[[assay_name]][sel_genes, ])    
        }else{
            rna_dt = as.data.frame(seurat_obj@assays[[assay_name]]@counts[sel_genes, ])    
        }
        
        
    }
    rna_dt$gene_name = rownames(rna_dt)
    rna_dt = melt(as.data.table(rna_dt), variable.name = "id", value.name = "expression", id.vars = "gene_name")
    rna_dt
}

my_ggsave = function(file, p, width = dev.size()[1], height = dev.size()[2]){
    file = basename(file)
    file = sub(".png$", "", file)
    file = sub(".pdf$", "", file)
    file = res_file(file)
    ggsave(paste0(file, ".png"), p, width = width, height = height)
    ggsave(paste0(file, ".pdf"), p, width = width, height = height, device = cairo_pdf)
}

my_saveWidget = function(p, file){
    file = basename(file)
    htmlwidgets::saveWidget(plotly::ggplotly(p), file, selfcontained = TRUE)
    file.copy(file, res_file(file), overwrite = TRUE)
    file.remove(file)
}

plot_module_scores = function(seu_obj, module_genes, cluster_var = "seurat_clusters", facet_var = "orig.ident"){
    seu_obj = readRDS("Bcell_combined.refixed_120919.rds")
    module_genes = fread("Bcell_inactive_UP_DEG.txt")[[1]]
    cluster_var = "seurat_clusters"
    mod_res = Seurat::AddModuleScore(seu_obj, list(module_genes), name = "myscore")
    meta_dt = get_meta_dt(mod_res)
    rna_dt = get_rna_dt(seu_obj, module_genes)
    
    ggplot(rna_dt[gene_name %in% module_genes[1:6]], aes(x = expression)) + geom_density() + facet_wrap(~gene_name)
    fraction_dt = rna_dt[, .(fraction = sum(expression > 0)/.N) , .(gene_name)][order(fraction, decreasing = TRUE)]
    fraction_dt$fraction %>% plot
    
    if(!is.null(meta_dt$module_score)) meta_dt$module_score = NULL
    setnames(meta_dt, "myscore1", "module_score")
    p_umap = ggplot(meta_dt, aes(x = UMAP_1, y = UMAP_2, color = module_score)) +
        geom_point(size = .3) +
        scale_color_viridis_c(option = "magma") +
        theme(panel.background = element_rect(fill = "gray50"))
    if(facet_var != "" && is.character(facet_var)){
        p_umap + facet_wrap(paste0(".~", facet_var))
    }
    
    p_bar = ggplot(meta_dt, aes_string(x = cluster_var, fill = facet_var)) + 
        geom_bar(position = "dodge")
    
    p_violin = ggplot(meta_dt, aes_string(x = cluster_var, y = "module_score", fill = facet_var)) + 
        geom_violin() 
    return(list(umap = p_umap, bar= p_bar, violin = p_violin))
}

gene_name_gs2mm = function(x){
    x = tolower(x)
    substr(x, 1, 1) = toupper(substr(x, 1, 1))
    x = sub("(?<=[0-9])rik$", "Rik", x, perl = TRUE)
    is_rik = grepl("Rik$", x)
    substr(x[is_rik], 8, 8) = toupper(substr(x[is_rik], 8, 8))
    x
}


my_harmonize = function(seurat_obj, imm_mat, gois, cluster_var = "seurat_clusters", meta_dt = NULL, 
                        # row_agg_FUN = function(x){(x + .1)/mean(x + .1)}
                        row_agg_FUN = function(x){x}
){
    DefaultAssay(seurat_obj) = "RNA"
    sc_dt = get_rna_dt(seurat_obj, gois)
    if(is.null(meta_dt)){
        meta_dt = get_meta_dt(seurat_obj)    
    }
    
    sc_dt = merge(sc_dt, meta_dt[, .(id, orig.ident, grouping = get(cluster_var))], by = "id")
    
    sc_wide = dcast(sc_dt, gene_name~grouping, fun.aggregate = mean, value.var = "expression")
    sc_mat = as.matrix(sc_wide[,-1])
    rownames(sc_mat) = sc_wide$gene_name
    sc_mat = sc_mat[rowSums(sc_mat) > 0,]
    sc_mat.norm =t(apply(sc_mat, 1, row_agg_FUN))
    
    imm_mat =  imm_mat[rownames(sc_mat.norm), ]
    imm_mat.norm =t(apply(imm_mat, 1, row_agg_FUN))
    stopifnot(all(rownames(sc_mat) == rownames(imm_mat)))
    
    harmony_mat = limma::normalizeQuantiles(cbind(sc_mat.norm, imm_mat.norm))
    f = tempfile()
    png(f)
    hres_harmony = double_clustering(harmony_mat)
    dev.off()
    file.remove(f)
    harmony_dt = my_melt(harmony_mat, 
                         rownames.var = "gene_name", 
                         colnames.var = "group", 
                         value.var = "expression", 
                         rownames.order = colnames(hres_harmony$carpet),
                         colnames.order = rownames(hres_harmony$carpet)
    )
    return(list(data.table = harmony_dt, matrix = harmony_mat, cluster_var = cluster_var, heatmap2_res = hres_harmony))
}

my_melt = function(mat, rownames.var = "row", colnames.var = "column", value.var = "value", rownames.order = NULL, colnames.order = NULL){
    if(is.null(rownames.order)){
        rownames.order = rownames(mat)
    }
    if(is.null(colnames.order)){
        colnames.order = colnames(mat)
    }
    dt = melt(as.data.table(mat, keep.rownames = TRUE), id.vars = "rn")
    setnames(dt, c("rn", "variable", "value"), c(rownames.var, colnames.var, value.var)) 
    dt[[rownames.var]] = factor(dt[[rownames.var]], levels = rownames.order)
    dt[[colnames.var]] = factor(dt[[colnames.var]], levels = colnames.order)
    dt
}

symbol2uniprot = function(x, db = org.Mm.eg.db::org.Mm.eg.db){
    bres = clusterProfiler::bitr(x, fromType = "SYMBOL",
                toType = c("UNIPROT"),
                OrgDb = db)
    bres[!duplicated(bres$SYMBOL),]$UNIPROT
    
}

my_cprof = function(gls, bg = NULL, ont = c("BP", "MF", "CC", "Kegg")[1]){
    if(ont == "BP"){
        rname = digest::digest(list(gls, bg, "V1"))    
    }else{
        rname = digest::digest(list(gls, bg, ont, "V1"))    
    }
    
    cprof = ssvRecipes::bfcif(bfc, rname, function(){
        if(is.null(bg)){
            if(ont != "Kegg"){
                clusterProfiler::compareCluster(gls, 
                                                # universe      = unlist(hres@cluster_members),
                                                fun = "enrichGO", 
                                                OrgDb = org.Mm.eg.db::org.Mm.eg.db, 
                                                keyType = "SYMBOL", ont = ont, 
                                                pAdjustMethod = "BH",
                                                pvalueCutoff  = 1,
                                                qvalueCutoff  = 1)  
            }else{
                gls = lapply(gls, symbol2uniprot)
                clusterProfiler::compareCluster(gls, 
                                                # universe      = unlist(hres@cluster_members),
                                                fun = "enrichKEGG", 
                                                organism     = 'mmu',
                                                pAdjustMethod = "BH",
                                                pvalueCutoff  = 1,
                                                qvalueCutoff  = 1)  
            }
            
        }else{
            if(ont != "Kegg"){
                clusterProfiler::compareCluster(gls, 
                                                bg,
                                                fun = "enrichGO", 
                                                OrgDb = org.Mm.eg.db::org.Mm.eg.db, 
                                                keyType = "SYMBOL", ont = ont, 
                                                pAdjustMethod = "BH",
                                                pvalueCutoff  = 1,
                                                qvalueCutoff  = 1)  
            }else{
                gls = lapply(gls, symbol2uniprot)
                bg = symbol2uniprot(bg)
                clusterProfiler::compareCluster(gls, 
                                                bg,
                                                fun = "enrichKEGG", 
                                                organism     = 'mmu',
                                                pAdjustMethod = "BH",
                                                pvalueCutoff  = 1,
                                                qvalueCutoff  = 1)  
            }
            
        }
        
    })
    raw = cprof
    cprof@compareClusterResult = subset(cprof@compareClusterResult, p.adjust < .05)
    
    cprof = ssvRecipes::bfcif(bfc, digest::digest(list(cprof, "V1")), function(){
        clusterProfiler::simplify(cprof)
    })
    p = clusterProfiler::dotplot(cprof)
    plot(p)
    invisible(list(plot = p, raw = raw, final = cprof)) 
}

setup_IK_np = function(){
    ik_peak_files = dir("~/R/SF_AutoImmune_ssv/data/chromatin/IKAROS/peaks/macs2/", full.names = TRUE)
    ik_gr = easyLoad_narrowPeak(ik_peak_files)
    names(ik_gr) = sub("_peaks.narrowPeak", "", names(ik_gr))
    k = grepl("IK[1-3]$", names(ik_gr))
    ik_gr = ik_gr[k]
    ik_olaps = ssvOverlapIntervalSets(ik_gr)
    p_peak_venn = ssvFeatureVenn(ik_olaps) + labs(title = "Ikaros peak call overlap")
    p_peak_olap = ssvFeatureBinaryHeatmap(ik_olaps, raster_approximation = TRUE)
    pg_olap = cowplot::plot_grid(p_peak_venn, p_peak_olap, rel_widths = c(1.5, 1))
    length(ik_olaps)
    qgr = resize(ik_olaps, 6e3, fix = "center")
    names(qgr) = paste0("region_", names(qgr))
    
    return(list(ikaros_np_gr = ik_gr, ikaros_np_overlaps = ik_olaps, ikaros_query_gr = qgr,
                plot_peak_venn = p_peak_venn, plot_peak_beatmap = p_peak_olap))
}

setup_IK_bams = function(){
    ik_bam_files = dir("~/R/SF_AutoImmune_ssv/data/chromatin/IKAROS/bam/", pattern = ".bam$", full.names = TRUE)
    qdt = data.table(file = ik_bam_files)
    qdt[, sample := basename(file) %>% sub("\\..+", "", .)]
    qdt[, sample := sub("Ikaros", "IK", sample)]
    qdt[, sample := sub("mature", "mat", sample)]  
    return(list(query_bam_table = qdt))
}

# setup_HM_np = function(){
#     np_files = dir("data/chromatin/histones/peaks/", full.names = TRUE)
# }

setup_HM_bw = function(){
    bw_files = dir("~/R/SF_AutoImmune_ssv/data/chromatin/histones/FE_bigwig/", full.names = TRUE, pattern = ".bw$")
    bw_files = bw_files[!dir.exists(bw_files)]
    
    bw_files2 = dir("~/R/SF_AutoImmune_ssv/data/chromatin/histones/FE_bigwig/published/", full.names = TRUE, pattern = ".bw$")
    bw_files3 = c("~/R/SF_AutoImmune_ssv/data/chromatin/IKAROS/FE_bigwig/WT_unstim_IK_FE.bw")
    
    qdt = data.table(file = c(bw_files, bw_files2, bw_files3))
    # qdt = qdt[normalizePath(file) != 
    #               normalizePath("~/R/SF_AutoImmune_ssv/data/chromatin/histones/FE_bigwig/published//Ikaros_proB_FE.bw")]#bad file
    
    qdt[, source := ifelse(grepl("published", file), "public", "inhouse")]
    
    qdt[, status := ifelse(grepl("[dD]F4", file), "DF4", "WT")]
    qdt[, stage := ifelse(grepl("proB", file), "proB", "matB")]
    qdt[, mark := tstrsplit(basename(file), split = "_", keep = 1)]
    qdt[!grepl("(^Ik)|(H3)", mark), mark := tstrsplit(basename(file), split = "_", keep = 2)]
    qdt[mark == "H34me3", mark := "H3K4me3"]
    # qdt[, .N, .(status, stage, mark)]
    
    qdt[, status := toupper(status)]
    qdt[, mark := toupper(mark)]
    qdt[, treatment := "0"]
    qdt[, rep := "rep1"]
    qdt[grepl("rep2", file), rep := "rep2"]
    
    qdt[file == "~/R/SF_AutoImmune_ssv/data/chromatin/IKAROS/FE_bigwig/WT_unstim_IK_FE.bw", mark := "IKAROS"]
    
    qdt[, name := paste(source, stage, status, treatment, mark, rep, sep = "\n")]
    
    qdt$stage = factor(qdt$stage, levels = c("proB", "matB"))
    qdt$status = factor(qdt$status, levels = c("WT", "DF4"))
    qdt$treatment = factor(qdt$treatment, levels = c("0", "6", "M", "M40"))
    qdt$mark = factor(qdt$mark, levels = c("H3K4ME3", "H3K27AC", "H3K9AC", "H3K27ME3", 'ATAC', "IKAROS"))
    qdt = qdt[order(rep)][order(treatment)][order(mark)][order(status)][order(stage)]
    qdt$name = factor(qdt$name, levels = qdt$name %>% unique)
    qdt[mark == "IKAROS" & stage == "matB", file := "/slipstream/home/sethandy/project_AI/chromatin/IKAROS/FE_bigwig/IK_merge_FE.bw"]
    
    
    qdt[is.na(mark) & grepl("K4", file), mark := "H3K4ME3"]
    
    qdt[is.na(mark) & grepl("27ac", file), mark := "H3K27AC"]
    
    qdt[grepl("high", file), rep := paste(rep, "high")]
    qdt[grepl("low", file), rep := paste(rep, "low")]
    
    stopifnot(any(!duplicated(qdt$name)))
    stopifnot(any(!is.na(qdt$mark)))
    stopifnot(any(!is.na(qdt$rep)))
    stopifnot(any(!is.na(qdt$treatment)))
    return(qdt[])
}

setup_ATAC_bw = function(){
    # bw_files = dir("~/R/SF_AutoImmune_ssv/data/chromatin/ATACseq/round2_ATAC/", full.names = TRUE)
    bw_files = dir("/slipstream/home/sethandy/project_AI/atac/merged_atac", full.names = TRUE, pattern = "rpkm.+bw")
    
    qdt = data.table(file = bw_files)
    qdt[, c("status", "treatment", "rep") := tstrsplit(basename(file), split = "[_\\.]", keep = 1:3)] 
    qdt$mark = "ATAC"
    qdt[, status := toupper(status)]
    qdt[, mark := toupper(mark)]
    qdt[ rep == "bw", rep := treatment]
    qdt[, rep := toupper(rep)]
    qdt[grepl("^r[12]$", treatment), treatment := sub("(DF4)|(WT)", "", status)]
    # qdt[treatment == "0", treatment := "t0"]
    qdt[treatment == "U", treatment := "6"]
    qdt[grepl("DF4", status), status := "DF4"]
    qdt[grepl("WT", status), status := "WT"]
    
    qdt[, name := paste(status, treatment, mark, rep, sep = "\n")]
    
    qdt$status = factor(qdt$status, levels = c("WT", "DF4"))
    qdt$treatment = factor(qdt$treatment, levels = c("0", "6", "M", "M40"))
    qdt$mark = factor(qdt$mark, levels = c("H3K9AC", "H3K27ME3", 'ATAC'))
    qdt = qdt[order(rep)][order(treatment)][order(mark)][order(status)]
    qdt$name = factor(qdt$name, levels = qdt$name %>% unique)
    qdt[grepl("merge", file), rep := "merged"]
    qdt[, treatment := "naive"]
    return(list(query_bw_table = qdt))
}

library(openxlsx)
setup_gene_lists = function(src_file = "~/R/SF_AutoImmune_ssv/Bcell_relevant_gene_lists.xlsx"){
    sn = openxlsx::getSheetNames(src_file)
    bcell_gl = lapply(sn, function(s){
        df = read.xlsx(src_file, sheet = s)
        out = toupper(df[[1]])
        out = gsub(" ", "", out)
        out
    })
    
    names(bcell_gl) = gsub(" ", "_", sn)
    
    RTK_module_259 = fread("~/R/SF_AutoImmune_ssv/GSEA_MODULE_259.RTK_signaling.txt")[[1]][-1]
    
    bcell_gl$RTK_module_259 = RTK_module_259
    
    surface_markers = openxlsx::read.xlsx("~/R/SF_AutoImmune_ssv/ENTREZ_gene_name_CFSA_validated_surface_proteins_121219.xlsx", colNames = FALSE)[[1]]
    tf = openxlsx::read.xlsx("~/R/SF_AutoImmune_ssv/gene_lists/complete_list_TF_atlas.xlsx", sheet = 2)[[1]]
    
    bcell_gl[["Surface_Markers"]] = surface_markers
    bcell_gl[["Transcription_Factors"]] = tf
    
    return(bcell_gl)    
}

setup_signature_lists = function(src_file = "supplementaryTable4_MOESM4_ESM.xlsx"){
    sn = openxlsx::getSheetNames(src_file)
    # bcell_gl = lapply(sn, function(s){
    df = read.xlsx(src_file, sheet = sn[2], startRow = 2)
    df = subset(df, !is.na(mouse.gene))
    # df$`Supplemental.Table.4..Top.20.marker.genes.for.each.single-cell.RNA-seq.cluster.`
    valid = which(df$cluster != "NA")
    froms = valid
    tos = c(froms[-1]-1, nrow(df))
    df$cluster[froms]
    
    for(i in seq_along(froms)){
        df$cluster[seq(froms[i]+1, tos[i])] = df$cluster[froms[i]]
    }
    
    
    gl = split(df$mouse.gene, df$cluster)
    names(gl) = sub(" \\(.+", "", names(gl))
    gl
}

setup_surface_markers = function(src_file = "ENTREZ_gene_name_CFSA_validated_surface_proteins_121219.xlsx"){
    gois = read.xlsx(src_file)[[1]]
    return(gois)
}

setup_smale_list = function(fc_cutoff = 3){
    smale_de = openxlsx::read.xlsx("13072_2015_12_MOESM3_ESM.xlsx", startRow = 4)
    smale_de = subset(smale_de, grepl("BCR", cluster))
    smale_de$fc = log2(smale_de$`bcr/rest`)
    smale_de = subset(smale_de, abs(fc) >= fc_cutoff)
    smale_de$group = ifelse(smale_de$fc > 0, "bcr high / rest low", "bcr low / rest high")
    
    split(smale_de$gene_name, smale_de$group)
}

setup_volcano_genes = function(files = dir("gene_lists", full.names = TRUE)){
    gl = lapply(files, function(f){
        fread(f, header = FALSE)[[1]]
    })
    names(gl) = basename(files) %>% sub("\\..+", "", .)
    
    lapply(gl, gene_name_gs2mm)
}

setup_DE_results = function(dir_DE = "~/R/SF_AutoImmune_ssv/data/RNAseq/DESeq2_outputs"){
    de_files = dir(dir_DE, pattern = "^[(DF)|(WT)].+vs", full.names = TRUE)
    names(de_files) = basename(de_files)
    is_by_status = grepl("WTvsDF4", de_files)
    
    dt_by_status = lapply(de_files[is_by_status], function(x){
        dt = fread(x)
        dt[, file := basename(x)]
        val_cn = setdiff(colnames(dt), c("ID", "padj", "log2FoldChange", "pvalue", "foldChange", "log10padj", "file"))
        dt$baseMean = rowMeans(dt[, val_cn, with = FALSE])
        
        dt[, c("from", "to", "treatment") := tstrsplit(file, "(vs)|(_)", keep = 1:3)]
        dt[treatment == "N", treatment := "6"]
        dt[, vs := paste0(to, "vs", from, "_", treatment)]
        
        dt[, .(file, gene_name = ID, baseMean, log2FoldChange, pvalue, padj, vs, to, from, treatment)]
        
    }) %>% rbindlist
    dt_by_treatment = lapply(de_files[!is_by_status], function(x){
        dt = fread(x)
        dt[, file := basename(x)]
        val_cn = setdiff(colnames(dt), c("ID", "padj", "log2FoldChange", "pvalue", "foldChange", "log10padj", "file"))
        dt$baseMean = rowMeans(dt[, val_cn, with = FALSE])
        
        dt[, c("status", "from", "to") := tstrsplit(file, "(vs)|(_)", keep = 1:3)]
        dt[, vs := paste0(status, "_", to, "vs", from)]
        dt[from == "N", from := "U"]
        dt[to == "N", to := "U"]
        dt[status == "DF4", status := "DF"]
        dt[, .(file, gene_name = ID, baseMean, log2FoldChange, pvalue, padj, vs, status, to, from)]
        
    }) %>% rbindlist
    
    de_prob = setup_DE_results.proB()
    dt_by_status = rbind(dt_by_status, de_prob)
    
    return(list(DE_by_treatment = dt_by_treatment, DE_by_status = dt_by_status))
}

setup_DE_results.mirrored = function(){
    de_res = setup_DE_results()
    de_res$DE_by_treatment = rbind(de_res$DE_by_treatment, 
                                   de_res$DE_by_treatment[, .(file, gene_name, baseMean, log2FoldChange = -log2FoldChange, pvalue, padj, vs, status, to = from, from = to)]
    )
    de_res$DE_by_status = rbind(de_res$DE_by_status,
                                de_res$DE_by_status[, .(file, gene_name, baseMean, log2FoldChange = -log2FoldChange, pvalue, padj, vs, to = from, from = to, treatment)])
    de_res
}

setup_DE_results.proB = function(de_prob_file = "~/R/SF_AutoImmune_ssv/data/RNAseq/DESeq2_outputs/proB_wtvsdf4_RLE_normalized_DEGs.csv"){
    
    de_dt = fread(de_prob_file)
    de_dt$file = basename(de_prob_file)
    de_dt = de_dt[, .(file, 
                      gene_name = ID, 
                      baseMean = ( proB_wt_rep1 + proB_wt_rep2 + proB_df4_rep1 + proB_df4_rep2)/4, 
                      log2FoldChange, 
                      pvalue, 
                      padj, 
                      vs = file, 
                      to = "DF4", 
                      from = "WT", 
                      treatment = "proB")]
    de_dt
}

setup_DE_results.old = function(de_dir = "~/R/SF_AutoImmune_ssv/data/RNAseq/DESeq2_outputs/old.with_WT_M40_rep1/"){
    warning("deprecacted")
    de_files = dir(de_dir, full.names = TRUE, pattern = "^deseq")
    # de_files = de_files[!grepl("batch", de_files)]
    
    is_treatment_comp = grepl("(DF_MvsM40)|(DF_UvM40)|(DF_UvM)|(WT_MvM40)|(WT_UvM40)|(WT_UvM)", de_files)
    k = is_treatment_comp
    de_files[is_treatment_comp]
    de_files[!is_treatment_comp]
    
    de_dt = lapply(de_files[k], fread)
    de_dt = lapply(de_dt, function(x)x[,1:7])
    names(de_dt) = basename(de_files[k])
    de_dt = rbindlist(de_dt, use.names = TRUE, idcol = "file")
    de_dt[, vs := file %>% sub("deseq_results_", "", .) %>% 
              sub(".csv", "", .) %>% sub("_batch_0.05", "", .)]
    de_dt[, c("status", "to", "from") := tstrsplit(vs, "(vs?)|(_)")]
    de_dt[, .N, .(status, from, to)]
    setnames(de_dt, "V1", "gene_name")
    de_dt.trt = de_dt
    
    k = !k
    de_dt = lapply(de_files[k], fread)
    de_dt = lapply(de_dt, function(x)x[,1:7])
    names(de_dt) = basename(de_files[k])
    de_dt = rbindlist(de_dt, use.names = TRUE, idcol = "file")
    de_dt[, vs := file %>% sub("deseq_results_", "", .) %>% 
              sub(".csv", "", .) %>% sub("_batch_0.05", "", .)]
    de_dt$vs %>% unique
    de_dt$to = "WT"
    de_dt$from = "DF4"
    de_dt$treatment = "-"
    de_dt[vs == "WT06vDF06", treatment := "6"]
    de_dt[vs == "WTMvDFM", treatment := "M"]
    de_dt[vs == "WTvDF_M40", treatment := "M40"]
    
    de_dt[, .N, .(treatment, from, to)]
    setnames(de_dt, "V1", "gene_name")
    de_dt.status = de_dt

    de_prob = setup_DE_results.proB()
    
    de_dt.status = rbind(de_dt.status[, colnames(de_prob), with = FALSE], de_prob)
    # message("DE genes loaded as 'de_dt'")
    return(list(DE_by_treatment = de_dt.trt, DE_by_status = de_dt.status))
}

setup_DE_results.mirrored_old = function(){
    de_res = setup_DE_results()
    de_res$DE_by_treatment = rbind(de_res$DE_by_treatment, 
          de_res$DE_by_treatment[, .(file, gene_name, baseMean, log2FoldChange = -log2FoldChange, lfcSE, stat, pvalue, padj, vs, status, to = from, from = to)]
    )
    de_res$DE_by_status = rbind(de_res$DE_by_status,
                                de_res$DE_by_status[, .(file, gene_name, baseMean, log2FoldChange = -log2FoldChange, pvalue, padj, vs, to = from, from = to, treatment)])
    de_res
}

setup_transcript = function(){
    bfc = BiocFileCache()   
    ssvRecipes::bfcif(bfc, "gen_M14_transcript", function(){
        rtracklayer::import.gff("~/gencode.vM14.annotation.gtf.gz", 
                                feature.type = "transcript", format = "gtf")
    })
}

setup_transcript_tss = function(){
    setup_transcript() %>% promoters(., 1, 1)
}

setup_RNA_dt = function(){
    .rna_cnt = setup_RNA_counts()
    
    
    #manually remove oddbal rep
    ko =  "WT_M40_matB_r1"
    cnt_dt = .rna_cnt$tidy
    cnt_dt = cnt_dt[!file %in% ko]
    cnt_mat = .rna_cnt$wide
    cnt_mat[, setdiff(colnames(cnt_mat), ko)]
    norm_factors = edgeR::calcNormFactors(cnt_mat)
    norm_dt = data.table(file = names(norm_factors), norm_factor = norm_factors)
    
    cnt_dt = merge(cnt_dt, norm_dt, by = "file")
    cnt_dt[, norm_count := count / norm_factor]
    cnt_dt[, toss := max(count) < 5, .(gene_name)]
    cnt_dt[, c("status", "treatment", "rep") := tstrsplit(file, "_", keep = c(1, 2, 4))]
    
    zscore = function(x){
        (x - mean(x)) / sd(x)
    }
    
    cnt_dt
    
    cnt_dt[, z := zscore(norm_count), .(gene_name)]
    cnt_dt$status = factor(cnt_dt$status, levels = c("WT", "dF4"))
    table(cnt_dt$treatment)
    
    cnt_dt[treatment == "6", rep := paste0(rep, "b")]
    cnt_dt
    
}

setup_RNA_counts = function(){
    cnt_files = dir("~/R/SF_AutoImmune_ssv/data/RNAseq/counts/", pattern = "counts$", full.names = TRUE)
    cnt_files
    qdt = data.table(file = cnt_files)
    qdt[, c("condition", "treatment", "cell", "rep", "tmp") := tstrsplit(sub(".counts", "", basename(file)), "_")]
    qdt
    k = grepl("(^WT)|(dF4)", qdt$condition)
    k
    qdt = qdt[k,]
    qdt[is.na(rep), rep := "r1"]
    qdt$tmp = NULL
    #remove bad replicate
    qdt = qdt[file != "/slipstream/home/joeboyd/R/SF_AutoImmune_ssv/data/RNAseq/counts//WT_M40_matB.counts"]
    cnt_dt = lapply(qdt$file, function(f){
        dt = fread(f)
        setnames(dt, c("gene_name", "count"))
        dt$file = sub(".counts", "", basename(f))
        dt
    }) %>% rbindlist()
    cnt_dt = cnt_dt[!grepl("^__", gene_name)]
    cnt_dt[!grepl("r[0-3]", file), file := paste0(file, "_r1")]
    cnt_dt.wide = dcast(cnt_dt, gene_name~file, value.var = "count", fill = 0)
    cnt_mat = as.matrix(cnt_dt.wide[,-1])
    rownames(cnt_mat) = cnt_dt.wide$gene_name
    return(list(tidy = cnt_dt, wide = cnt_mat))
}

setup_RNA_counts_with_proB = function(){
    cnt_res = setup_RNA_counts()
    cnt_res$tidy
    proB_cnt = fread("~/R/SF_AutoImmune_ssv/data/RNAseq/BC_remix_norm_counts.csv")
    proB_cnt = melt(proB_cnt, id.vars = "gene_name", variable.name = "file", value.name = "count")
    proB_cnt[, c("stage", "status", "rep") := tstrsplit(file, "_")]
    proB_cnt[, status := toupper(status)]
    proB_cnt[status == "DF4", status := "dF4"]
    proB_cnt[, rep := sub("rep", "r", rep)]
    proB_cnt[, file := paste(status, "0", stage, rep, sep = "_")]
    cnt_dt = rbind(cnt_res$tidy, proB_cnt[, .(gene_name, count, file)])
    cnt_dt[, c("status", "treatment", "stage", 'rep') := tstrsplit(file, "_")]
    table(cnt_dt$stage)
    table(cnt_dt$status)
    table(cnt_dt$treatment)
    table(cnt_dt$rep)
    cnt_dt$status = factor(cnt_dt$status, levels = c("WT", "dF4"))
    cnt_dt$stage = factor(cnt_dt$stage, levels = c("proB", "matB"))
    cnt_dt$treatment = factor(cnt_dt$treatment, levels = c("0", "6", "M", "M40"))
    cnt_dt$rep = factor(cnt_dt$rep, levels = sort(unique(cnt_dt$rep)))
    cnt_dt = cnt_dt[order(rep)][order(status)][order(treatment)][order(stage)]
    cnt_dt$file = factor(cnt_dt$file, levels = unique(cnt_dt$file))
    cnt_dt$file %>% levels
    cnt_res$tidy = cnt_dt
    cnt_res$wide = dcast(cnt_dt, gene_name~file, value.var = "count")
    cnt_res
}



setup_combined_clusters = function(meta_dt, td = list("west_island" = c(15, 6), 
                                                      "west_projection" = c(13, 4), 
                                                      "lonely_island" = c(11), 
                                                      "north_archipelago" = c(9, 8, 14, 12), 
                                                      "main_continent" = c(7, 0, 2, 16, 2, 1, 3, 5, 10))){
    td_dt = lapply(td, function(x){
        data.table(seurat_clusters = factor(x, levels = levels(meta_dt$seurat_clusters)))
    }) %>% rbindlist(idcol = "meta_cluster")
    meta_dt = merge(meta_dt, td_dt, by = 'seurat_clusters', allow.cartesian = TRUE) %>% unique
    meta_dt
}

setup_hematopoetic_markers = function(){
    features = as.list(fread("hematopoesis_markers.csv"))
    features = lapply(features, function(g){
        g = g[g != ""]
        unique(g)
    })
    features
}

setup_PR_gene_lists = function(){
    
    flip_dir = c("up" = "down", "down" = "up")
    files = dir("proB_matB_TSS_analysis/", full.names = TRUE)
    names(files) = basename(files)
    lapply(files, fread, header = FALSE) %>% sapply(., nrow)
    pg_genes_dt = lapply(files, fread, header = FALSE) %>% rbindlist(., idcol = "file") %>% setnames(., "V1", "gene_name")
    pg_genes_dt[, .N, file]
    pg_genes_dt[, c("stage", "direction", "group") := tstrsplit(sub("roB_mat", "roB&mat", file), "[_\\.]", keep = 1:3)]
    pg_genes_dt[, stage := sub("Pro", "pro", stage)]
    pg_genes_dt[, direction := sub("UP", "up", direction)]
    
    pg_genes_dt[file == "matB_UP_proB_down_common.txt", stage := "proB!matB"]
    pg_genes_dt[file == "matB_UP_proB_down_common.txt", group := "common"]
    
    pg_genes_dt[file == "proB_up_matB_down_common.txt", stage := "proB!matB"]
    pg_genes_dt[file == "proB_up_matB_down_common.txt", group := "down"]
    pg_genes_dt[file == "proB_up_matB_down_common.txt", group := "common"]
    
    pg_genes_dt[file == "matB_UP_proB_down_common.txt", direction := flip_dir[direction]]
    
    pg_genes_dt[, name := paste(stage, direction, group, sep = "_")]
    table(pg_genes_dt$name)
    
    table(pg_genes_dt$stage)
    table(pg_genes_dt$direction)
    table(pg_genes_dt$group)
    
    
    gls = split(pg_genes_dt$gene_name, pg_genes_dt$name)
    return(list(dt = pg_genes_dt, list = gls))
}

write_geneList_matrix = function(gl, file){
    if(is.null(names(gl))) stop("gene lists must be named")
    gmat = matrix("", nrow = max(lengths(gl)), ncol = length(gl))
    for(i in seq_along(gl)){
        x = gl[[i]]
        gmat[seq_along(x), i] = as.character(x)
    }
    colnames(gmat) = names(gl)
    write.csv(gmat, file, row.names = FALSE, quote = FALSE)
}


pca_clust = function(pbmc, gl, cells = NULL, min_detect_rate = 0, n_dim = min(10, length(gl)), value.var = c("expression", "zscore")[1]){
    
    rna_dt.full = get_rna_dt(pbmc, gl)
    detect_dt = rna_dt.full[, .(fraction = sum(expression > 0) / .N), .(gene_name)]
    detect_dt$group = "all"
    if(is.null(cells)){
        rna_dt = rna_dt.full#[id %in% names(ps)]
        rna_dt$id = factor(rna_dt$id)
    }else{
        rna_dt = rna_dt.full[id %in% cells]
        rna_dt$id = factor(rna_dt$id, levels = cells)
    }
    
    detect_dt.sel = rna_dt[, .(fraction = sum(expression > 0) / .N), .(gene_name)]
    detect_dt.sel$group = "selected"
    detect_dt.sel[order(fraction)]
    
    # detect_dt = rbind(detect_dt, detect_dt.sel)
    
    rna_dt = rna_dt[gene_name %in% detect_dt.sel[fraction > min_detect_rate]$gene_name]
    
    
    
    
    rna_dt[, zscore := (expression - mean(expression)) / sd(expression), .(gene_name)]
    
    # ggplot(rna_dt[gene_name %in% sample(unique(gene_name), 5)], aes(x = gene_name, y = expression)) + geom_boxplot()
    # ggplot(rna_dt[gene_name %in% sample(unique(gene_name), 5)], aes(x = gene_name, y = zscore)) + geom_boxplot()
    # 
    if(value.var == "expression"){
        rna_dt.wide = dcast(rna_dt, gene_name~id, value.var = "expression")
    }else if(value.var == "zscore"){
        rna_dt.wide = dcast(rna_dt, gene_name~id, value.var = "zscore")
    }else{
        stop(value.var, " : value.var not OK")
    }
    
    rna_mat = as.matrix(rna_dt.wide[,-1])
    dim(rna_mat)
    rownames(rna_mat) = rna_dt.wide$gene_name
    
    pca_res = prcomp(rna_mat)
    
    pc = pca_res$x[, 1:n_dim]
    
    hc_res = hclust(dist(pc))
    gene_lev = rownames(pc)[hc_res$order]
    rna_dt$gene_name = factor(rna_dt$gene_name, levels = gene_lev)
    detect_dt$gene_name = factor(detect_dt$gene_name, levels = gene_lev)
    
    rna_dt
}


plot_heatmap_meta_clust = function(rna_dt, meta_dt, n_per = 100, cluster_colors = get_clusters_colors()){
    meta_dt.downsample = meta_dt[id %in% rna_dt$id]
    meta_dt.downsample = meta_dt.downsample[, .(id = ssvRecipes::sampleCap(id, n_per)), .(meta_cluster)]
    meta_dt.downsample = merge(meta_dt.downsample, meta_dt[, .(id, orig.ident)], by= "id")
    meta_dt.downsample$id = factor(meta_dt.downsample$id, levels = cell_o)
    
    rna_dt.downsample = merge(rna_dt, meta_dt.downsample[, .(id, meta_cluster)], by = "id")
    rna_dt.downsample$id = factor(rna_dt.downsample$id, levels =  cell_o)
    
    # meta_dt.downsample = meta_dt.downsample[id %in% rna_dt$id]
    # rna_dt.downsample = rna_dt.downsample[id %in% rna_dt$id]
    
    p_marker_heatmap = ggplot(rna_dt.downsample, aes(x = id, y = gene_name, fill= expression)) + geom_raster() + scale_fill_viridis_c() +
        theme(axis.ticks.x = element_blank(), 
              axis.text.x = element_blank(),
              axis.text.y = element_text(size = 6)) +
        labs(x = "", y = "")
    
    # ggsave("tmp.png", p_marker_heatmap, width = 6, height = 12)
    
    # ggplot(rna_dt.downsample[grepl("df", id)], aes(x = id, y = gene_name, fill= expression)) + geom_raster() + scale_fill_viridis_c() +
    #     theme(axis.ticks.x = element_blank(), axis.text.x = element_blank())
    # 
    # ggplot(rna_dt.downsample[grepl("wt", id)], aes(x = id, y = gene_name, fill= expression)) + geom_raster() + scale_fill_viridis_c() +
    #     theme(axis.ticks.x = element_blank(), axis.text.x = element_blank())
    
    
    p_clust = ggplot(meta_dt.downsample, aes(x = id, y = 1, fill = meta_cluster)) + geom_tile() +
        scale_fill_manual(values = cluster_colors) +
        theme_void() +
        coord_cartesian(ylim = c(0,2)) +
        theme(legend.key.size = grid::unit(.3, units = "cm")) +
        guides(fill= guide_legend(override.aes = list(size = .5), ncol = 6))
    p_genotype = ggplot(meta_dt.downsample, aes(x = id, y = 1, fill = orig.ident)) + geom_tile() + 
        scale_fill_manual(values = c("df4" = "red", "wt" = "black")) +
        theme_void()+
        coord_cartesian(ylim = c(0,2)) +
        theme(legend.key.size = grid::unit(.3, units = "cm")) +
        guides(fill= guide_legend(override.aes = list(size = .5), ncol = 3))
    
    g = ssvRecipes::sync_width(list(p_marker_heatmap, p_clust, p_genotype))
    pg = cowplot::plot_grid(plotlist = g, ncol = 1, rel_heights = c(16, 1, 1))
    parts = list(heatmap = p_marker_heatmap, clusters = p_clust, genotype = p_genotype)
    parts$assembled = plot_heatmap_meta_clust.assemble(parts)
    return(parts)
}

plot_heatmap_meta_clust.assemble = function(plots){
    g = ssvRecipes::sync_width(list(plots$heatmap, plots$clusters, plots$genotype))
    pg = cowplot::plot_grid(plotlist = g, ncol = 1, rel_heights = c(16, 1, 1))
    pg
}


#' simple counting at qgr
#'
#' @param bam_file 
#' @param qgr 
#' @param id_prefix 
#'
#' @return
#' @export
#'
#' @examples
sc_fetch_count = function(bam_file, qgr, id_prefix = NULL){
    what = c("rname", "strand", "pos", "qwidth", "cigar", "CR")
    what = Rsamtools::scanBamWhat()
    sbParam = Rsamtools::ScanBamParam(
        which = qgr,
        what = what, 
        tag = c("CR", "UB"))
    bam_raw = Rsamtools::scanBam(bam_file, param = sbParam)
    
    bam_dt = lapply(bam_raw, function(x){
        data.table(seqnames = x$rname, strand = x$strand,
                   start = x$pos, width = x$qwidth, cigar = x$cigar, id = x$tag$CR, umi = x$tag$UB)
    })
    bam_dt = data.table::rbindlist(bam_dt,
                                   use.names = TRUE,
                                   idcol = "which_label")
    
    bam_dt = bam_dt[!is.na(width)]
    toflip = sub(":-", "", as.character(subset(qgr, strand == "-")))
    target_strand = "+"
    if(target_strand %in% c("+", "-")){
        bam_dt = bam_dt[strand == target_strand & !(which_label %in% toflip) |
                            strand != target_strand & which_label %in% toflip]
    }
    if(!is.null(id_prefix)){
        bam_dt$id = paste0(id_prefix, bam_dt$id)
    }
    (bam_dt[, .(id, umi)] %>% unique)[, .(count = .N), .(id)]
}

sc_fetch_count.list = function(bam_list, qgr){
    stopifnot(!is.null(names(bam_list)))
    lapply(names(bam_list), function(nam){
        f = bam_list[[nam]]
        sc_fetch(f, qgr, id_prefix = nam)
    }) %>% rbindlist
}

#' pileup across qgr
#'
#' @param bam_file 
#' @param qgr 
#' @param id_prefix 
#'
#' @return
#' @export
#'
#' @examples
sc_fetch_pileup = function(bam_file, qgr, id_prefix = NULL){
    what = c("rname", "strand", "pos", "qwidth", "cigar", "CR")
    what = Rsamtools::scanBamWhat()
    sbParam = Rsamtools::ScanBamParam(
        which = qgr,
        what = what, 
        tag = c("CR", "UB"))
    bam_raw = Rsamtools::scanBam(bam_file, param = sbParam)
    
    bam_dt = lapply(bam_raw, function(x){
        data.table(seqnames = x$rname, strand = x$strand,
                   start = x$pos, width = x$qwidth, cigar = x$cigar, id = x$tag$CR, umi = x$tag$UB)
    })
    bam_dt = data.table::rbindlist(bam_dt,
                                   use.names = TRUE,
                                   idcol = "which_label")
    
    bam_dt = bam_dt[!is.na(width)]
    toflip = sub(":-", "", as.character(subset(qgr, strand == "-")))
    target_strand = "+"
    if(target_strand %in% c("+", "-")){
        bam_dt = bam_dt[strand == target_strand & !(which_label %in% toflip) |
                            strand != target_strand & which_label %in% toflip]
    }
    if(!is.null(id_prefix)){
        bam_dt$id = paste0(id_prefix, bam_dt$id)
    }
    (bam_dt[, .(id, umi)] %>% unique)[, .(count = .N), .(id)]
}


