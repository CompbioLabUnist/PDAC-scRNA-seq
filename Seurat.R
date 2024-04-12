

seurat.obj <- NormalizeData(object = seurat.obj, normalization.method = "LogNormalize", scale.factor = 10000)
seurat.obj <- FindVariableFeatures(object = seurat.obj, selection.method = "vst")
seurat.obj <- CellCycleScoring(seurat.obj, s.features = cc.genes$s.genes, g2m.features = cc.genes$g2m.genes, set.ident = T)
seurat.obj <- ScaleData(object = seurat.obj, features = rownames(seurat.obj), vars.to.regress = c("nCount_RNA", "percent.mt", "S.Score", "G2M.Score"))

seurat.obj <- RunPCA(object = seurat.obj, features = VariableFeatures(object = seurat.obj), npcs = 17)
seurat.obj <- JackStraw(object = seurat.obj, num.replicate = 100, dims = 17)
seurat.obj <- ScoreJackStraw(object = seurat.obj, dims = 1:17)
seurat.obj <- FindNeighbors(object = seurat.obj, dims = 1:17, force.recalc = T)
seurat.obj <- FindClusters(object = seurat.obj, resolution = 0.6)
seurat.obj <- RunUMAP(object = seurat.obj, dims = 1:17)
seurat.obj <- RunTSNE(object = seurat.obj, dims = 1:17)

seurat.degs <- FindAllMarkers(object = seurat.obj, min.pct = 0.25, logfc.threshold = 0.25)
