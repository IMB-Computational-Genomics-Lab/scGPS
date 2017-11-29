clusters <- readRDS("/Users/quan.nguyen/Downloads/APC_clusters.RDS")
DEgenes <-read.table("/Users/quan.nguyen/Downloads/APC_DEG_C4vsOthers.txt")
Exprs <- readRDS("/Users/quan.nguyen/Downloads/APC_exp.RDS")
allDat <-NewscGPS(ExpressionMatrix = Exprs, GeneMetadata = row.names(Exprs),
                  CellMetadata = clusters)

genes <-as.vector(DEgenes$V1)

test <- bootstrap_scGPS(nboots = 5,mixedpop1 = allDat, mixedpop2 = allDat, genes =genes, c_selectID = 1, listData =list() )

