# Entrega 3, paso 3.2: enriquecimiento GO de los NAC de E. globulus (como en el artículo 2:
# clusterProfiler, análisis de sobrerrepresentación).
#
# Entradas
#   tablas/T41_homologos_Ath_NAC189.csv          EglNAC -> mejor homólogo de Arabidopsis (paso 3.1)
#   datos/entrega3/go/Egl_mejor_homologo_Ath.tsv  mejor homólogo de cada gen de E. globulus (paso 3.1)
#
# Análisis principal (unidad = gen de Arabidopsis)
#   lista    = AGI distintos de los NAC con homólogo
#   universo = AGI distintos de todo el proteoma con homólogo
#   enrichGO con org.At.tair.db, ontologías BP, MF y CC; BH; p.adjust <= 0.05 y q <= 0.2;
#   tamaño de conjunto 10-500 (valores por defecto); simplify (Wang, corte 0.7) en BP.
#
# Contraste 1 (unidad = gen de E. globulus)
#   Cada gen de E. globulus hereda los términos GO (con ancestros) de su homólogo.
#   lista = genes NAC con homólogo; universo = genes de E. globulus con homólogo.
#   Se prueban exactamente los mismos términos que en el análisis principal, así que la
#   corrección BH cubre el mismo número de pruebas y solo cambia la unidad de conteo.
#
# Salidas en datos/entrega3/go/
#   enrichGO_<ONT>_todos.tsv, enrichGO_<ONT>_signif.tsv, enrichGO_BP_simplify.tsv
#   enricherEgl_<ONT>_todos.tsv, enricherEgl_<ONT>_signif.tsv
#   anotacion_GOALL_universo.tsv.gz  pares AGI-término (con ancestros) para el recálculo en Python
#   go_objetos.rds                 objetos enrichResult para las figuras
#   resumen_go.json

suppressPackageStartupMessages({
  library(clusterProfiler); library(org.At.tair.db); library(AnnotationDbi); library(GO.db)
  library(enrichplot); library(jsonlite)
})

A   <- paste0(Sys.getenv("NAC_ROOT", "."), "/ANALISIS_Eglobulus")
GOD <- file.path(A, "datos/entrega3/go")
t41 <- read.csv(file.path(A, "tablas/T41_homologos_Ath_NAC189.csv"), stringsAsFactors = FALSE)
mej <- read.delim(file.path(GOD, "Egl_mejor_homologo_Ath.tsv"), stringsAsFactors = FALSE)

ok    <- t41$pasa_filtros %in% c(TRUE, "True")
genes <- unique(t41$AGI[ok])
univ  <- unique(mej$AGI)
stopifnot(all(genes %in% univ), all(t41$locus[ok] %in% mej$gen_Egl))
cat("lista AGI:", length(genes), "| universo AGI:", length(univ), "| genes Egl NAC:", sum(ok),
    "| universo Egl:", nrow(mej), "\n")

P_CORTE <- 0.05; Q_CORTE <- 0.2; MIN_GS <- 10; MAX_GS <- 500
ONTS <- c("BP", "MF", "CC")
set.seed(20261013)

# anotación con ancestros, construida igual que clusterProfiler:::get_GO_data
go_ont <- AnnotationDbi::Ontology(GO.db::GOTERM)
ann <- suppressMessages(AnnotationDbi::select(org.At.tair.db, keys = univ, keytype = "TAIR", columns = "GOALL"))
ann <- unique(ann[!is.na(ann$GOALL), c("TAIR", "GOALL")])
ann$ONTOLOGIA <- unname(go_ont[ann$GOALL])
ann <- ann[!is.na(ann$ONTOLOGIA), ]
gz <- gzfile(file.path(GOD, "anotacion_GOALL_universo.tsv.gz"), "w")
write.table(ann, gz, sep = "\t", quote = FALSE, row.names = FALSE); close(gz)
nombres_go <- AnnotationDbi::Term(GO.db::GOTERM)

objetos <- list(); resumen <- list()
for (ont in ONTS) {
  ego <- enrichGO(gene = genes, universe = univ, OrgDb = org.At.tair.db, keyType = "TAIR", ont = ont,
                  pAdjustMethod = "BH", pvalueCutoff = P_CORTE, qvalueCutoff = Q_CORTE,
                  minGSSize = MIN_GS, maxGSSize = MAX_GS, readable = FALSE)
  todo <- enrichGO(gene = genes, universe = univ, OrgDb = org.At.tair.db, keyType = "TAIR", ont = ont,
                   pAdjustMethod = "BH", pvalueCutoff = 1, qvalueCutoff = 1,
                   minGSSize = MIN_GS, maxGSSize = MAX_GS, readable = FALSE)
  ego_r <- setReadable(ego, OrgDb = org.At.tair.db, keyType = "TAIR")
  sig <- as.data.frame(ego_r)
  write.table(todo@result, file.path(GOD, sprintf("enrichGO_%s_todos.tsv", ont)), sep = "\t", quote = FALSE, row.names = FALSE)
  write.table(sig, file.path(GOD, sprintf("enrichGO_%s_signif.tsv", ont)), sep = "\t", quote = FALSE, row.names = FALSE)
  objetos[[paste0("ego_", ont)]] <- ego_r

  # contraste 1: genes de E. globulus como unidades. La anotación es la de toda la ontología
  # (así el universo anotado se define igual que en enrichGO); después se conservan solo los
  # términos probados en el análisis principal y se rehace BH sobre ellos.
  probados <- todo@result$ID
  a_ont <- ann[ann$ONTOLOGIA == ont, c("GOALL", "TAIR")]
  t2g <- merge(a_ont, mej[, c("gen_Egl", "AGI")], by.x = "TAIR", by.y = "AGI")[, c("GOALL", "gen_Egl")]
  t2n <- data.frame(GOALL = unique(t2g$GOALL), nombre = unname(nombres_go[unique(t2g$GOALL)]))
  egl_lista <- t41$locus[ok]
  egl_univ  <- unique(mej$gen_Egl)
  enr <- enricher(gene = egl_lista, universe = egl_univ, TERM2GENE = t2g, TERM2NAME = t2n,
                  pAdjustMethod = "BH", pvalueCutoff = 1, qvalueCutoff = 1, minGSSize = 1, maxGSSize = Inf)
  enr_todo <- enr@result[enr@result$ID %in% probados, ]
  stopifnot(setequal(enr_todo$ID, probados))
  enr_todo$p.adjust <- p.adjust(enr_todo$pvalue, method = "BH")
  enr_todo$qvalue <- tryCatch(qvalue::qvalue(enr_todo$pvalue, lambda = 0.05, pi0.method = "bootstrap")$qvalues,
                              error = function(e) NA_real_)
  enr_todo <- enr_todo[order(enr_todo$pvalue), ]
  enr@result <- enr_todo
  enr_sig <- enr_todo[enr_todo$pvalue <= P_CORTE & enr_todo$p.adjust <= P_CORTE &
                        (is.na(enr_todo$qvalue) | enr_todo$qvalue <= Q_CORTE), ]
  write.table(enr_todo, file.path(GOD, sprintf("enricherEgl_%s_todos.tsv", ont)), sep = "\t", quote = FALSE, row.names = FALSE)
  write.table(enr_sig, file.path(GOD, sprintf("enricherEgl_%s_signif.tsv", ont)), sep = "\t", quote = FALSE, row.names = FALSE)
  objetos[[paste0("enrEgl_", ont)]] <- enr

  resumen[[ont]] <- list(terminos_probados = nrow(todo@result), significativos = nrow(sig),
                         probados_contraste_Egl = nrow(enr_todo), significativos_contraste_Egl = nrow(enr_sig),
                         compartidos = length(intersect(sig$ID, enr_sig$ID)),
                         N_universo_anotado = length(todo@universe), n_lista_anotada = sum(genes %in% todo@universe))
  cat(ont, ": probados", nrow(todo@result), "| signif", nrow(sig), "| contraste Egl: probados", nrow(enr_todo),
      "signif", nrow(enr_sig), "| compartidos", resumen[[ont]]$compartidos, "\n")
}

ego_bp_s <- simplify(objetos$ego_BP, cutoff = 0.7, by = "p.adjust", select_fun = min, measure = "Wang")
write.table(as.data.frame(ego_bp_s), file.path(GOD, "enrichGO_BP_simplify.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
objetos$ego_BP_simplify <- ego_bp_s
resumen$BP$tras_simplify <- nrow(as.data.frame(ego_bp_s))
cat("BP tras simplify:", resumen$BP$tras_simplify, "\n")

saveRDS(objetos, file.path(GOD, "go_objetos.rds"))
info_go  <- GO_dbInfo();           info_org <- org.At.tair_dbInfo()
resumen$fuentes <- list(
  GO.db = as.character(packageVersion("GO.db")), org.At.tair.db = as.character(packageVersion("org.At.tair.db")),
  clusterProfiler = as.character(packageVersion("clusterProfiler")),
  GO_fecha = info_go$value[info_go$name == "GOSOURCEDATE"],
  TAIR_GO_fecha = paste(info_org$value[grepl("SOURCEDATE", info_org$name)], collapse = "; "))
resumen$parametros <- list(pAdjustMethod = "BH", pvalueCutoff = P_CORTE, qvalueCutoff = Q_CORTE,
                           minGSSize = MIN_GS, maxGSSize = MAX_GS, simplify = "Wang, 0.7, p.adjust mínimo",
                           lista_AGI = length(genes), universo_AGI = length(univ),
                           lista_Egl = sum(ok), universo_Egl = nrow(mej))
write_json(resumen, file.path(GOD, "resumen_go.json"), auto_unbox = TRUE, pretty = TRUE)
print(resumen$fuentes)
sessionInfo()
