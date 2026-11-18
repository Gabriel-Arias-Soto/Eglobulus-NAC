# Entrega 3, paso 3.3 (datos de figuras y tablas de GO). Se corre en nac3, después de
# e3_05_go_enriquecimiento.R. Mismos parámetros que el análisis principal: universo = AGI con
# homólogo, BH, p y p.adjust <= 0.05, q <= 0.2, tamaño de conjunto 10-500.
#
#   1. compareCluster por subfamilia (BP). Unidad = AGI distinto. Una subfamilia entra si sus
#      genes con homólogo dan al menos MIN_AGI AGI distintos: con uno o dos genes, cualquier
#      término de 10 genes al que estén anotados sale "enriquecido" y el resultado solo
#      repite la anotación de ese gen. Se repite por grupo (I, II, ... XI) como contraste.
#      Para la figura se quitan además los términos con un solo gen (Count < MIN_COUNT) y la
#      redundancia (simplify, Wang 0,7, se queda el de p.adjust mínimo, igual que en el
#      análisis principal): cc_<nivel>_BP_simplify.tsv.
#   2. Sensibilidad sin límite superior de tamaño (maxGSSize = Inf) en BP, MF y CC: los
#      términos propios de un factor de transcripción (actividad de unión a ADN, núcleo)
#      tienen más de 500 genes y el análisis principal no los prueba. Sirven de control
#      positivo de la transferencia de anotación.
#   3. Similitud semántica de Wang entre los términos BP significativos y su agrupamiento
#      jerárquico, igual que enrichplot::treeplot (1 - similitud, hclust "ward.D", 5 grupos).
#   4. Anotación GO directa (sin ancestros), con código de evidencia, de los AGI de la lista.
#
# Salidas en datos/entrega3/go/
#   cc_subfamilia_BP_todos.tsv, cc_grupo_BP_todos.tsv     todos los términos probados
#   cc_subfamilia_BP_simplify.tsv, cc_grupo_BP_simplify.tsv  significativos con Count >= 2, sin redundancia
#   cc_listas.tsv                                        AGI de cada subfamilia y grupo
#   sens_maxInf_<ONT>_todos.tsv
#   termsim_BP_signif.tsv, arbol_BP_signif.nwk, arbol_BP_grupos.tsv
#   anotacion_directa_lista.tsv
#   resumen_go_figuras.json

suppressPackageStartupMessages({
  library(clusterProfiler); library(org.At.tair.db); library(AnnotationDbi); library(GO.db)
  library(GOSemSim); library(ape); library(jsonlite)
})

A   <- paste0(Sys.getenv("NAC_ROOT", "."), "/ANALISIS_Eglobulus")
GOD <- file.path(A, "datos/entrega3/go")
t41 <- read.csv(file.path(A, "tablas/T41_homologos_Ath_NAC189.csv"), stringsAsFactors = FALSE)
t01 <- read.csv(file.path(A, "tablas/T01_censo_maestro_189.csv"), stringsAsFactors = FALSE)
mej <- read.delim(file.path(GOD, "Egl_mejor_homologo_Ath.tsv"), stringsAsFactors = FALSE)
stopifnot(nrow(t41) == 189, setequal(t41$simbolo, t01$simbolo))

P_CORTE <- 0.05; Q_CORTE <- 0.2; MIN_GS <- 10; MAX_GS <- 500; MIN_AGI <- 3; N_GRUPOS_ARBOL <- 5
MIN_COUNT <- 2
PARCIAL <- "no evaluada (dominio parcial)"

ok   <- t41$pasa_filtros %in% c(TRUE, "True")
univ <- unique(mej$AGI)
z    <- merge(t41[ok, c("simbolo", "locus", "AGI")], t01[, c("simbolo", "subfamilia_final")], by = "simbolo")
stopifnot(nrow(z) == sum(ok))
z$grupo <- ifelse(z$subfamilia_final == PARCIAL, PARCIAL, sub("^([IVX]+).*$", "\\1", z$subfamilia_final))

signif_ <- function(d) d$pvalue <= P_CORTE & d$p.adjust <= P_CORTE & (is.na(d$qvalue) | d$qvalue <= Q_CORTE)

listas_por <- function(columna) {
  l <- lapply(split(z$AGI, z[[columna]]), unique)
  l <- l[names(l) != PARCIAL]
  list(entran = l[lengths(l) >= MIN_AGI], fuera = l[lengths(l) < MIN_AGI])
}

sd_bp <- suppressMessages(godata(annoDb = "org.At.tair.db", keytype = "TAIR", ont = "BP", computeIC = FALSE))

resumen <- list(parametros = list(MIN_AGI = MIN_AGI, MIN_COUNT_figura = MIN_COUNT,
                                  simplify = "Wang, 0.7, p.adjust mínimo",
                                  pvalueCutoff = P_CORTE, qvalueCutoff = Q_CORTE,
                                  minGSSize = MIN_GS, maxGSSize = MAX_GS, universo_AGI = length(univ),
                                  arbol = sprintf("Wang; hclust ward.D sobre 1 - similitud; %d grupos", N_GRUPOS_ARBOL)))
filas_listas <- list()
for (columna in c("subfamilia_final", "grupo")) {
  et <- if (columna == "grupo") "grupo" else "subfamilia"
  L <- listas_por(columna)
  for (nm in names(L$entran)) filas_listas[[length(filas_listas) + 1]] <-
    data.frame(nivel = et, grupo = nm, entra = TRUE, n_AGI = length(L$entran[[nm]]),
               n_Egl = sum(z[[columna]] == nm), AGI = paste(L$entran[[nm]], collapse = ";"))
  for (nm in names(L$fuera)) filas_listas[[length(filas_listas) + 1]] <-
    data.frame(nivel = et, grupo = nm, entra = FALSE, n_AGI = length(L$fuera[[nm]]),
               n_Egl = sum(z[[columna]] == nm), AGI = paste(L$fuera[[nm]], collapse = ";"))
  cc <- compareCluster(geneClusters = L$entran, fun = "enrichGO", universe = univ, OrgDb = org.At.tair.db,
                       keyType = "TAIR", ont = "BP", pAdjustMethod = "BH", pvalueCutoff = 1, qvalueCutoff = 1,
                       minGSSize = MIN_GS, maxGSSize = MAX_GS, readable = FALSE)
  d <- as.data.frame(cc)
  d$significativo <- signif_(d)
  write.table(d, file.path(GOD, sprintf("cc_%s_BP_todos.tsv", et)), sep = "\t", quote = FALSE, row.names = FALSE)
  por <- tapply(d$significativo, d$Cluster, sum)
  cc_s <- cc
  cc_s@compareClusterResult <- d[d$significativo & d$Count >= MIN_COUNT, setdiff(names(d), "significativo")]
  cc_s <- simplify(cc_s, cutoff = 0.7, by = "p.adjust", select_fun = min, measure = "Wang", semData = sd_bp)
  ds <- as.data.frame(cc_s)
  write.table(ds, file.path(GOD, sprintf("cc_%s_BP_simplify.tsv", et)), sep = "\t", quote = FALSE, row.names = FALSE)
  resumen[[paste0("compareCluster_", et)]] <- list(
    entran = as.list(lengths(L$entran)), fuera = as.list(lengths(L$fuera)),
    probados = as.list(table(d$Cluster)), significativos = as.list(por),
    significativos_count_min = as.list(tapply(d$significativo & d$Count >= MIN_COUNT, d$Cluster, sum)),
    tras_simplify = as.list(table(factor(ds$Cluster, levels = levels(factor(d$Cluster))))),
    terminos_signif_distintos = length(unique(d$ID[d$significativo])))
  cat(et, "| entran:", paste(sprintf("%s(%d)", names(L$entran), lengths(L$entran)), collapse = " "),
      "\n   signif por grupo:", paste(sprintf("%s=%d", names(por), por), collapse = " "),
      "| términos distintos:", length(unique(d$ID[d$significativo])), "\n")
}
write.table(do.call(rbind, filas_listas), file.path(GOD, "cc_listas.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

# 2. sensibilidad sin límite superior de tamaño
genes <- unique(t41$AGI[ok])
for (ont in c("BP", "MF", "CC")) {
  s <- enrichGO(gene = genes, universe = univ, OrgDb = org.At.tair.db, keyType = "TAIR", ont = ont,
                pAdjustMethod = "BH", pvalueCutoff = 1, qvalueCutoff = 1, minGSSize = MIN_GS,
                maxGSSize = Inf, readable = FALSE)@result
  s$significativo <- signif_(s)
  s$supera_maxGSSize <- as.integer(sub("/.*", "", s$BgRatio)) > MAX_GS
  write.table(s, file.path(GOD, sprintf("sens_maxInf_%s_todos.tsv", ont)), sep = "\t", quote = FALSE, row.names = FALSE)
  top <- head(s[s$supera_maxGSSize & s$significativo, ], 3)
  resumen[[paste0("sensibilidad_", ont)]] <- list(
    probados = nrow(s), significativos = sum(s$significativo),
    significativos_mayores_500 = sum(s$significativo & s$supera_maxGSSize),
    primeros_mayores_500 = lapply(seq_len(nrow(top)), function(i)
      list(ID = top$ID[i], termino = top$Description[i], GeneRatio = top$GeneRatio[i],
           BgRatio = top$BgRatio[i], p.adjust = signif(top$p.adjust[i], 3))))
  cat("sensibilidad", ont, ": probados", nrow(s), "| signif", sum(s$significativo),
      "| signif con más de 500 genes", sum(s$significativo & s$supera_maxGSSize), "\n")
}

# 3. similitud de Wang y árbol de los términos BP significativos
sig <- read.delim(file.path(GOD, "enrichGO_BP_signif.tsv"), stringsAsFactors = FALSE)
S <- termSim(sig$ID, sig$ID, semData = sd_bp, method = "Wang")
stopifnot(isSymmetric(unname(S)), all(abs(diag(S) - 1) < 1e-12))
write.table(data.frame(ID = rownames(S), S, check.names = FALSE), file.path(GOD, "termsim_BP_signif.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
hc <- hclust(as.dist(1 - S), method = "ward.D")
g0 <- cutree(hc, k = N_GRUPOS_ARBOL)
# cutree numera por orden de aparición en los datos; se renumera por el orden de las hojas
# para que el grupo 1 sea el primero del dendrograma y la figura lea 1, 2, 3... seguidos
grupo_arbol <- setNames(match(g0, unique(g0[hc$order])), names(g0))
stopifnot(identical(as.vector(table(g0)[unique(g0[hc$order])]), as.vector(table(grupo_arbol))))
write.tree(as.phylo(hc), file.path(GOD, "arbol_BP_signif.nwk"))
ga <- data.frame(ID = names(grupo_arbol), grupo_arbol = unname(grupo_arbol), orden_hojas = match(names(grupo_arbol), hc$labels[hc$order]))
ga <- merge(ga, sig[, c("ID", "Description", "p.adjust", "Count")], by = "ID")
ga <- ga[order(ga$orden_hojas), ]
write.table(ga, file.path(GOD, "arbol_BP_grupos.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
resumen$arbol <- list(terminos = nrow(S), grupos = as.list(table(grupo_arbol)),
                      similitud_mediana_fuera_diagonal = round(median(S[upper.tri(S)]), 3))

# 4. anotación GO directa de los AGI de la lista
ad <- suppressMessages(AnnotationDbi::select(org.At.tair.db, keys = genes, keytype = "TAIR",
                                             columns = c("GO", "EVIDENCE", "ONTOLOGY")))
ad <- unique(ad[!is.na(ad$GO), ])
ad$termino <- unname(AnnotationDbi::Term(GO.db::GOTERM)[ad$GO])
write.table(ad, file.path(GOD, "anotacion_directa_lista.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
resumen$anotacion_directa <- list(AGI = length(unique(ad$TAIR)), pares_AGI_GO = nrow(unique(ad[, c("TAIR", "GO")])),
                                  por_ontologia = as.list(table(unique(ad[, c("TAIR", "GO", "ONTOLOGY")])$ONTOLOGY)))
resumen$versiones <- list(clusterProfiler = as.character(packageVersion("clusterProfiler")),
                          GOSemSim = as.character(packageVersion("GOSemSim")),
                          org.At.tair.db = as.character(packageVersion("org.At.tair.db")),
                          GO.db = as.character(packageVersion("GO.db")), ape = as.character(packageVersion("ape")))
write_json(resumen, file.path(GOD, "resumen_go_figuras.json"), auto_unbox = TRUE, pretty = TRUE)
cat("árbol: grupos", paste(table(grupo_arbol), collapse = " "), "| anotación directa:", nrow(ad), "filas\n")
sessionInfo()
