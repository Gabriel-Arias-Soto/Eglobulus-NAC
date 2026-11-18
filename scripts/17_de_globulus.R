# 17_de_globulus.R -> expresion diferencial en HOJA contra el genoma
# conspecifico GCA_046226855.1 (E. globulus), a partir de la cuantificacion que
# deja 16_quant_hoja.sh.
#
# POR QUE SE REHACE. La tabla anterior (DE_NAC_modelo_conjunto.csv) esta en
# loci Eucgr de E. grandis, asi que para pintarla junto al censo habia que
# cruzarla por el NOMBRE del simbolo. Y el nombre no identifica: "EglNAC33"
# es ACJRO7_012424 en el censo y ACJRO7_023127 en el manuscrito. Ese cruce
# metio en T10 las cifras de otro gen. Cuantificando contra E. globulus la
# tabla sale ya en ACJRO7_*, que es el identificador del censo y el del
# transcrito que amplifican los cebadores (T11.transcrito -> t2g.tsv), de modo
# que no queda ningun emparejamiento por texto en toda la cadena.
#
# Se ejecuta desde la RAIZ:  Rscript scripts/17_de_globulus.R
if (!exists("ANCHO_COLUMNA")) source("scripts/00_estilo.R")
suppressPackageStartupMessages({
  library(tximport); library(edgeR); library(jsonlite); library(patchwork)
})

RNASEQ <- file.path(RAIZ, "rnaseq_udec")
QUANT  <- file.path(RNASEQ, "quant")
# el control se llama "NA" en el articulo original (non-acclimated) y readr lo
# leeria como ausente: hay que desactivar la conversion
dis <- suppressMessages(read_csv(file.path(TABLAS, "T33_diseno_rnaseq_hoja.csv"),
                                 na = character(), show_col_types = FALSE))
NIVELES <- c("NA", "CABF", "CAAF", "DS")     # control primero
dis$condicion <- factor(dis$condicion, levels = NIVELES)
stopifnot(!any(is.na(dis$condicion)), all(table(dis$condicion) == 3))

# LO QUE EL DEPOSITO PERMITE AFIRMAR. Las 12 corridas son 12 experimentos SRX
# distintos -12 preparaciones de libreria- y los autores las rotularon R1, R2 y
# R3, pero las tres de cada condicion cuelgan de UNA SOLA biomuestra (BioSample).
# Tal como esta depositado, la replicacion es de libreria, no de planta: no hay
# constancia de tres individuos por condicion. El diseno tampoco se puede
# corregir desde fuera, asi que se analiza como lo rotularon los autores y se
# deja que el dato responda: la BCV comun separa las dos situaciones sin
# ambiguedad (tecnica ~0,02-0,10; biologica en planta ~0,20-0,40). La cifra se
# registra abajo con su lectura y tiene que ir en Metodos.
bio_por_cond <- tapply(dis$biomuestra, dis$condicion, function(z) length(unique(z)))
stopifnot(!anyNA(bio_por_cond))
registrar_cifra("de_biomuestras_por_condicion",
                paste(sprintf("%s=%d", names(bio_por_cond), bio_por_cond), collapse = ", "),
                "T33.biomuestra (BioSample de ENA)")
registrar_cifra("de_estructura_replicas",
                sprintf("%d corridas, %d experimentos SRX, %d biomuestras",
                        nrow(dis), length(unique(dis$experimento)),
                        length(unique(dis$biomuestra))),
                "T33_diseno_rnaseq_hoja.csv")
if (all(bio_por_cond == 1))
  message("  AVISO DE: una sola biomuestra por condicion; las 3 replicas son de libreria.")

archivos <- file.path(QUANT, dis$corrida, "quant.sf")
names(archivos) <- dis$corrida
faltan <- dis$corrida[!file.exists(archivos)]
if (length(faltan)) stop("sin cuantificar: ", paste(faltan, collapse = ", "),
                         "\ncorre antes  bash scripts/16_quant_hoja.sh")

# --- control de calidad ------------------------------------------------------
# num_processed de salmon frente a read_count de ENA: es la comprobacion de que
# el fastq entro entero. Siete de las doce descargas previas estaban truncadas
# y salmon no siempre lo denuncia, asi que la barrera se queda en el script.
meta <- lapply(dis$corrida, function(s)
  fromJSON(file.path(QUANT, s, "aux_info", "meta_info.json")))
names(meta) <- dis$corrida
qc <- dis |>
  mutate(lecturas_salmon = vapply(meta, function(z) as.numeric(z$num_processed), 0),
         pct_mapeo       = vapply(meta, function(z) as.numeric(z$percent_mapped), 0),
         tipo_libreria   = vapply(meta, function(z) paste(z$library_types, collapse = ","), ""),
         integra         = lecturas_salmon == lecturas_ENA)
if (!all(qc$integra))
  stop("librerias truncadas (salmon != ENA): ",
       paste(qc$corrida[!qc$integra], collapse = ", "))
write_csv(qc, file.path(TABLAS, "T34_qc_rnaseq_hoja.csv"))

# --- lectura y agregacion a gen ---------------------------------------------
t2g <- suppressMessages(read_tsv(file.path(RNASEQ, "t2g.tsv"),
                                 col_names = c("transcrito", "locus"), show_col_types = FALSE))
txi <- tximport(archivos, type = "salmon", tx2gene = t2g)

# Offsets de longitud efectiva segun la vineta de tximport para edgeR: corrigen
# que el mismo gen tenga longitud efectiva distinta en cada libreria (las
# profundidades van de 12 a 108 M de lecturas).
cts <- txi$counts
nm  <- txi$length / exp(rowMeans(log(txi$length)))
lib <- calcNormFactors(cts / nm) * colSums(cts / nm)
y   <- DGEList(cts, group = dis$condicion, samples = dis)
y   <- scaleOffset(y, log(sweep(nm, 2, lib, "*")))

diseno <- model.matrix(~ 0 + condicion, data = dis)
colnames(diseno) <- NIVELES
mantener <- filterByExpr(y, design = diseno)
y <- y[mantener, , keep.lib.sizes = FALSE]

y   <- estimateDisp(y, diseno, robust = TRUE)
fit <- glmQLFit(y, diseno, robust = TRUE)

# Los nombres son los mismos que usa 07_expresion.R en CONTRASTES, para que la
# figura no tenga que cambiar de constante. Se construye la matriz a mano: con
# makeContrasts el nivel "NA" rompe el parseo de la formula.
CONTR <- list(Leaf_NA_DS     = c("DS",   "NA"),
              Leaf_NA_CABF   = c("CABF", "NA"),
              Leaf_NA_CAAF   = c("CAAF", "NA"),
              Leaf_CABF_CAAF = c("CAAF", "CABF"))
cm <- matrix(0, nrow = length(NIVELES), ncol = length(CONTR),
             dimnames = list(NIVELES, names(CONTR)))
for (k in names(CONTR)) { cm[CONTR[[k]][1], k] <- 1; cm[CONTR[[k]][2], k] <- -1 }

de <- lapply(names(CONTR), function(k) {
  tt <- topTags(glmQLFTest(fit, contrast = cm[, k]), n = Inf, sort.by = "none")$table
  tibble(locus = rownames(tt), contraste = k, logFC = tt$logFC, logCPM = tt$logCPM,
         PValue = tt$PValue, FDR = tt$FDR)
}) |> bind_rows()
write_csv(de, file.path(DATOS, "DE_todos_globulus.csv"))

# --- atribuibilidad de la senal ---------------------------------------------
# Con la referencia de E. grandis, "senal atribuible" queria decir que el locus
# era el unico del censo asignado a su co-ortologo. Contra el conspecifico esa
# pregunta desaparece y la sustituye una que se mide: que fraccion de las
# lecturas del locus salmon pudo asignar de forma unica, en vez de repartirla
# por el EM entre transcritos indistinguibles (tipico entre duplicados en
# tandem recientes). ambig_info.tsv va en el MISMO orden que quant.sf.
tx_ord <- suppressMessages(read_tsv(file.path(QUANT, dis$corrida[1], "quant.sf"),
                                    show_col_types = FALSE))$Name
g <- factor(t2g$locus[match(tx_ord, t2g$transcrito)])
stopifnot(!anyNA(g))
frac <- vapply(dis$corrida, function(s) {
  a <- suppressMessages(read_tsv(file.path(QUANT, s, "aux_info", "ambig_info.tsv"),
                                 show_col_types = FALSE))
  stopifnot(nrow(a) == length(g))
  u <- rowsum(a$UniqueCount, g)[, 1]; m <- rowsum(a$AmbigCount, g)[, 1]
  ifelse(u + m > 0, u / (u + m), NA_real_)
}, numeric(nlevels(g)))
atrib <- tibble(locus = levels(g),
                frac_unica = apply(frac, 1, function(z)
                  if (all(is.na(z))) NA_real_ else median(z, na.rm = TRUE)))
write_csv(atrib, file.path(TABLAS, "T35_atribuibilidad_globulus.csv"))

# --- subconjunto NAC ---------------------------------------------------------
censo <- leer_censo()
nac <- censo |>
  select(locus, simbolo, subfamilia_final, es_focal, simbolo_manuscrito_2026) |>
  left_join(atrib, by = "locus") |>
  inner_join(de, by = "locus") |>
  arrange(locus, contraste)
write_csv(nac, file.path(TABLAS, "T36_DE_NAC_globulus.csv"))

sin_de <- setdiff(censo$locus, unique(nac$locus))
if (length(sin_de))
  message("  AVISO: ", length(sin_de), " loci del censo sin estimacion (no pasan filterByExpr)")

# --- figura de control -------------------------------------------------------
COL_COND <- c("NA" = "grey40", "CABF" = "#0072B2", "CAAF" = "#56B4E9", "DS" = "#D55E00")
pa <- ggplot(qc, aes(lecturas_salmon / 1e6, pct_mapeo, colour = condicion)) +
  geom_point(size = 1.8) +
  ggrepel::geom_text_repel(aes(label = sub("^SRR", "", corrida)), size = 2.1,
                           seed = 1, min.segment.length = 0.2, max.overlaps = Inf) +
  scale_colour_manual(values = COL_COND, name = NULL) +
  labs(x = "reads (millions)", y = "mapping (%)") +
  theme_nac(8) + theme(legend.position = "none")

mds <- limma::plotMDS(y, plot = FALSE)
pb <- ggplot(tibble(x = mds$x, y = mds$y, condicion = dis$condicion,
                    rep = dis$replica),
             aes(x, y, colour = condicion)) +
  geom_point(size = 2) +
  ggrepel::geom_text_repel(aes(label = rep), size = 2.1, seed = 1, max.overlaps = Inf) +
  scale_colour_manual(values = COL_COND, name = NULL) +
  labs(x = sprintf("dim 1 (%.0f %%)", 100 * mds$var.explained[1]),
       y = sprintf("dim 2 (%.0f %%)", 100 * mds$var.explained[2])) +
  theme_nac(8) + theme(legend.position = "bottom", legend.margin = margin(0, 0, 0, 0))

pc <- ggplot(tibble(cpm = y$AveLogCPM, bcv = sqrt(y$tagwise.dispersion),
                    nac = rownames(y) %in% censo$locus),
             aes(cpm, bcv)) +
  geom_point(data = ~ subset(.x, !nac), colour = "grey80", size = 0.25) +
  geom_point(data = ~ subset(.x, nac), colour = COL_DESTACADO, size = 0.7) +
  geom_hline(yintercept = sqrt(y$common.dispersion), linetype = 2, linewidth = 0.3) +
  labs(x = expression(mean~log[2]~CPM), y = "BCV") +
  theme_nac(8)

guardar_figura((pa | pb | pc) + plot_annotation(tag_levels = "a"),
               "F12_qc_rnaseq_globulus", alto = 2.5)

# --- cifras ------------------------------------------------------------------
sig <- nac |> filter(FDR < 0.05, abs(logFC) >= 1)
registrar_cifra("de_referencia", "GCA_046226855.1 (E. globulus), 42.218 transcritos",
                "16_quant_hoja.sh")
registrar_cifra("de_estudio",
                sprintf("PRJNA486291, hoja, 4 condiciones x 3 %s",
                        if (all(bio_por_cond == 1)) "replicas de libreria" else "replicas"),
                "T33_diseno_rnaseq_hoja.csv")
registrar_cifra("de_lecturas_rango",
                sprintf("%.1f a %.1f M", min(qc$lecturas_salmon) / 1e6, max(qc$lecturas_salmon) / 1e6),
                "T34.lecturas_salmon")
registrar_cifra("de_mapeo_rango", sprintf("%.1f a %.1f %%", min(qc$pct_mapeo), max(qc$pct_mapeo)),
                "T34.pct_mapeo")
registrar_cifra("de_genes_tras_filtro", nrow(y), "filterByExpr sobre el diseno completo")
registrar_cifra("de_nac_con_estimacion", length(unique(nac$locus)),
                sprintf("de %d loci del censo", nrow(censo)))
bcv <- sqrt(y$common.dispersion)
registrar_cifra("de_bcv_comun", sprintf("%.3f", bcv), "edgeR estimateDisp")
# La lectura de la BCV es la que decide como se puede citar la significacion.
# Si sale en rango tecnico, los p y los FDR miden reproducibilidad de libreria,
# no variacion entre plantas: el log2FC sigue siendo la magnitud estimada, pero
# la inferencia sobre la especie la sostiene el qRT-PCR, que si tiene replicas
# biologicas. Esta frase va tal cual a Metodos.
registrar_cifra("de_bcv_lectura",
                if (bcv < 0.10)
                  sprintf(paste("BCV comun %.3f, en rango de replicacion tecnica: los FDR miden",
                                "reproducibilidad entre librerias de una misma biomuestra y no",
                                "variacion biologica; el log2FC es descriptivo"), bcv)
                else if (bcv < 0.20)
                  sprintf(paste("BCV comun %.3f, intermedia: por debajo de lo habitual entre",
                                "plantas, compatible con una sola biomuestra por condicion"), bcv)
                else
                  sprintf(paste("BCV comun %.3f, en rango de variacion biologica en planta pese",
                                "a la biomuestra unica"), bcv),
                "edgeR estimateDisp + T33.biomuestra")
registrar_cifra("de_nac_significativos", length(unique(sig$locus)),
                "FDR < 0,05 y |log2FC| >= 1 en al menos un contraste de hoja")
registrar_cifra("de_atribuibilidad_focal",
                sprintf("%.3f", atrib$frac_unica[atrib$locus == GEN_FOCAL_LOCUS]),
                "T35: fraccion de lecturas asignadas de forma unica")
volcar_cifras("de_globulus")
message("  T36_DE_NAC_globulus.csv: ", nrow(nac), " filas")
