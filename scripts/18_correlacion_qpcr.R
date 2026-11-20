# 18_correlacion_qpcr.R -> reconstruye T10 (qPCR frente a RNA-seq) sin ningun
# emparejamiento por nombre de gen.
#
# La T10 anterior se armo a mano y su columna rnaseq_log2FC salio de buscar el
# simbolo del manuscrito dentro de una tabla de DE en loci de E. grandis. Como
# el mismo simbolo designa loci distintos en cada nomenclatura, tres de las
# doce celdas traian cifras de otro gen y una cuarta traia un 0,000 que en
# realidad significaba "no evaluado". Aqui la cadena es:
#
#   T08.gen (simbolo del manuscrito)
#     -> T11.simbolo_manuscrito -> T11.locus_tag   (puente curado, verificado
#        cebador a cebador contra los 42.218 transcritos de la referencia)
#     -> T36.locus                                  (DE contra E. globulus)
#
# Ningun paso compara cadenas de simbolos entre nomenclaturas distintas.
#
# El contraste C_DS del qPCR (control regado -> deficit hidrico) se enfrenta a
# Leaf_NA_DS del RNA-seq. No son el mismo control: el RNA-seq no incluye
# librerias de planta regada y usa como referencia el NA del ensayo de frio,
# que ademas viene de otra remesa de secuenciacion (SRR9887586 y siguientes
# frente a SRR77264xx). Es la equivalencia que asumia el manuscrito y se
# conserva, pero hay que declararla en Metodos.
#
# Se ejecuta desde la RAIZ:  Rscript scripts/18_correlacion_qpcr.R
if (!exists("ANCHO_COLUMNA")) source("scripts/00_estilo.R")

t08 <- suppressMessages(read_csv(file.path(TABLAS, "T08_qpcr_expresion_relativa.csv"), show_col_types = FALSE))
t11 <- suppressMessages(read_csv(file.path(TABLAS, "T11_especificidad_primers_qpcr.csv"), show_col_types = FALSE))
t36 <- suppressMessages(read_csv(file.path(TABLAS, "T36_DE_NAC_globulus.csv"), show_col_types = FALSE))

EQUIV <- c(C_DS = "Leaf_NA_DS", NA_CABF = "Leaf_NA_CABF", NA_CAAF = "Leaf_NA_CAAF")

# --- lado qPCR: media de log2(RQ) por gen y contraste ------------------------
qpcr <- t08 |>
  filter(!is.na(contraste), contraste != "", RQ > 0) |>
  group_by(simbolo_manuscrito = gen, contraste) |>
  summarise(qpcr_log2FC = mean(log2(RQ)), ee = sd(log2(RQ)) / sqrt(n()), n = n(),
            .groups = "drop")

# --- puente y lado RNA-seq ---------------------------------------------------
puente <- t11 |> select(simbolo_manuscrito, locus = locus_tag, simbolo = simbolo_censo_189)
stopifnot(all(qpcr$simbolo_manuscrito %in% puente$simbolo_manuscrito),
          !anyDuplicated(puente$locus))

t10 <- qpcr |>
  left_join(puente, by = "simbolo_manuscrito") |>
  mutate(contraste_rnaseq = unname(EQUIV[contraste])) |>
  left_join(t36 |> select(locus, contraste_rnaseq = contraste,
                          rnaseq_log2FC = logFC, rnaseq_FDR = FDR),
            by = c("locus", "contraste_rnaseq")) |>
  transmute(simbolo_unico_EglNAC = simbolo, contraste, qpcr_log2FC, ee, n,
            rnaseq_log2FC, rnaseq_FDR, locus, simbolo_manuscrito, contraste_rnaseq) |>
  arrange(match(contraste, names(EQUIV)), simbolo_unico_EglNAC)

# si existe una T10 previa (cuantificada contra E. grandis) se guarda una copia
previo <- file.path(TABLAS, "T10_correlacion_qpcr_rnaseq.csv")
if (file.exists(previo) && !file.exists(file.path(DATOS, "T10_previo_Egrandis.csv")))
  file.copy(previo, file.path(DATOS, "T10_previo_Egrandis.csv"))
write_csv(t10, previo)

sin_de <- filter(t10, is.na(rnaseq_log2FC))
if (nrow(sin_de))
  message("  AVISO: ", nrow(sin_de), " pares sin estimacion de RNA-seq (no pasan filterByExpr): ",
          paste(unique(sin_de$simbolo_unico_EglNAC), collapse = ", "))

registrar_cifra("t10_pares", nrow(t10), "4 genes x 3 contrastes")
registrar_cifra("t10_pares_evaluables", sum(!is.na(t10$rnaseq_log2FC)),
                "con estimacion de RNA-seq contra E. globulus")
registrar_cifra("t10_genes_qpcr",
                paste(sprintf("%s (%s)", t10$simbolo_unico_EglNAC, t10$simbolo_manuscrito) |>
                        unique() |> sort(), collapse = "; "),
                "simbolo del censo (simbolo del manuscrito)")
volcar_cifras("correlacion_qpcr")
message("  T10 reescrita: ", nrow(t10), " pares, ",
        sum(!is.na(t10$rnaseq_log2FC)), " con RNA-seq")
