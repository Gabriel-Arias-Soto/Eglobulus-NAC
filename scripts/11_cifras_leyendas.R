# 11_cifras_leyendas.R -> tablas/T23_cifras_leyendas.csv
# Une los volcados que cada script de figura dejo en tablas/_cifras/ y escribe la
# tabla unica que lee el .Rmd. Ninguna leyenda lleva un numero teclado a mano:
# esa es la barrera contra lo que encontro la revision (F06 imprimia IVa 47 donde
# el censo da 53; F07 afirmaba 163 loci con cobertura completa y no se reproduce).
if (!exists("ANCHO_COLUMNA")) source("scripts/00_estilo.R")

censo <- leer_censo()
registrar_cifra("global_genes", nrow(censo), "T01_censo_maestro_189.csv")
registrar_cifra("global_subfamilias",
                dplyr::n_distinct(censo$subfamilia_final[censo$subfamilia_final != SUBFAM_NO_EVAL]),
                "T01.subfamilia_final")
registrar_cifra("global_gen_focal", GEN_FOCAL, "NOMENCLATURA_decision.md")
registrar_cifra("global_gen_focal_locus", GEN_FOCAL_LOCUS, "T01")
registrar_cifra("global_ancho_figura", sprintf("%.2f in", ANCHO_COLUMNA), "00_estilo.R")
volcar_cifras("global")

ar <- list.files(CIFRAS_DIR, pattern = "\\.csv$", full.names = TRUE)
if (!length(ar)) stop("no hay volcados en ", CIFRAS_DIR, ": corre antes los scripts de figura")
z <- dplyr::bind_rows(lapply(ar, function(f)
  suppressMessages(readr::read_csv(f, show_col_types = FALSE, col_types = readr::cols(.default = "c")))))
dup <- z$clave[duplicated(z$clave)]
if (length(dup)) message("  AVISO: claves repetidas: ", paste(unique(dup), collapse = ", "))
z <- z[!duplicated(z$clave), ]
z <- z[order(z$clave), ]
readr::write_csv(z, file.path(TABLAS, "T23_cifras_leyendas.csv"))
message(sprintf("  T23_cifras_leyendas.csv: %d cifras de %d scripts", nrow(z), length(ar)))
print(as.data.frame(table(sub("_.*", "", z$clave))))
