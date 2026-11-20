# correr_todo.R -> regenera el juego de figuras completo desde las tablas.
# Uso, desde la raiz de "Revision NAC":   Rscript scripts/correr_todo.R
# Cada script corre en SU PROPIO proceso de R a proposito: cargar GenomicRanges o
# ComplexHeatmap enmascara verbos de dplyr (S4Vectors exporta rename) y en una
# sola sesion esa contaminacion se arrastra al script siguiente. El traspaso de
# informacion entre scripts va por disco: figuras/_objetos/*.rds para los paneles
# y tablas/_cifras/*.csv para las cifras de las leyendas.
# 17 y 18 rehacen la expresion contra el genoma conspecifico y necesitan que
# antes haya corrido  bash scripts/16_quant_hoja.sh  (descarga las 12 librerias
# de hoja de PRJNA486291, verifica md5 y las cuantifica con salmon; ~2 h la
# primera vez, reanudable). Si no hay cuantificacion, 17 se detiene y lo dice.
# 14/15 (ortologia SNAC) van aparte: 14 corre en el entorno conda con mafft,
# trimal e iqtree, no en R.
PASOS <- c("01_filogenia", "02_subfamilias", "03_arquitectura", "04_motivos",
           "05_cromosomas", "06_tandem", "17_de_globulus", "18_correlacion_qpcr",
           "07_expresion", "08_qpcr",
           "09_estructura", "15_arbol_SNAC",
           "12_paneles", "11_cifras_leyendas", "10_auditoria")

RSCRIPT <- file.path(R.home("bin"), "Rscript")
# Las cifras de 13 y 14 las escriben scripts de Python que no estan en PASOS
# (necesitan pymol / mafft-iqtree y sus propios entornos). Borrarlas aqui dejaria
# a 15_arbol_SNAC sin insumo y a T23 sin las filas de ortologia y estructura.
# f10.csv lo escribe 13_estructura3D.py, que se lanza a mano con el PyMOL del
# entorno struct.
CIFRAS_EXTERNAS <- c("ortologia.csv", "f10.csv")
.cif <- list.files("ANALISIS_Eglobulus/tablas/_cifras", full.names = TRUE)
unlink(.cif[!basename(.cif) %in% CIFRAS_EXTERNAS])
t0 <- Sys.time(); fallos <- character()
for (s in PASOS) {
  message("\n== ", s, " ==")
  cod <- system2(RSCRIPT, c(shQuote(file.path("scripts", paste0(s, ".R")))))
  if (cod != 0) { fallos <- c(fallos, s); message("  FALLO en ", s, " (codigo ", cod, ")") }
}
# Las versiones de Metodos se sondean del entorno, no se teclean. Va fuera del
# bucle porque es un .sh y porque no alimenta ninguna figura: alimenta el
# manuscrito. Si falla no se aborta la corrida, solo se avisa.
message("\n== 19_versiones ==")
if (system2("bash", c("scripts/19_versiones.sh")) != 0)
  message("  AVISO: no se pudo regenerar manuscrito/versiones_software.csv")

message(sprintf("\nlisto en %.1f min | %d fallos%s",
                as.numeric(difftime(Sys.time(), t0, units = "mins")), length(fallos),
                if (length(fallos)) paste0(": ", paste(fallos, collapse = ", ")) else ""))
