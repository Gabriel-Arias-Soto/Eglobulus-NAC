# 12_paneles.R -> figuras compuestas del manuscrito
# La panelizacion vive SOLO aqui: cada script de figura guarda su ggplot en
# figuras/_objetos/ y este script los compone. Asi cambiar el reparto de paneles
# no obliga a tocar los scripts de figura.
# F02 (ComplexHeatmap) y F08 (ya compuesta en 09_estructura.R) van aparte porque
# no son ggplots componibles con patchwork.
if (!exists("ANCHO_COLUMNA")) source("scripts/00_estilo.R")
suppressPackageStartupMessages(library(patchwork))

etq <- function(p) p + plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(face = "bold", size = 10))

# --- Figura 1: clasificacion filogenetica ------------------------------------
# Un solo panel: el arbol. El recuento por subfamilia (P_subfamilias) sale de la
# figura del cuerpo y queda como figura suelta (F05), porque el color de rama del
# arbol ya da la lectura por subfamilia que aportaba ese panel.
guardar_figura(leer_objeto("P_filogenia"), "FIG01_clasificacion", alto = 5.30)

# --- Figura 2: arquitectura del dominio --------------------------------------
fig2 <- leer_objeto("P_arquitectura") | leer_objeto("P_motivos")
guardar_figura(etq(fig2), "FIG02_arquitectura", alto = ALTO_PANELES)

message("  nota: F01, F05, F09, F02 y F08 van como figuras de un solo panel")
