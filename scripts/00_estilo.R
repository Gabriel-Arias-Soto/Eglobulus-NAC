# 00_estilo.R -- reglas de estilo comunes a todas las figuras del manuscrito NAC.
# Se carga con source() al principio de cada script. No dibuja nada por si mismo.
# Convenciones que fija este fichero:
#   1. un unico mapa de color para subfamilia, usado en todas las figuras
#   2. ningun titulo dentro de la imagen: el texto explicativo va en el fig.cap del .Rmd
#   3. los dos anchos que admite Trees: 174 mm (caja completa) y 84 mm (columna)
#   4. nomenclatura vigente: el gen protagonista es EglNAC66 (NUNCA EglNAC34)

suppressPackageStartupMessages({
  library(ggplot2); library(dplyr); library(tidyr); library(readr)
  library(scales); library(grid)
})

if (!exists("RAIZ")) RAIZ <- getwd()
ANALISIS <- file.path(RAIZ, "ANALISIS_Eglobulus")
TABLAS   <- file.path(ANALISIS, "tablas")
FIGURAS  <- file.path(ANALISIS, "figuras")
DATOS    <- file.path(ANALISIS, "datos")
stopifnot(dir.exists(TABLAS), dir.exists(FIGURAS))

# --- anchos canonicos: norma de Trees (Springer) -----------------------------
# Springer fija DOS anchos de figura y un alto maximo, en milimetros. No son
# los del cuerpo de A4 que usabamos antes: las figuras se maquetan a la caja de
# la revista, no a la del borrador. Las medidas van en pulgadas porque es lo
# que toman ggsave() y png(), pero la cifra de referencia es la de mm.
MM <- 1 / 25.4
ANCHO_COLUMNA <- 174 * MM   # 6,850 in -- ancho completo de caja
ANCHO_MEDIA   <-  84 * MM   # 3,307 in -- una columna
ALTO_MAXIMO   <- 234 * MM   # 9,213 in -- alto maximo admitido
ALTO_PANELES  <- 8.60       # alto de los paneles de 189 filas (F07, F04 y FIG02)
# ALTO_PANELES no sube al nuevo techo a proposito: la geometria de esos paneles
# esta calibrada por alto_barra() para que las barras no se invadan, y el limite
# de la revista es un techo, no una altura exigida.

# --- tipografia --------------------------------------------------------------
# Trees pide Helvetica o Arial. Si el equipo las tiene (Windows, macOS) se usan
# tal cual; si no, se recurre a Liberation Sans, que es METRICAMENTE identica a
# Arial, y se registra desde el paquete fontLiberation. No se renombra a
# "Arial": el PDF debe declarar la fuente que lleva de verdad.
#
# fontconfig (14-10-2026). register_font() solo lo ven ragg y svglite. cairo_pdf()
# y png(type = "cairo") -los que usan guardar_figura() para el PDF y guardar_base()
# para todo- buscan la familia en fontconfig, que no conocia Liberation Sans y
# ponia DejaVu Sans en su lugar sin avisar: los PNG de ggplot salian en Liberation
# Sans y sus PDF en DejaVu Sans, mas ancha. Solo en Linux se anade a fontconfig la
# carpeta de fontLiberation, ANTES de la primera llamada a systemfonts, que es la
# que inicializa fontconfig en el proceso; despues ya no tendria efecto.
if (Sys.info()[["sysname"]] == "Linux" && nzchar(system.file(package = "fontLiberation")) &&
    basename(Sys.getenv("FONTCONFIG_FILE")) != "fonts_nac.conf") local({
  cand <- c(Sys.getenv("FONTCONFIG_FILE"),
            file.path(R.home(), "..", "..", "etc", "fonts", "fonts.conf"), "/etc/fonts/fonts.conf")
  base <- cand[nzchar(cand) & file.exists(cand)][1]
  if (is.na(base)) return(invisible())        # sin configuracion que heredar: no se toca
  conf <- file.path(tempdir(), "fonts_nac.conf")
  writeLines(c('<?xml version="1.0"?>', '<!DOCTYPE fontconfig SYSTEM "fonts.dtd">', "<fontconfig>",
               sprintf('  <include ignore_missing="yes">%s</include>', normalizePath(base)),
               sprintf("  <dir>%s</dir>",
                       system.file("fonts", "liberation-fonts", package = "fontLiberation")),
               "</fontconfig>"), conf)
  Sys.setenv(FONTCONFIG_FILE = conf)
})
FAMILIA_FIG <- local({
  hay <- unique(systemfonts::system_fonts()$family)
  if ("Arial" %in% hay) return("Arial")
  if ("Helvetica" %in% hay) return("Helvetica")
  # con el bloque anterior Liberation Sans ya es fuente del sistema, y volver a
  # registrarla da error
  if ("Liberation Sans" %in% hay) return("Liberation Sans")
  d <- system.file("fonts", "liberation-fonts", package = "fontLiberation")
  if (!nzchar(d)) stop("sin Arial/Helvetica ni fontLiberation: instala una de ellas")
  systemfonts::register_font(
    name       = "Liberation Sans",
    plain      = file.path(d, "LiberationSans-Regular.ttf"),
    bold       = file.path(d, "LiberationSans-Bold.ttf"),
    italic     = file.path(d, "LiberationSans-Italic.ttf"),
    bolditalic = file.path(d, "LiberationSans-BoldItalic.ttf"))
  "Liberation Sans"
})
# Comprobacion: la fuente que incrusta cairo_pdf() tiene que ser FAMILIA_FIG. Es la
# que fallo sin que nada lo detectara, asi que se mira en cada carga.
FUENTE_PDF <- local({
  if (!isTRUE(capabilities("cairo"))) return(NA_character_)
  f <- tempfile(fileext = ".pdf")
  grDevices::cairo_pdf(f, width = 1, height = 1, family = FAMILIA_FIG)
  grid::grid.text("Ag"); grDevices::dev.off()
  b <- readBin(f, "raw", file.info(f)$size); unlink(f)
  m <- regmatches(rawToChar(b[b != as.raw(0)]),
                  regexpr("/FontName */[A-Z]{6}\\+[A-Za-z0-9-]+", rawToChar(b[b != as.raw(0)]),
                          useBytes = TRUE))
  v <- if (length(m)) sub(".*\\+", "", m) else NA_character_
  if (is.na(v) || !grepl(gsub(" ", "", FAMILIA_FIG), v, fixed = TRUE))
    warning(sprintf("cairo_pdf() incrusta %s y no %s: los PDF no llevan la fuente de la revista",
                    v, FAMILIA_FIG), call. = FALSE)
  v
})

# --- altura del cuerpo de barra en los paneles de gggenes --------------------
# gggenes dibuja el cuerpo de cada gen con una altura ABSOLUTA (unit "mm"), no
# proporcional al paso entre filas. Con 189 filas el paso ronda 1 mm, asi que
# el valor heredado de 1,5 mm hacia que cada barra invadiera a sus vecinas y el
# panel se leyera como una mancha continua. La altura tiene que salir del
# numero de filas y del alto util, no fijarse a ojo.
#   reservado: pulgadas que no son panel (eje x, leyenda, separacion entre
#              facetas y margenes), medido sobre la figura compuesta.
#   ocupacion: fraccion del paso que ocupa la barra; el resto es el aire que
#              separa una secuencia de la siguiente.
alto_barra <- function(n_filas, alto_fig = ALTO_PANELES,
                       reservado = 1.45, ocupacion = 0.62) {
  paso_mm <- (alto_fig - reservado) * 25.4 / n_filas
  unit(paso_mm * ocupacion, "mm")
}
# Trees pide 300 ppp para medio tono, 600 para arte combinado (trazo mas
# rotulacion) y 1200 para trazo puro. Todas estas figuras llevan rotulacion
# sobre trazo o sobre medio tono, asi que la categoria que aplica es la
# combinada y el minimo exigible es 600.
DPI           <- 600

# --- nomenclatura -----------------------------------------------------------
GEN_FOCAL        <- "EglNAC66"          # simbolo vigente (NOMENCLATURA_decision.md)
GEN_FOCAL_LOCUS  <- "ACJRO7_023128"
GEN_FOCAL_SUBFAM <- "Vb"

# Codigos de condicion tal como los define Materiales y metodos: CON control
# regado y DS deficit hidrico en el ensayo de sequia; NA no aclimatada, CABF
# aclimatacion al frio antes de la helada y CAAF despues de la helada en el
# ensayo de frio. Las etiquetas legibles se usan en F02 y F03.
ETQ_CONTRASTE <- c(
  "C_DS"      = "CON \u2192 DS\nwater deficit",
  "NA_DS"     = "NA \u2192 DS\nwater deficit",
  "NA_CABF"   = "NA \u2192 CABF\ncold, before frost",
  "NA_CAAF"   = "NA \u2192 CAAF\ncold, after frost",
  "CABF_CAAF" = "CABF \u2192 CAAF\nfrost effect")
etq_contraste <- function(x) {
  y <- unname(ETQ_CONTRASTE[as.character(x)])
  ifelse(is.na(y), as.character(x), y)
}

# --- paleta unica de subfamilia ---------------------------------------------
# 19 subfamilias con asignacion + un nivel gris para los dominios parciales.
# Las etiquetas son las de Hussey et al. (2015), Fig. S1, transferidas por el
# ancla mas proxima (scripts/20_subfamilias_Hussey.py); de ahi Va(1)/Va(2) y
# VIIIa. Colores escogidos sobre las paletas de Okabe-Ito y Tol, discriminables
# en deuteranopia y en impresion en gris. Va(1) conserva el color de la antigua
# Va; Va(2) va en ciruela oscuro (misma familia de tono, dE >= 24 frente a todos
# los demas tambien con deuteranopia simulada). VIIIa va en el salmon de la
# paleta clara de Tol (dE >= 13,5 frente a todo lo fijo, incluidos los grises y
# el rojo focal; el minimo es IVa, al otro lado del arbol). Se descarto el negro
# porque los nodos con UFBoot >= 95 de F06 se marcan en grey25, y el marino
# porque se confundia con Va(2) en deuteranopia (dE 4,6).
SUBFAM_NIVELES <- c("Ia","Ib","Ic","II","IIIa/b","IIIc","IVa","IVb","IVc","IVd",
                    "Va(1)","Va(2)","Vb","VIa","VIb","VIc","VII","VIIIa","XI")
SUBFAM_COLORES <- c("Ia"="#332288","Ib"="#6699CC","Ic"="#88CCEE","II"="#44AA99",
                    "IIIa/b"="#117733","IIIc"="#999933","IVa"="#DDCC77","IVb"="#E69F00",
                    "IVc"="#CC6677","IVd"="#882255","Va(1)"="#AA4499","Va(2)"="#4B0055",
                    "Vb"="#CC79A7","VIa"="#0072B2","VIb"="#56B4E9","VIc"="#009E73",
                    "VII"="#D55E00","VIIIa"="#EE8866","XI"="#661100")
# Las tablas se escribieron con esta categoria en espanol y las figuras van en
# ingles. Se traduce al leer y no se reescribe T01: la tabla es el dato que se
# publica y su contenido no cambia porque cambie el idioma de una figura.
SUBFAM_NO_EVAL_CSV <- "no evaluada (dominio parcial)"
SUBFAM_NO_EVAL     <- "not assessed (partial domain)"
ETQ_ANCLA  <- "E. grandis (anchor)"   # las puntas de la especie de referencia
ETQ_YA_MS    <- "already in the manuscript"
ETQ_NUEVO_MS <- "new in the 189-gene census"
# El gris del grupo sin evaluar es MEDIO, no claro: es una categoria con datos
# (6 genes del censo) y a #BDBDBD se leia como fondo o como dato ausente.
PAL_SUBFAM <- c(SUBFAM_COLORES, setNames("#8C8C8C", SUBFAM_NO_EVAL),
                setNames("#9E9E9E", ETQ_ANCLA))

escala_subfam <- function(...) list(
  scale_colour_manual(values = PAL_SUBFAM, na.value = "#BDBDBD",
                      breaks = c(SUBFAM_NIVELES, SUBFAM_NO_EVAL), name = "subfamily", ...),
  scale_fill_manual(values = PAL_SUBFAM, na.value = "#BDBDBD",
                    breaks = c(SUBFAM_NIVELES, SUBFAM_NO_EVAL), name = "subfamily", ...))

# --- paletas de dominio y subdominio ---------------------------------------
# La region C-terminal va en gris claro a proposito: ocupa la mayor parte del
# area de F07 y cualquier color saturado ahi domina la figura. El salmon queda
# reservado a la banda de desorden de F08, para que no se confundan.
# Las clases de dominio y de subdominio viajan como constantes y no como texto
# suelto: F07 y F04 las nombran en la leyenda y en los datos, y con literales
# repartidos por dos scripts bastaba una tilde de diferencia para que la clase
# cayera en na.value y se quedara sin clave.
DOM_NAM <- "NAM domain (PF02365)"
DOM_TR  <- "C-terminal region (TR)"
DOM_NTL <- "transmembrane helix (NTL)"
DOM_TM  <- "TM helix outside NTL criterion"
PAL_DOMINIO <- setNames(c("#2166AC", "#E0E0E0", "#B2182B", "#E8A0A8"),
                        c(DOM_NAM, DOM_TR, DOM_NTL, DOM_TM))
SUBDOM_NA   <- "no subdomain assigned"
PAL_SUBDOM  <- c("A"="#08519C","B"="#3182BD","C"="#6BAED6","D"="#FD8D3C","E"="#E6550D",
                 setNames("#BDBDBD", SUBDOM_NA))
ETQ_EJE_AA  <- "position in protein (aa)"   # eje comun de F07 y F04
# T04.subdominio guarda esta clase en espanol: misma regla que SUBFAM_NO_EVAL,
# se traduce al leer y la tabla publicada no se toca.
MOTIVO_CTERM_CSV <- "C-terminal / espec\u00edfico de subfamilia"
MOTIVO_CTERM     <- "C-terminal / subfamily-specific"

# Rotulo corto del grupo sin evaluar para las BANDAS de faceta. La clave larga
# cabe en una leyenda pero no en una banda lateral, y dejarla caer a NA no es
# opcion: en este juego de figuras "NA" ya es el nombre de una condicion del
# diseno (NA -> DS, NA -> CABF), asi que una banda rotulada NA se lee mal.
ETQ_PARCIAL <- "partial domain"
subfam_facetas <- function(x) {
  x <- as.character(x)
  x[x %in% SUBFAM_NO_EVAL_CSV] <- SUBFAM_NO_EVAL   # T16 y T01 guardan el valor en espanol
  z <- factor(x, levels = c(SUBFAM_NIVELES, SUBFAM_NO_EVAL),
              labels = c(SUBFAM_NIVELES, ETQ_PARCIAL))
  stopifnot(!any(is.na(z)))
  z
}
COL_DESTACADO <- "#B2182B"   # para marcar el gen protagonista
COL_NEUTRO    <- "#BDBDBD"

# --- tema comun de las dos figuras genomicas (F01 y F09) --------------------
# F01 dice DONDE esta cada bloque en tandem y QUE EglNAC lo forman; F09 despliega
# la estructura interna de ese mismo bloque. El codigo B01-B21 de T15 es la clave
# que las une, asi que las dos comparten color, separador y orden del rotulo: el
# lector salta de una a otra por el codigo. Cambiar el tema aqui las cambia a las
# dos; era justo lo que fallaba cuando cada script fijaba su propio gris.
COL_IDEOGRAMA    <- "grey88"   # cuerpo del cromosoma en F01
COL_NAC          <- "grey45"   # un NAC cualquiera, en las dos figuras
COL_BLOQUE       <- "#2166AC"  # bloque en tandem: marca en F01, cabecera en F09
COL_BLOQUE_SUAVE <- "#DCE6F1"
ETQ_NO_NAC       <- "non-NAC gene within block"
SEP_BLOQUE       <- " \u00b7 "  # separador unico de los rotulos de bloque

tema_genomico <- function(base = 8) theme_nac(base) +
  theme(strip.background = element_rect(fill = COL_BLOQUE_SUAVE, colour = NA),
        strip.text = element_text(size = base - 0.6, face = "bold", colour = "grey15"),
        panel.grid = element_blank())

# Numeros de EglNAC de un bloque, colapsando las rachas consecutivas: c(45:55)
# da "45-55" y c(28,29,31) da "28-29, 31". Sin esto el rotulo de B01 serian once
# numeros separados por comas y no cabria en la anchura de una columna.
rango_compacto <- function(n, max_trozos = 3L) {
  n <- sort(unique(as.integer(n)))
  if (!length(n)) return("")
  trozos <- vapply(split(n, cumsum(c(1, diff(n) != 1))), function(v)
    if (length(v) == 1L) as.character(v) else sprintf("%d\u2013%d", min(v), max(v)),
    character(1))
  if (length(trozos) > max_trozos)
    trozos <- c(trozos[seq_len(max_trozos)], "\u2026")
  paste(trozos, collapse = ", ")
}

# --- tema -------------------------------------------------------------------
# plot.title y plot.subtitle van en blanco a proposito: hacen cumplir la regla 2.
# Trees exige rotulacion de 8 a 12 pt y "variacion minima de tamano dentro de
# una misma figura". theme_bw() deja axis.text en rel(0.8) = 7,2 pt con base 9,
# por debajo del suelo, asi que se fija explicitamente. Con base 9 el rango
# entero de la figura queda en 8-12 pt: 8 (ejes, leyenda), 8,5 (facetas),
# 9 (titulos de eje) y 12 solo la letra de panel.
theme_nac <- function(base = 9) {
  theme_bw(base_size = base, base_family = FAMILIA_FIG) +
    theme(panel.grid.minor  = element_blank(),
          panel.grid.major.x = element_blank(),
          panel.grid.major.y = element_line(colour = "grey92", linewidth = 0.25),
          panel.border = element_rect(colour = "grey30", linewidth = 0.4),
          axis.ticks   = element_line(colour = "grey30", linewidth = 0.3),
          axis.text    = element_text(size = base - 1),
          axis.title   = element_text(size = base),
          strip.background = element_rect(fill = "grey95", colour = NA),
          strip.text   = element_text(size = base - 0.5, face = "bold"),
          legend.key.size = unit(0.34, "cm"),
          legend.title = element_text(face = "bold", size = base - 1),
          legend.text  = element_text(size = base - 1),
          legend.background = element_rect(fill = alpha("white", 0.85), colour = NA),
          plot.title = element_blank(), plot.subtitle = element_blank(),
          plot.tag = element_text(face = "bold", size = base + 3),
          plot.margin = margin(2, 3, 2, 2))
}
# leyenda interna abajo-derecha (regla 5)
leyenda_interna <- function() theme(legend.position = c(0.99, 0.01),
                                    legend.justification = c(1, 0))

# --- guardado ---------------------------------------------------------------
guardar_figura <- function(p, nombre, alto, ancho = ANCHO_COLUMNA, pdf = TRUE) {
  if (alto > ALTO_MAXIMO)
    warning(sprintf("%s mide %.2f in de alto: desborda el alto util (%.2f in)",
                    nombre, alto, ALTO_MAXIMO), call. = FALSE)
  ggsave(file.path(FIGURAS, paste0(nombre, ".png")), p, width = ancho,
         height = alto, dpi = DPI, units = "in", bg = "white")
  if (pdf) {
    disp <- if (isTRUE(capabilities("cairo"))) grDevices::cairo_pdf else "pdf"
    ggsave(file.path(FIGURAS, paste0(nombre, ".pdf")), p, width = ancho,
           height = alto, units = "in", device = disp, bg = "white")
  }
  message(sprintf("  guardada %-34s %.2f x %.2f in", nombre, ancho, alto))
  invisible(p)
}

# --- lectores de datos ------------------------------------------------------
leer_censo <- function() {
  z <- suppressMessages(read_csv(file.path(TABLAS, "T01_censo_maestro_189.csv"),
                                 show_col_types = FALSE))
  z$subfamilia_final <- factor(ifelse(z$subfamilia_final == SUBFAM_NO_EVAL_CSV,
                                      SUBFAM_NO_EVAL, z$subfamilia_final),
                               levels = c(SUBFAM_NIVELES, SUBFAM_NO_EVAL))
  stopifnot(!any(is.na(z$subfamilia_final)))
  z$en_manuscrito <- as.logical(z$en_manuscrito)
  z$es_focal <- z$simbolo == GEN_FOCAL
  stopifnot(nrow(z) == 189, sum(z$es_focal) == 1)
  z
}

# Genes medidos por qRT-PCR. La fuente es T02.medido_por_qPCR y no la tabla de
# correlacion: asi el conjunto es el mismo en F01, F02 y F03. Se comprueba que
# T10 cubre exactamente ese conjunto.
genes_qpcr <- function() {
  z <- suppressMessages(read_csv(file.path(TABLAS, "T02_equivalencia_nomenclaturas_189.csv"),
                                 show_col_types = FALSE))
  v <- sort(z$simbolo[z$medido_por_qPCR %in% TRUE])
  t10 <- suppressMessages(read_csv(file.path(TABLAS, "T10_correlacion_qpcr_rnaseq.csv"),
                                   show_col_types = FALSE))
  if (!setequal(v, unique(t10$simbolo_unico_EglNAC)))
    message("  AVISO: T02.medido_por_qPCR y T10 no cubren el mismo conjunto de genes")
  v
}
VALIDADOS <- genes_qpcr()
stopifnot(length(VALIDADOS) >= 1, GEN_FOCAL %in% VALIDADOS)

# mapa de accesiones del ensamblado a numero de cromosoma (orden de accesion)
mapa_cromosomas <- function() {
  g <- suppressMessages(read_tsv(file.path(DATOS, "genes_todos.tsv"),
        col_names = c("scaf","ini","fin","hebra","locus","biotipo"), show_col_types = FALSE))
  cm <- sort(unique(g$scaf[grepl("^CM", g$scaf)]))
  list(genes = g, mapa = setNames(paste0("chr", seq_along(cm)), cm),
       largos = g |> filter(scaf %in% cm) |> group_by(scaf) |>
         summarise(largo = max(fin), .groups = "drop") |>
         mutate(cromosoma = paste0("chr", match(scaf, cm))))
}

# --- guardado de figuras que no son ggplot (ComplexHeatmap, grid, base) -----
# 'dibujar' es una funcion sin argumentos que pinta en el dispositivo activo.
# suplementaria = TRUE: la figura no va en el cuerpo, asi que el alto util de
# pagina no la ata y el desborde es una decision, no un descuido. Se avisa igual
# pero como mensaje, para que el aviso de verdad siga destacando.
guardar_base <- function(dibujar, nombre, alto, ancho = ANCHO_COLUMNA, pdf = TRUE,
                         suplementaria = FALSE) {
  if (alto > ALTO_MAXIMO) {
    aviso <- sprintf("%s mide %.2f in de alto: desborda el alto util (%.2f in)",
                     nombre, alto, ALTO_MAXIMO)
    if (suplementaria) message("  ", aviso, " [suplementaria]")
    else warning(aviso, call. = FALSE)
  }
  # family: ComplexHeatmap y el resto de graficos base no pasan por theme_nac(),
  # asi que la tipografia de la revista hay que darsela al dispositivo.
  png(file.path(FIGURAS, paste0(nombre, ".png")), width = ancho, height = alto,
      units = "in", res = DPI, bg = "white", family = FAMILIA_FIG,
      type = if (isTRUE(capabilities("cairo"))) "cairo" else "Xlib")
  dibujar(); dev.off()
  if (pdf && isTRUE(capabilities("cairo"))) {
    grDevices::cairo_pdf(file.path(FIGURAS, paste0(nombre, ".pdf")),
                         width = ancho, height = alto, bg = "white",
                         family = FAMILIA_FIG)
    dibujar(); dev.off()
  }
  message(sprintf("  guardada %-34s %.2f x %.2f in", nombre, ancho, alto))
  invisible(nombre)
}

# --- objetos intermedios para la panelizacion -------------------------------
# Cada script de figura guarda su ggplot aqui; 12_paneles.R los compone en las
# seis figuras del manuscrito. Asi la panelizacion vive en un solo sitio.
OBJETOS <- file.path(FIGURAS, "_objetos")
dir.create(OBJETOS, showWarnings = FALSE)
guardar_objeto <- function(p, nombre) {
  saveRDS(p, file.path(OBJETOS, paste0(nombre, ".rds"))); invisible(p)
}
leer_objeto <- function(nombre) {
  f <- file.path(OBJETOS, paste0(nombre, ".rds"))
  if (!file.exists(f)) stop("falta el objeto ", nombre, ": corre antes su script")
  readRDS(f)
}

# --- registro de cifras para las leyendas -----------------------------------
# Ninguna cifra se teclea en una leyenda: cada script registra las suyas aqui y
# el .Rmd las lee de T23_cifras_leyendas.csv con cif("clave"). Es la barrera
# contra lo que encontro la revision: numeros impresos que ya no cuadran con las
# tablas (F06 daba IVa 47 donde el censo da 53; F07 afirmaba 163 loci).
.CIFRAS <- new.env(parent = emptyenv())
registrar_cifra <- function(clave, valor, fuente, nota = NA_character_) {
  assign(clave, list(clave = clave, valor = as.character(valor),
                     fuente = fuente, nota = nota), envir = .CIFRAS)
  invisible(valor)
}
# Cada script vuelca SUS cifras a tablas/_cifras/<prefijo>.csv y 11_cifras_leyendas.R
# las une. Asi cada script puede correr en su propio proceso de R: cargar un
# paquete de Bioconductor enmascara verbos de dplyr (S4Vectors exporta rename)
# y esa contaminacion no debe cruzar de un script al siguiente.
CIFRAS_DIR <- file.path(TABLAS, "_cifras")
dir.create(CIFRAS_DIR, showWarnings = FALSE)
volcar_cifras <- function(prefijo) {
  k <- ls(.CIFRAS)
  if (!length(k)) return(invisible(NULL))
  z <- do.call(rbind, lapply(k, function(i)
    as.data.frame(get(i, envir = .CIFRAS), stringsAsFactors = FALSE)))
  readr::write_csv(z[order(z$clave), ], file.path(CIFRAS_DIR, paste0(prefijo, ".csv")))
  message(sprintf("  %d cifras -> _cifras/%s.csv", nrow(z), prefijo))
  invisible(z)
}

escribir_cifras <- function(archivo = file.path(TABLAS, "T23_cifras_leyendas.csv")) {
  k <- ls(.CIFRAS)
  if (!length(k)) { warning("no hay cifras registradas", call. = FALSE); return(invisible(NULL)) }
  z <- do.call(rbind, lapply(k, function(i)
    as.data.frame(get(i, envir = .CIFRAS), stringsAsFactors = FALSE)))
  z <- z[order(z$clave), ]
  write_csv(z, archivo)
  message(sprintf("  escritas %d cifras en %s", nrow(z), basename(archivo)))
  invisible(z)
}

message(sprintf("00_estilo.R cargado | Trees: %.0f/%.0f mm, %d ppp, %s | gen focal %s",
                ANCHO_COLUMNA/MM, ANCHO_MEDIA/MM, DPI, FAMILIA_FIG, GEN_FOCAL))
