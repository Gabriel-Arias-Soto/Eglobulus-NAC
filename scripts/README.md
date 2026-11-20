# scripts/ — generacion de las figuras del manuscrito NAC

Cada figura del manuscrito se genera con un script de R, numerado en el orden en
que la figura aparece en el texto. Ningun script depende del estado de otro: todos
leen de `ANALISIS_Eglobulus/tablas/` y escriben en `ANALISIS_Eglobulus/figuras/`.

## Uso

```bash
cd "$NAC_ROOT"   # carpeta del proyecto (ver README.md de la raíz)
Rscript scripts/correr_todo.R          # el juego completo
Rscript scripts/01_filogenia.R         # una sola figura
```

Los scripts se ejecutan **desde la raiz del proyecto**: `00_estilo.R` toma
`getwd()` como raiz. `correr_todo.R` respalda los PNG anteriores en
`figuras/_v1_python/` antes de sobrescribir.

## Ficheros

| script | salida | insumos |
|---|---|---|
| `00_estilo.R` | — (se carga con `source`) | — |
| `20_subfamilias_Hussey.py` | `Hussey_FigS1_subfamilias_189.csv`; reescribe `subfamilia_final`, `ancla`, `dist_ancla` y `fiable` en `T01` y la subfamilia en `T02`, `T11`, `T14` y demas copias | `suppinfohussey/nph13139-sup-0001-figs1.pdf`, `Hussey_TablaS4_189simbolos.csv`, `arbol_NAC_189.treefile` |
| `01_filogenia.R` | `F06_filogenia_NAC` | `arbol_NAC_189.treefile`, `T01` |
| `02_subfamilias.R` | `F05_subfamilias_NAC` | `T01` |
| `03_arquitectura.R` | `F07_arquitectura_dominios` | `T16` |
| `04_motivos.R` | `F04_motivos_NAC` | `T04`, `T05`, `T01` |
| `05_cromosomas.R` | `F01_cromosomas_NAC` | `T01`, `T15`, `T10`, `datos/genes_todos.tsv` |
| `06_tandem.R` | `F09_vecindades_tandem` | `T15`, `T14`, `T01` |
| `16_quant_hoja.sh` | `rnaseq_udec/quant/SRR*/` | ENA `PRJNA486291`, `rnaseq_udec/salmon_idx` |
| `17_de_globulus.R` | `T33`, `T34`, `T35`, `T36`, `F12_qc_rnaseq_globulus` | `rnaseq_udec/quant/`, `t2g.tsv`, `T01` |
| `18_correlacion_qpcr.R` | `T10` (regenera), `T10_previo_Egrandis.csv` | `T08`, `T11`, `T36`, `T01` |
| `07_expresion.R` | `F02_heatmap_expresion`, `S01_heatmap_expresion_completo` | `T36`, `T01` |
| `08_qpcr.R` | `F03_correlacion_qpcr_rnaseq` | `T10` |
| `09_estructura.R` | `F08_monomero_EglNAC66` | `T17`, `figuras/_render_EglNAC66.png` |
| `10_auditoria.R` | `T21_auditoria_figuras_R.csv` | los PNG generados |

Cada script guarda **PNG a 300 dpi y PDF vectorial** con el mismo nombre. El
`.Rmd` inserta el PDF cuando compila a LaTeX y el PNG cuando compila a Word.

## Convenciones que fija `00_estilo.R`

1. **Un solo mapa de color para subfamilia** (`PAL_SUBFAM`, 19 niveles, los de la
   Fig. S1 de Hussey et al. 2015, mas gris para dominio parcial), usado en F05, F06 y F09. Antes la subfamilia Vb era
   azul marino en F05 y rosa en F06.
2. **Ningun titulo dentro de la imagen.** `theme_nac()` pone `plot.title` y
   `plot.subtitle` en blanco a proposito. El texto explicativo va en el `fig.cap`
   del `.Rmd`.
3. **Un solo ancho canonico:** `ANCHO_COLUMNA = 6.30` in (cuerpo de A4 con
   margenes de 2,5 cm) o `ANCHO_MEDIA = 3.05` in. Las figuras se generan ya al
   ancho final y no se reescalan al insertarlas, para que el cuerpo de letra sea
   el mismo en todas.
4. **Nomenclatura vigente:** `GEN_FOCAL = "EglNAC66"`. El simbolo `EglNAC34`
   esta retirado (`NOMENCLATURA_decision.md`); ningun script debe emitirlo.
5. **Leyenda interna abajo-derecha** (`leyenda_interna()`) o externa a la derecha
   en las figuras anchas.
6. **Rotulos de panel** `a`, `b` en negrita via `patchwork::plot_annotation`.

`guardar_figura()` avisa por consola si una figura supera el alto util de la
pagina (`ALTO_MAXIMO = 8.70` in).

## Dependencias

R 4.5 con `tidyverse`, `patchwork`, `ggrepel`, `scales`, `png`, `ape` y, de
Bioconductor, `ggtree` y `treeio`. En el proyecto estan en el entorno conda
`r-nac`.

## Notas

- Subfamilias: `python scripts/20_subfamilias_Hussey.py` muestra
  los cambios sin escribir; con `--escribir` respalda las tablas en
  `ANALISIS_Eglobulus/_respaldo_pre_subfamilias_Hussey/` y las reescribe. Se
  ejecuta antes de `correr_todo.R`. Las etiquetas de las anclas salen de la
  Fig. S1 de Hussey et al. (2015), no de la Tabla S10: donde ambas discrepan
  (EgrNAC137-139) se usa la figura.
- `01_filogenia.R` dibuja 341 puntas (185 de *E. globulus* con punta en el arbol
  de los 189 del censo, mas 156 anclas de *E. grandis*). Los valores de soporte
  estan en el `.treefile` y se omiten del grafico; para mostrarlos, poner
  `MOSTRAR_SOPORTE <- TRUE`.
- `07_expresion.R` parte el mapa de calor en dos: la figura principal lleva solo
  los genes con respuesta significativa (`UMBRAL_FDR`, `UMBRAL_LOGFC` editables
  en la cabecera) y el mapa completo va como `S01` a material suplementario.
  El mapa completo no cabe legible en una pagina.
- `09_estructura.R` anade el panel del dimero con ADN en cuanto exista
  `figuras/_render_dimero_ADN.png`. Ese panel paso a ser opcional: la imagen de
  la proteina sobre el ADN ya existe en F10 y no sale de una prediccion, sino de
  superponer el dominio sobre el co-cristal ANAC019-ADN (3SWP). Un
  co-plegamiento no mejoraria esa evidencia; como mucho la reproduciria.
- `13_estructura3D.py` es el unico script de figuras en Python y el unico que no
  lanza `correr_todo.R`, porque necesita PyMOL, que esta en el
  entorno conda `struct`, no en `r-nac`. Se ejecuta a mano desde la raiz:
  `pymol -cq scripts/13_estructura3D.py`. Produce
  `F10_estructura3D_<GEN_FOCAL>.png`, deja los renders sin componer en
  `figuras/_render3D/` y sus cifras en `tablas/_cifras/f10.csv`. Hereda las dos
  matrices de camara de `estructura_EglNAC34/vista_*.pml` para fijar la
  orientacion.
- La ortologia con el trio de estres de Arabidopsis va en dos pasos, porque el
  analisis necesita MAFFT/trimAl/IQ-TREE (entorno conda `filo`) y la figura
  necesita ggtree (entorno `r-nac`). Desde la raiz:
  `python scripts/14_ortologia_SNAC.py` y despues
  `Rscript scripts/15_arbol_SNAC.R`. El primero baja las anclas de UniProt,
  alinea, recorta al dominio NAM, infiere dos arboles (dominio y proteina
  completa), contrasta con la prueba AU la topologia que esperaria la filogenia
  de especies, y escribe `T30`-`T32` mas `tablas/_cifras/ortologia.csv`. El
  segundo dibuja `F11_ortologia_SNAC`.
- Dos cosas de `14_ortologia_SNAC.py` que no se deben tocar sin saber por que:
  `KMP_AFFINITY=disabled` (sin eso IQ-TREE puede abortar al iniciar OpenMP,
  p. ej. en WSL2) y `SEMILLA`, que se pasa a las tres invocaciones de IQ-TREE: sin
  semilla fija el soporte de ramas cambia de una corrida a otra. La p de la
  prueba AU es ademas una estimacion Monte Carlo y se mueve en la ultima cifra
  aunque la semilla este fija; por eso se corre con 100.000 replicas RELL y se
  reporta como aproximada.
- Las etiquetas de contraste estan en la constante `CONTRASTES` de
  `07_expresion.R` y en las de `08_qpcr.R`: si el manuscrito cambia de
  nomenclatura de contrastes, se editan ahi y en ningun otro sitio.
- **La expresion se cuantifica contra el genoma conspecifico.** Cuantificar
  contra *E. grandis* (loci `Eucgr.*`) obligaria a cruzar con el censo por el
  texto del simbolo, y un mismo simbolo nombra loci distintos en numeraciones
  distintas. `16` + `17` cuantifican las 12 librerias de hoja de `PRJNA486291`
  contra `GCA_046226855.1` y devuelven la expresion en loci `ACJRO7_*`, que es la
  clave primaria del censo; `18` rehace `T10` por la cadena cebador -> transcrito
  -> locus de `T11`. Ningun cruce se hace por nombre.
- `16_quant_hoja.sh` verifica el tamano de cada FASTQ contra el que declara ENA
  antes de cuantificar, reanuda la descarga de lo que falte y comprueba el md5.
  No es celo: siete de las doce descargas iniciales estaban truncadas, y un
  flujo gzip cortado limpiamente se cuantifica sin dar error, devolviendo un
  subconjunto de lecturas distinto en cada grupo. `17_de_globulus.R` repite la
  comprobacion por el otro lado (`num_processed` de salmon contra `read_count`
  de ENA) y se detiene si alguna libreria no cuadra.
- La banda de senal propia de `F02` es la fraccion de lecturas que salmon asigna
  de forma unica al locus, mediana sobre las 12 librerias (`T35`, umbral
  `UMBRAL_UNICA`).
