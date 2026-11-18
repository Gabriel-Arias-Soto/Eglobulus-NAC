# NAC transcription factor family of *Eucalyptus globulus*

Analysis and figure scripts for the article *Genome-wide analysis of the NAC family in Eucalyptus globulus identifies
EglNAC66 as a candidate regulator of drought and cold responses* (Aguayo et al., submitted to *Tree Genetics & Genomes*).

The repository contains code only. Comments inside the scripts are in Spanish.

## Layout

The folders mirror the layout of the project directory, so that the relative paths used by the scripts keep working.

| Folder | Content |
|---|---|
| `scripts/` | Figures and tables of the main analysis (phylogeny, subfamilies, domain architecture, motifs, chromosomes, tandem blocks, leaf RNA-seq quantification and differential expression, qRT-PCR, structure, SNAC orthology, subfamily assignment from Hussey et al. 2015). `scripts/correr_todo.R` runs the R figure set; `scripts/README.md` maps each script to its inputs and outputs. |
| `ANALISIS_Eglobulus/datos/entrega3/scripts/` | Genome-wide analyses: representative isoforms, 2-kb promoters, PlantCARE and JASPAR scans against a GC-matched background, A. thaliana homologs and GO enrichment (clusterProfiler), collinearity and duplication modes (MCScanX, JCVI), Ka/Ks and dating, expansion by subfamily, and the qRT-PCR figures recomputed from raw Cq values. |
| `rnaseq_udec/` | Download of the leaf RNA-seq libraries (ENA, BioProject PRJNA486291), HISAT2 alignment in batches, the RNA-seq check of the NAC genes with a partial domain (`alineamiento_parciales/`) and the review of the 189 gene models with miniprot and splice junctions (`revision_modelos_189/`). |
| `arbol_185/` | Corrected gene models of EglNAC37 and EglNAC87, their addition to the domain alignment, the maximum-likelihood tree with 185 *E. globulus* proteins, the read check of the possible pseudogenes and the motif scan of the corrected proteins. |
| `docs/` | Software versions used for the article. |

## Running

Scripts are run from the project directory. Two environment variables replace the paths of the original workstation:

- `NAC_ROOT`: project directory (default: the current directory).
- `CONDA_ENVS`: directory that holds the conda environments (default: `~/conda/envs`).

Programs are expected on the `PATH` of the active environment. Software and versions are listed in
`docs/software_versions.csv` and in Methods S1 of the article; the analyses also use DIAMOND, MCScanX, JCVI, PAML
(yn00), HISAT2, SAMtools, miniprot and the MEME suite.

## Data

No data are included. Public inputs:

- *E. globulus* genome assembly and annotation GCA_046226855.1 (NCBI).
- *E. grandis* v2 (GCA_000612305.2), *Populus trichocarpa* v4.1 (GCF_000002775.5) and *Arabidopsis thaliana* TAIR10 (GCF_000001735.4).
- Leaf RNA-seq reads, BioProject PRJNA486291 (NCBI SRA / ENA).
- Pfam profile PF02365 (NAM).
- Subfamily designations and supplementary tables of Hussey et al. (2015).

The raw qRT-PCR Cq files belong to the laboratory that ran the assays and are available from the authors on request.

## Notes

- The 339-sequence domain alignment of the first tree (183 *E. globulus* proteins and 156 *E. grandis* anchors) was built
  interactively with hmmsearch, MAFFT L-INS-i and trimAl, as described in Methods. `arbol_185/` documents how the two
  corrected domains were added to it.
- `ANALISIS_Eglobulus/datos/entrega3/scripts/e3_19_controles.py` includes checks specific to the original workstation.

## License

MIT (see `LICENSE`).
