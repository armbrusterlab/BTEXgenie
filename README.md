# BTEXgenie: A curated and user-friendly tool for profile HMM-based substrate-specific annotation of BTEX degradation genes

## Table of Contents
* [Installation](#installation)
* [KOfam database download](#kofam-database-download)
* [Usage](#usage)
* [KEGG Pathway visualizations](#kegg-Pathway-visualizations)
* [Circos plot visualizations](#circos-plot-visualizations)

## Installation

This toolkit uses hmmscan from the HMMER suite to query BTEXgenie profiles.
The following packages are needed for visualization and analysis:

- python 3.9 or newer
- pip
- Rscript
- hmmer 3.3 or newer
- prodigal
- kofamscan
- pandas 2.0 or newer
- numpy 1.24 or newer
- matplotlib 3.7 or newer
- biopython 1.81 or newer
- circos 0.69 or newer

If your system already satisfies these requirements, you can move directly to running the BTEXgenie scripts. Otherwise, you may install everything through Conda as shown below.

1. Conda Install
Confirm that a working Conda installation is available. See [Conda installation](https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html) for more details.

<!-- ### Development install
Create the BTEXgenie environment and install all dependencies plus an editable install of the local repo via *btex_env.yml*: -->

```bash
cd /path/to/BTEXgenie
conda env create -n BTEXgenie -f btex_env.yml
```

Activate the environment:

```bash
conda activate BTEXgenie
```

Confirm that the KOfamScan executable is installed:

```bash
command -v exec_annotation
```

KOfamScan should be available on the active conda environment PATH when the environment is set up correctly.

## KOfam Database download:

The KOfam HMM database can be installed for users interested in the broad metabolic or degradation potential of their genomes.

Run the following command to download KOfam database:

```bash
btex-build-db --db-dir /path/to/databases/kofam
```

This command will:

- download `ko_list.gz` and `profiles.tar.gz`
- extract them into the requested database directory
- verify that `profiles/` and `ko_list` were created successfully
- write `$CONDA_PREFIX/etc/conda/activate.d/kofam.sh` so `KOFAM_DB` is exported automatically when Conda environment is activated

If the target Conda environment is not currently active, provide its prefix with `--conda-prefix`:

```bash
btex-build-db --db-dir /path/to/databases/kofam --conda-prefix /path/to/conda/env
```

Confirm the database path and executable:
```bash
conda deactivate
conda activate BTEXgenie
echo $KOFAM_DB
command -v exec_annotation
```

`KOFAM_DB` should point to the database directory containing `profiles/` and `ko_list`.

## Usage
To run BTEXgenie, input can be either a directory or single file containing genome DNA FASTA or protein FASTA input.

### Example with protein files in *test_genomes*
```bash
btex-annotate -g btexhmm/test_genomes/protein_fastas \
              -o path/to/output_dir
```
> [!NOTE]
> For proper parsing of genomic coordinates, protein files produced from Prodigal are needed.

### Example with genome FASTA files
```bash
btex-annotate -g btexhmm/test_genomes/dna_fastas \
              -o path/to/output_dir \
              --kofam
```

> [!NOTE]
> By default, `btex-annotate` runs the BTEX HMM search together with the putative-target BLASTP search. Use `--kofam` to additionally run the KOfam annotation.
> The bundled HMM database used by default is `btexhmm/hmms/BTEXgenie_updated.hmm`; the selected path is recorded in `btex_annotate.log` at the start of every run.
>
> For FASTA sequence inputs, the program will run gene-calling with Prodigal in `--single` mode as default, unless the `--meta` specified as input.

**Main outputs**

1. `btex_genie_summary.csv`
   Reports individual BTEX HMM hits, including the matched HMM, the threshold used, the domain evalue, the hit score, and the protein sequence header for the corresponding gene.

2. `btex_genie_summary_counts.csv`
   Summarizes BTEX HMM hits by HMM, reporting hit counts instead of individual protein sequence headers.

3. 'btex_putative_blastp_unfiltered_summary.csv'
   Reports BLASTp alignments between input protein sequences and putative enzyme targets with an E-value below 1 x 10^-5

4. 'btex_putative_blastp_filtered_summary.csv'
   Reports BLASTp alignments to putative enzyme targets that meet the user defined E-value, sequence identity, and alignment coverage thresholds. 

5. `prodigal_output/`
   Generated when genome DNA sequences are used as input. This directory contains one subdirectory per genome and contains:
   ` {genome}_prodigal.gbk ` (not produced when the input consists of protein FASTAs)
   ` {genome}.faa `
   ` {genome}_kofam_abv_thres.tsv `, produced when `--kofam` is enabled and the KOfam step runs successfully for that sample.

6. `hmmscan_output/`
   Contains one subdirectory per input file with raw `.domtblout` results generated before and after filtering by GA thresholds.

7. `btex_annotate.log`
   Records run progress as well as detailed warning and error messages.

## KEGG Pathway visualizations
BTEXgenie supports visualization of hits to BTEXgenie or KOfam HMMs on KEGG pathways with HTML files.

**For visualization of BTEXgenie hits on BTEX-associated KEGG pathways:**

```bash
btex-vis --hmmscan /path/to/output_dir/btex_hmm_summary.csv \
         -o /path/to/btex-vis-outputs
```

> [!Note]
> Using -s {genome} generates a visualization for hits from a specific genome. The {genome} value must exactly match either the sample name in btex_hmm_summary.csv or the genome prefix of the corresponding {genome}.faa file.

**Example of BTEXgenie hits using genomes in test_genomes/protein_fastas:**

**(BTEXgenie genes are outlined in red; samples are color-coded)**
<p align="center">
  <img src="img/ethylbenzene_degra_BTEXgenie.png" width="800">
</p>


**For visualization of KOfam hits on a KEGG pathway:**

```bash
btex-vis -g /path/to/prodigal_output \
         -o /path/to/btex-vis-outputs \
         --pathways 00642
```
> [!Note]
> btex-vis takes the `prodigal_output` directory with `{genome}_kofam_abv_thres.tsv` per input genome instead of `btex_hmm_summary.csv` for visualizing all KOfam hits on an interested pathway.
>
> The `prodigal_output` directory is created by`btex-annotate` with the `--kofam` flag.

**Example of KOfam hits on ethylbezene degradation pathway:**
<p align="center">
  <img src="img/ethylbenzene_degra_KOfam.png" width="800">
</p>


**Outputs:**
1. `{output_dir}/{pathway_name}_{pathway_ID}.html`

   The main HTML file that can be opened to view hits on KEGG pathways.

   By default, if `--pathways` is not supplied, BTEXgenie generates HTML files for the following BTEX-associated KEGG pathways:

   * `xylene_degradation_00642.html`
   * `toluene_degradation_00623.html`
   * `ethylbenzene_degradation_00622.html`
   * `benzoate_degradation_00362.html`

   > [!NOTE]
   > An HTML file is not generated for a pathway if no hits are detected for that pathway.

   Users can also provide additional KEGG pathway IDs separated by comma with `--pathways`.
   A full list of KEGG pathways is available [here](https://www.genome.jp/kegg/pathway.html#energy).

2. {output_dir}/`sample_color_legend.tsv`
   Contains the color assigned to each input genome for KEGG pathway visualization.


## Circos plot visualizations


**BTEXgenie hits:**
```bash
btex-run-circos \
  --hmmscan /path/to/btex_hmm_summary.csv \
  -g test_genomes/dna_fastas/Aromatoleum_bremense_PbN1T.fna \
  -o /path/to/output_dir \
  -s "Aromatoleum_bremense_PbN1T"
```

**KOfam and BTEXgenie hits:**

```bash
  btex-run-circos \
  --hmmscan /path/to/btex_hmm_summary.csv \
  -g test_genomes/dna_fastas/Aromatoleum_bremense_PbN1T.fna \
  -o /path/to/outdir \
  -s sample_name \
  --prodigal-gbk /path/to/sample_prodigal.gbk \
  --kofam-output /path/to/kofam_abv_thres.tsv
```

> [!Note]
> Optionally, `btex-run-circos` takes `--prodigal-gbk`, which specifies the Prodigal GenBank file used to parse genomic coordinates of genes, and `--kofam-output`, which contains all hits to the KOfam HMM database. These files are produced by `btex-annotate` with the `--kofam` flag with DNA FASTA inputs.


**Input:**
- `btex-run-circos` takes `btex_hmm_summary.csv` together with the genome sequence file for a single sample. In the example above, the genome sequence file for Aromatoleum bremense PbN1T in the test_genomes folder is used.
- `--window-size` can be use to adjust the window size used to calculate GC-skew for better visualization (default: 5000 bp).
- Provide a sample name with `-s` that exactly matches the prefix of the corresponding {genome}.faa file.

**Output**

1. `{output_dir}/circos_plot.pdf`
   Main visualization showing the genome track, GC skew track, and genomic distribution of BTEXgenie hits. Optionally includes a KOfam track displaying hits to xenobiotic degradation pathways on KEGG.

2. `{output_dir}/kofam_density_track_windows.tsv`
   Table of pathway density values across genomic windows for xenobiotic degradation pathways.

3. `{output_dir}/btex_genie_hits.gbk`
   GenBank formatted file listing genes identified as BTEXgenie hits.

4. `{output_dir}/kofam_category_hits.tsv`
   Table of all KOfam hits with their corresponding KO identifiers.

5. `gene_hits.tsv`, `karyotype.tsv`, `hmm_colors.tsv`, `contig_length.tsv`
   Configuration files used to generate the Circos plot.


**Example output using the Aromatoleum bremense PbN1T genome:**
<p align="center">
  <img src="img/Aromatoleum_bremense_PbN1T_circos.png" width="800">
</p>
