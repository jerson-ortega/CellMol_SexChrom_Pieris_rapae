# Sex-Chromosome Dosage in a Selected Animal Nuclear Genome Using NCBI, Galaxy, and R

**Activity:** Sex-Chromosome Dosage Analysis (`CellMol_SexChrom_Pieris_rapae`)  
**Name:** Jerson Lloyd T. Ortega  
**Program:** BS Biology (Cell and Molecular Biology)  

---

## Selected Species & Study Question

* **Species:** *Pieris rapae* (Small white butterfly)
* **Sex-Determination System:** Female-heterogametic **ZZ/ZW** diploid system ($2n = 50$; 24 autosome pairs + Z/W)
* **Study Question:** Can whole-genome sequencing (WGS) read depth normalized against autosomal median coverage reliably infer the chromosomal sex of *Pieris rapae* (`ZW` female vs. `ZZ` male)?

---

## Accession IDs & Workflow Links

* **Reference Genome Assembly Accession (RefSeq):** `GCF_905147795.1` (`ilPieRapa1.1`)
* **Assembly BioSample Accession:** Darwin Tree of Life reference specimen (distinct from WGS donor; verify exact `SAMEA` ID on the `GCF_905147795.1` NCBI Assembly header)
* **WGS Run Accession (SRA):** `ERR6002744` (Illumina paired-end sequencing; 43,214,315 read pairs, 151 bp, 13.1 Gbp total)
* **WGS BioSample Accession:** `SAMEA7523164` (Recorded sex: **Female**)
* **Galaxy Shared History (`Sex_Dosage_Pieris_rapae`):** [Paste Your Shareable Galaxy History URL Here]

---

## Chromosome Groups & Selection (`chromosome_groups.csv`)

* **Autosomes (`A`):** 24 RefSeq nuclear autosomes (`NC_059509.1` through `NC_059532.1`), used to compute the diploid median autosomal sequencing depth baseline.
* **W Chromosome (`W`):** RefSeq accession `NC_059533.1`, representing the female-limited sex chromosome (`W:A ≈ 0.50` in females, `0.00` in males).
* **Z Chromosome (`Z`):** RefSeq accession `NC_059534.1`, representing the shared sex chromosome (`Z:A ≈ 0.50` in hemizygous `ZW` females, `1.00` in homogametic `ZZ` males).
* **Exclusions:** The mitochondrial genome (`NC_060081.1`) and unplaced scaffolds (`NW_*`) were excluded from `chromosome_groups.csv` because organelle copy number and repetitive unplaced contigs distort nuclear diploid depth normalization.

---

## Bioinformatic Workflow & Quantitative Results

1. **Quality Control & Read Mapping (Galaxy):** Raw paired-end reads (`ERR6002744`) were retrieved via `fasterq-dump` and evaluated with **FastQC** (all read positions exceeded Phred Q30 with zero poor-quality sequences and no adapter contamination, requiring no trimming). Reads were aligned to `GCF_905147795.1` using **Map with BWA-MEM**.
2. **Coverage Calculation (Galaxy):** Per-chromosome mean read depth was computed from the coordinate-sorted BAM alignment using **mosdepth** with a minimum mapping quality filter of `MAPQ >= 20` (`chromosome_coverage.tsv`).
3. **Dosage Normalization & Plotting (Base R):** Relative read depths (`meandepth / median(A)`) were calculated and plotted using `sex_depth_analysis.R`.

| Metric | Observed Value | Expected Female (`ZW`) | Expected Male (`ZZ`) |
| :--- | :---: | :---: | :---: |
| **Median Autosome Depth (`A`)** | **`40.28×`** | `1.00×` baseline | `1.00×` baseline |
| **`Z / A` Relative Depth Ratio** | **`0.539`** | `~0.50` (`<= 0.65`) | `~1.00` (`>= 0.80`) |
| **`W / A` Relative Depth Ratio** | **`0.509`** | `~0.50` | `~0.00` |
| **Dosage Interpretation** | **`ZW-like` (Female)** | `ZW-like` | `ZZ-like` |

![Autosomal and sex-chromosome depth](results/chromosome_depth.png)

*(Step-by-step Galaxy execution screenshots are archived in the `screenshots/` directory.)*

---

## Questions for the Report

1. **What species did you choose, and what evidence supports a diploid XX/XY or ZZ/ZW system?**
   * **Answer:** I selected *Pieris rapae* (the small white butterfly). Cytogenetic karyotyping ($2n = 50$) and chromosome-level genome assemblies (`ilPieRapa1.1` / `GCF_905147795.1`) demonstrate that lepidopterans use a female-heterogametic **ZZ/ZW** sex-chromosome system, in which males carry two Z chromosomes (`ZZ`) and females carry one Z and one W chromosome (`ZW`).

2. **What are your assembly accession, WGS run, and both BioSample accessions? Are the donors the same?**
   * **Answer:** The RefSeq reference assembly is **`GCF_905147795.1`** (`ilPieRapa1.1`), and the analyzed WGS dataset is SRA run **`ERR6002744`** from BioSample **`SAMEA7523164`**. Within the Darwin Tree of Life *Pieris rapae* project, the primary reference assembly (`GCF_905147795.1`) and the short-read WGS run (`ERR6002744` / `SAMEA7523164`) were generated from separate collected individuals rather than the exact same donor organism.

3. **Which sequences were used as autosomes and as X or Z? How did you verify their identities?**
   * **Answer:** Twenty-four RefSeq chromosomes (**`NC_059509.1` through `NC_059532.1`**) were labeled as autosomes (`A`), **`NC_059534.1`** was labeled as the **`Z`** chromosome, and **`NC_059533.1`** was labeled as the **`W`** chromosome. Their identities were verified directly from the NCBI RefSeq Assembly sequence report and FASTA headers for `GCF_905147795.1`, while mitochondrial DNA and unplaced `NW_*` scaffolds were excluded.

4. **What are the median autosomal depth, X:A or Z:A ratio, and inferred sex-chromosome complement?**
   * **Answer:** The median autosomal read depth is **`40.28×`**, the **`Z / A` ratio is `0.539`**, and the **`W / A` ratio is `0.509`**. Because the `Z / A` ratio (`0.539`) falls well below the `0.65` single-copy threshold and `W / A` (`0.509`) confirms single-copy W presence, the inferred sex-chromosome complement is **`ZW-like` (Female)**.

5. **Does your result agree with the recorded sex of the WGS donor, if reported? Explain any discrepancy or uncertainty.**
   * **Answer:** Yes. The **`ZW-like`** dosage inference agrees directly with the NCBI BioSample record for `SAMEA7523164`, which lists the donor as **female**. The slight elevation of `Z / A` (`0.539`) above the ideal `0.50` expectation reflects normal biological and technical variation, including pseudoautosomal or shared homologous sequence regions between Z and W, GC-content coverage bias, and repetitive sequence filtering at `MAPQ >= 20`.

6. **Why is a chromosome-level assembly alone insufficient for a reliable read-depth sex inference?**
   * **Answer:** A chromosome-level reference assembly (`GCF_905147795.1`) only provides a static haploid template of chromosome coordinates and carries no quantitative copy-number information for a query specimen. Reliable sex inference requires aligning actual whole-genome sequencing reads from the target individual to that reference so that relative read depths across Z, W, and autosomes can be empirically measured; additionally, repetitive heterochromatin on W chromosomes or incomplete W assembly can bias mapping depth without quality filtering and autosomal normalization.

---

## Repository Structure

```text
CellMol_SexChrom_Pieris_rapae/
├── README.md
├── chromosome_coverage.tsv
├── chromosome_groups.csv
├── sex_depth_analysis.R
├── results/
│   ├── chromosome_depth.png
│   ├── chromosome_relative_depth.csv
│   └── sex_inference.txt
└── screenshots/
