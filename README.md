# LargeTarget

Scripts for the Large Neurogenomic Target Hypothesis project.

See:
Stanley CE Jr, Kulathinal RJ. Neurogenomics and the role of a large mutational target on rapid behavioral change. Biol Direct. 2016 Nov 8;11(1):60. doi: 10.1186/s13062-016-0162-1. PMID: 27825385; PMCID: PMC5101817.

Directories:

- **Dmel_genomic_coverage**: Sciprts for exploring the hypothesis on *Drosophila melanogaster* data
    - **Dmel_DGN_vcf_maker**: making SNP vcfs from Drosophila Genome Nexus Consensus Sequence data
    - **Dmel_genomic_architecture_site_counting**: Genomic architecture of *D. melanogaster*
    - **Dmel_pi_Fst_per_gene**: Gene-level pi for *D. melanogaster*
    - **Dmel_snpEff_snpSift**: Comparing synonymous and nonsynonymous sites per GO class

- **Homo_gene_expression_genomic_architecture**: Scripts for exploring the hyupothesis on *Homo sapiens* data
    - **Gene_Ontology_genomic_coverage**: Genomic architecture of H. sapiens using Gene ontology
    - **gene_expression_genomic_coverage**: testing the hypothesis on gene expression and tissue specificity data

- **fitness4gene_size**: Python tool for simulating evolution on genes of different CDS length and TFBS length and number
    - see README.md in directory for documentation and usuage instructions of the tool

