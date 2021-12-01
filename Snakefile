#!/usr/bin/env snakemake

configfile: "config.yaml"

import os
import pandas as pd

# make sure the tmp directory exists
os.makedirs(config["tmpdir"], exist_ok=True)

org2genome = pd.read_csv(config["genome_map"], sep="\t", index_col="organism").genome

rule all:
    input:
        expand("{species}/{annotation}.txcutr.w{width}.kdx", species=["homo_sapiens"],
               annotation=["Homo_sapiens.GRCh38.104",
                           "Homo_sapiens.GRCh38.93.cr_3_0_0",
                           "gencode.v38.annotation",
                           "gencode.v38.annotation.pc"], width=[500]),
        expand("{species}/{annotation}.txcutr.w{width}.kdx", species=["mus_musculus"],
               annotation=["Mus_musculus.GRCm38.93.cr_3_0_0",
                           "gencode.vM25.annotation"], width=[500])


################################################################################
## Downloading Rules
################################################################################

rule ensembl_hg38:
    output:
        gtf="homo_sapiens/Homo_sapiens.GRCh38.{release}.gtf.gz"
    params:
        url_base=lambda wcs: "ftp://ftp.ensembl.org/pub/release-%s/gtf/" % wcs.release
    wildcard_constraints:
        release='\d+'
    conda: "envs/samtools.yaml"
    shell:
        """
        wget -O {output.gtf} '{params.url_base}{output.gtf}'
        """

rule ensembl_mm10:
    output:
        gtf="mus_musculus/Mus_musculus.GRCm38.{release}.gtf.gz"
    params:
        url_base=lambda wcs: "ftp://ftp.ensembl.org/pub/release-%s/gtf/" % wcs.release
    wildcard_constraints:
        release='\d+'
    conda: "envs/samtools.yaml"
    shell:
        """
        wget -O {output.gtf} '{params.url_base}{output.gtf}'
        """

rule gencode_hg38:
    output:
        gtf="homo_sapiens/gencode.v{release}.annotation.gtf.gz"
    params:
        url_base=lambda wcs: "ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_%s/" % wcs.release
    wildcard_constraints:
        release='\d+'
    conda: "envs/samtools.yaml"
    shell:
        """
        file=$(basename {output.gtf})
        wget -O {output.gtf} \"{params.url_base}$file\"
        """

rule gencode_mm10:
    output:
        gtf="mus_musculus/gencode.v{release}.annotation.gtf.gz"
    params:
        url_base=lambda wcs: "ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_%s/" % wcs.release
    wildcard_constraints:
        release='M\d+'
    conda: "envs/samtools.yaml"
    shell:
        """
        file=$(basename {output.gtf})
        wget -O {output.gtf} \"{params.url_base}$file\"
        """
        
################################################################################
## Preprocessing Rules
################################################################################

rule clean_gtf_mRNA_ends:
    input:
        gtf="{species}/{annotation}.gtf.gz"
    output:
        gtf="{species}/{annotation}.mRNA_ends_found.gtf.gz"
    conda: "envs/samtools.yaml"
    shell:
        """
        gzip -cd {input.gtf} |
        awk '$0 !~ /mRNA_end_NF/' |
        bgzip -c > {output.gtf}
        """

rule filter_protein_coding:
    input:
        gtf="{species}/{annotation}.gtf.gz",
        awk="scripts/filter_protein_coding.awk"
    output:
        gtf="{species}/{annotation}.pc.gtf.gz"
    conda: "envs/samtools.yaml"
    shell:
        """
        gzip -cd {input.gtf} |
        awk -f {input.awk} |
        bgzip -c > {output.gtf}
        """

rule filter_cellranger_3_0_0:
    input:
        gtf="{species}/{annotation}.gtf.gz"
    output:
        gtf="{species}/{annotation}.cr_3_0_0.gtf.gz"
    # wildcard_constraints:
    #     annotation="^(?!.*\b\.(filtered_cr|mRNA_ends_found)[^.]*\b).*"
    params:
        cellranger=config["cellranger"],
        tmp=config["tmpdir"]
    conda: "envs/samtools.yaml"
    shell:
        """
        gtf_in=$(mktemp {params.tmp}/snakemake.filter_cr.XXXXXX.in.gtf)
        gtf_out=$(mktemp {params.tmp}/snakemake.filter_cr.XXXXXX.out.gtf)
        gzip -cd {input.gtf} > $gtf_in
        {params.cellranger} mkgtf $gtf_in $gtf_out \
            --attribute=gene_biotype:protein_coding \
            --attribute=gene_biotype:lincRNA \
            --attribute=gene_biotype:antisense \
            --attribute=gene_biotype:IG_LV_gene \
            --attribute=gene_biotype:IG_V_gene \
            --attribute=gene_biotype:IG_V_pseudogene \
            --attribute=gene_biotype:IG_D_gene \
            --attribute=gene_biotype:IG_J_gene \
            --attribute=gene_biotype:IG_J_pseudogene \
            --attribute=gene_biotype:IG_C_gene \
            --attribute=gene_biotype:IG_C_pseudogene \
            --attribute=gene_biotype:TR_V_gene \
            --attribute=gene_biotype:TR_V_pseudogene \
            --attribute=gene_biotype:TR_D_gene \
            --attribute=gene_biotype:TR_J_gene \
            --attribute=gene_biotype:TR_J_pseudogene \
            --attribute=gene_biotype:TR_C_gene 
        bgzip -c $gtf_out > {output.gtf}
        rm -f $gtf_in $gtf_out
        """

################################################################################
## txcutr Steps
################################################################################

rule txcutr_ensembl_gtf:
    input:
        gtf="{species}/{annotation}.mRNA_ends_found.gtf.gz"
    output:
        gtf="{species}/{annotation}.txcutr.w{width}.gtf.gz",
        fa="{species}/{annotation}.txcutr.w{width}.fa.gz",
        tsv="{species}/{annotation}.txcutr.w{width}.merge.tsv"
    params:
        mergeDist=200,
        genome=lambda wcs: org2genome[wcs.species]
    wildcard_constraints:
        width="\d+"
    threads: 20
    resources:
        mem_mb=4000
    conda: "envs/txcutr.yaml"
    script: "scripts/txcutr_gtf.R"

################################################################################
## Indexing Steps
################################################################################

rule kallisto_index_txcutr:
    input:
        "{species}/{annotation}.txcutr.w{width}.fa.gz"
    output:
        "{species}/{annotation}.txcutr.w{width}.kdx"
    conda: "envs/kallisto.yaml"
    resources:
        mem_mb=16000
    shell:
        """
        kallisto index -i {output} {input}
        """
