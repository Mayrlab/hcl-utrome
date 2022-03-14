#!/usr/bin env snakemake --snakefile

configfile: "config.yaml"

import pandas as pd
from sys import stderr
import os

# print to stderr
def message(*args, **kwargs):
    print(*args, file=stderr, **kwargs)
    
# ensure tmpdir exists
os.makedirs(config['tmpdir'], exist_ok=True)


message("[INFO] Loading metadata...")
metadata = pd.read_csv(config['metadataFile'], usecols=["sample_id", "Run"], index_col="sample_id")
message("[INFO] Found %d SRA runs." % len(metadata.index))

message("[INFO] Loading samples list...")
with open(config['samplesFile']) as f:
    samples_subset = f.read().splitlines()
message("[INFO] Found %d sample IDs." % len(samples_subset))

message("[INFO] Loading celltypes-samples list...")
celltype_sample_map = pd.read_csv(config['celltypesSamplesFile'])
message("[INFO] Found %d celltypes." % celltype_sample_map.celltype_id.nunique())
message("[INFO] Found %d celltype-sample pairs." % len(celltype_sample_map))


rule all:
    input:
        expand("data/bam/samples/{sample_id}.tagged.bam",
               sample_id=samples_subset),
        expand("data/bam/celltypes/{celltype_id}.bam",
               celltype_id=list(celltype_sample_map.celltype_id.unique()))
        # expand("data/kdx/utrome.e{epsilon}.t{threshold}.gc{version}.pas{tpm}.f{likelihood}.w{width}.kdx",
        #        epsilon=[3,5,10], threshold=[10,100,200,1000], version=[39], tpm=[3], likelihood=[0.9999], width=[500]),
        # expand("data/gff/utrome.e{epsilon}.t{threshold}.gc{version}.pas{tpm}.f{likelihood}.w{width}.m{merge}.tsv",
        #        epsilon=[3,5,10], threshold=[10,100,200,1000], version=[39], tpm=[3], likelihood=[0.9999], width=[500], merge=[200])

rule download_fastq:
    output:
        r1=temp("data/fastq/raw/{sample_id}_R1.fastq.gz"),
        r2=temp("data/fastq/raw/{sample_id}_R2.fastq.gz")
    params:
        srr=lambda wcs: metadata.Run[wcs.sample_id],
        tmpdir=config['tmpdir']
    conda: "envs/sratools.yaml"
    threads: 8
    shell:
        """
        ## expected output files
        tmp_r1="{params.tmpdir}/{params.srr}_1.fastq"
        tmp_r2="{params.tmpdir}/{params.srr}_2.fastq"
        
        ## dump raw fastq
        fasterq-dump -e {threads} -O {params.tmpdir} -t {params.tmpdir} -sS {params.srr}
        
        ## compress
        bgzip -@ {threads} -c $tmp_r1 > {output.r1}
        bgzip -@ {threads} -c $tmp_r2 > {output.r2}

        ## clean up
        rm $tmp_r1 $tmp_r2
        """

ruleorder: merge_subsamples > pear_merge
rule pear_merge:
    input:
        r1="data/fastq/raw/{sample_id}_R1.fastq.gz",
        r2="data/fastq/raw/{sample_id}_R2.fastq.gz"
    output:
        r12="data/fastq/assembled/{sample_id}.assembled.fastq.gz"
    params:
        tmpdir=config['tmpdir'] + "/pear",
        min_length=54+21,
        max_pval=0.0001
    conda: "envs/pear.yaml"
    threads: 12
    resources:
        mem_mb=1000
    shell:
        """
        mkdir -p {params.tmpdir}
        pear -j {threads} \\
          -n {params.min_length} -p {params.max_pval} \\
          -f {input.r1} -r {input.r2} \\
          -o {params.tmpdir}/{wildcards.sample_id}
        bgzip -@ {threads} -c {params.tmpdir}/{wildcards.sample_id}.assembled.fastq > {output.r12}
        rm {params.tmpdir}/{wildcards.sample_id}.*.fastq
        """

## NB: This method only applies to particular samples where the cell-type annotations
## do not directy correspond with a single sample_id as deposited on the SRA. To be
## processed by this rule, there must exist a corresponding file under:
##   metadata/merge/merge-{sample_id}.txt
def get_merge_targets (wcs):
    path = "metadata/merge/merge-%s.txt" % wcs.sample_id
    if os.path.isfile(path):
        with open(path) as f:
            sample_ids = f.read().splitlines()
        return expand("data/fastq/assembled/{sample_id}.assembled.fastq.gz",
                      sample_id=sample_ids)
    else:
        return ""

rule merge_subsamples:
    input:
        merge="metadata/merge/merge-{sample_id}.txt",
        fqs=get_merge_targets
    output:
        fq="data/fastq/assembled/{sample_id}.assembled.fastq.gz"
    shell:
        """
        cat {input.fqs} > {output.fq}
        """

rule extract_whitelist:
    input:
        tsv=config['annotsFile']
    output:
        txt="data/barcodes/{sample_id}.whitelist.txt"
    shell:
        """
        gzip -cd {input.tsv} |\\
          awk '{{ if ($2 ~ /^{wildcards.sample_id}$/) print $3 }}' > {output.txt}
        """

rule umitools_extract_assembled:
    input:
        bx="data/barcodes/{sample_id}.whitelist.txt",
        fq="data/fastq/assembled/{sample_id}.assembled.fastq.gz"
    output:
        temp("data/fastq/extracted/{sample_id}.assembled.bx.fastq.gz")
    resources:
        mem_mb=16000
    conda: "envs/umitools.yaml"
    shell:
        """
        umi_tools extract \\
          --filter-cell-barcode \\
          --extract-method=regex \\
          --bc-pattern='(?P<cell_1>.{{6}})(?P<discard_1>CGACTCACTACAGGG){{s<=1}}(?P<cell_2>.{{6}})(?P<discard_2>TCGGTGACACGATCG){{s<=1}}(?P<cell_3>.{{6}})(?P<umi_1>.{{6}})(T{{12}}){{s<=2}}.*' \\
          --whitelist={input.bx} \\
          --stdin={input.fq} \\
          --stdout={output}
        """

rule cutadapt_polyT_SE:
    input:
        "data/fastq/extracted/{sample_id}.assembled.bx.fastq.gz"
    output:
        "data/fastq/trimmed/{sample_id}.assembled.clean.fastq.gz"
    params:
        min_length=config['minReadLength']
    conda: "envs/cutadapt.yaml"
    threads: 8
    resources:
        mem_mb=1000
    shell:
        """
        cutadapt --cores={threads} --front='T{{100}}' \\
        --minimum-length={params.min_length} --length-tag='length=' \\
        --output={output} {input}
        """

rule hisat2_SE:
    input:
        "data/fastq/trimmed/{sample_id}.assembled.clean.fastq.gz"
    output:
        bam="data/bam/samples/{sample_id}.assembled.bam",
        log="qc/hisat2/{sample_id}.assembled.log"
    params:
        tmpdir=config['tmpdir'],
        idx=config['hisatIndex'],
        sam=config['tmpdir'] + "/{sample_id}.assembled.sam"
    conda: "envs/hisat2.yaml"
    threads: 16
    resources:
        mem_mb=1000
    shell:
        """
        hisat2 -p {threads} -x {params.idx} \\
          -U {input} -S {params.sam} \\
          --rna-strandness R \\
          --new-summary --summary-file {output.log}
        samtools sort -@ {threads} -T {params.tmpdir}/ -o {output.bam} {params.sam}
        rm -f {params.sam}
        """

rule tag_bx_umi:
    input:
        bam="data/bam/samples/{sample_id}.assembled.bam",
        script="scripts/tag_bx_umi.awk"
    output:
        bam="data/bam/samples/{sample_id}.tagged.bam",
        bai="data/bam/samples/{sample_id}.tagged.bam.bai"
    conda: "envs/hisat2.yaml"
    shell:
        """
        samtools view -h {input.bam} |\\
          awk -f {input.script} |\\
          samtools view -b > {output.bam}
        samtools index {output.bam}
        """

rule extract_celltype_sample_bxs:
    input:
        tsv=config['annotsFile']
    output:
        bxs="data/barcodes/celltypes/{cluster_id}-{celltype}/{cluster_id}-{celltype}.{sample_id}.bxs.txt"
    wildcard_constraints:
        cluster_id="\d+",
        celltype="[^.]+",
        sample_id="[^.]+(\.\d+)?"
    shell:
        """
        gzip -cd {input.tsv} |\\
          awk '{{ if ($2 ~ /^{wildcards.sample_id}$/ && $7 ~ /^{wildcards.cluster_id}$/) print $3 }}' > {output.bxs}
        """

rule filter_celltype_sample:
    input:
        bam="data/bam/samples/{sample_id}.tagged.bam",
        bxs="data/barcodes/celltypes/{cluster_id}-{celltype}/{cluster_id}-{celltype}.{sample_id}.bxs.txt"
    output:
        bam=temp("data/bam/celltypes/{cluster_id}-{celltype}/{cluster_id}-{celltype}.{sample_id}.bam")
    wildcard_constraints:
        cluster_id="\d+",
        celltype="[^.]+",
        sample_id="[^.]+(\.\d+)?"
    conda: "envs/hisat2.yaml"
    shell:
        """
        samtools view -D CB:{input.bxs} -o {output.bam} {input.bam}
        """

def get_celltype_samples (wcs):
    return expand("data/bam/celltypes/{celltype_id}/{celltype_id}.{sample_id}.bam",
                  celltype_id=wcs.celltype_id,
                  sample_id=celltype_sample_map.sample_id[celltype_sample_map.celltype_id == wcs.celltype_id])
    
rule merge_celltype_samples:
    input:
        bams=get_celltype_samples
    output:
        bam="data/bam/celltypes/{celltype_id}.bam",
        bai="data/bam/celltypes/{celltype_id}.bam.bai"
    wildcard_constraints:
        celltype_id="\d+-[^.]+"
    conda: "envs/hisat2.yaml"
    shell:
        """
        samtools merge -o {output.bam} {input.bams}
        samtools index {output.bam}
        """

rule cleavage_coverage:
    input:
        "data/bam/celltypes/{sample_id}.assembled.bam"
    output:
        neg="data/coverage/{sample_id}.assembled.negative.txt.gz",
        pos="data/coverage/{sample_id}.assembled.positive.txt.gz"
    conda: "envs/bedtools.yaml"
    shell:
        """
        bedtools genomecov -dz -5 -strand '-' -ibam {input} | gzip > {output.neg}
        bedtools genomecov -dz -5 -strand '+' -ibam {input} | gzip > {output.pos}
        """

rule sum_coverage:
    input:
        neg=expand("data/coverage/{sample_id}.assembled.negative.txt.gz",
                     sample_id=samples_subset),
        pos=expand("data/coverage/{sample_id}.assembled.positive.txt.gz",
                     sample_id=samples_subset)
    output:
        neg="data/coverage/utrome.assembled.negative.txt.gz",
        pos="data/coverage/utrome.assembled.positive.txt.gz"
    conda: "envs/bedtools.yaml"
    threads: 4
    resources:
        mem_mb=1000
    shell:
        """
        gzip -cd {input.neg} |\\
          datamash -s -g 1,2 sum 3 |\\
          sort -k 1,1 -k 2,2n |\\
          bgzip -@ {threads} > {output.neg}

        gzip -cd {input.pos} |\\
          datamash -s -g 1,2 sum 3 |\\
          sort -k 1,1 -k 2,2n |\\
          bgzip -@ {threads} > {output.pos}
        """

rule merge_coverage:
    input:
        cov="data/coverage/utrome.assembled.{strand}.txt.gz"
    output:
        cov="data/coverage/utrome.assembled.{strand}.e{epsilon,\d+}.txt.gz"
    threads: 16
    resources:
        mem_mb=1000
    conda: "envs/pandas.yaml"
    script:
        "scripts/merge_genomecov.py"

rule cov_to_bed:
    input:
        neg="data/coverage/utrome.assembled.negative.e{epsilon}.txt.gz",
        pos="data/coverage/utrome.assembled.positive.e{epsilon}.txt.gz"
    output:
        bed="data/bed/cleavage-sites/utrome.cleavage.e{epsilon}.t{threshold}.bed.gz",
        tbi="data/bed/cleavage-sites/utrome.cleavage.e{epsilon}.t{threshold}.bed.gz.tbi"
    params:
        tmpdir=config['tmpdir']
    wildcard_constraints:
        epsilon="\d+",
        threshold="\d+"
    conda: "envs/bedtools.yaml"
    shell:
        """
        tmpbed=$(mktemp -p {params.tmpdir})
        bgzip -cd {input.neg} |\\
          awk -v OFS='\\t' '$3>={wildcards.threshold}{{print $1, $2, $2, $1 \":\" $2 \":+\", $3, \"+\"}}' > $tmpbed
        bgzip -cd {input.pos} |\\
          awk -v OFS='\\t' '$3>={wildcards.threshold}{{print $1, $2, $2, $1 \":\" $2 \":-\", $3, \"-\"}}' >> $tmpbed
        sort -k1,1 -k2,2n $tmpbed | bgzip > {output.bed}
        tabix -0 -p bed {output.bed}
        rm -f $tmpbed
        """

rule cleanUpdTSeq_classify:
    input:
        bed="data/bed/cleavage-sites/utrome.cleavage.e{epsilon}.t{threshold}.bed.gz"
    output:
        probs="data/cleavage-sites/utrome.classification.e{epsilon}.t{threshold}.tsv.gz"
    wildcard_constraints:
        epsilon="\d+",
        threshold="\d+"
    conda: "envs/bioc_3_14.yaml"
    resources:
        mem_mb=4000
    threads: 24
    script:
        "scripts/classify_cleavage_sites.R"

rule cleanUpdTSeq_filter:
    input:
        bed="data/bed/cleavage-sites/utrome.cleavage.e{epsilon}.t{threshold}.bed.gz",
        probs="data/cleavage-sites/utrome.classification.e{epsilon}.t{threshold}.tsv.gz",
    output:
        bed="data/bed/cleavage-sites/utrome.cleavage.e{epsilon}.t{threshold}.f{likelihood}.bed.gz"
    wildcard_constraints:
        epsilon="\d+",
        threshold="\d+",
        likelihood="0.\d+"
    conda: "envs/bioc_3_14.yaml"
    script:
        "scripts/filter_cleavage_sites.R"

rule download_polyASite_atlas:
    output:
        atlas="data/cleavage-sites/polyAsite.atlas.tsv.gz"
    params:
        url="https://polyasite.unibas.ch/download/atlas/2.0/GRCh38.96/atlas.clusters.2.0.GRCh38.96.tsv.gz"
    shell:
        """
        wget -O {output.atlas} '{params.url}'
        """

rule download_gencode_gff:
    output:
        gff="data/gff/gencode.v{version}.annotation.gff3.gz"
    params:
        url_base="ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/",
        url_file=lambda wcs: "release_%s/gencode.v%s.annotation.gff3.gz" % (wcs.version, wcs.version)
    wildcard_constraints:
        version='\d+'
    shell:
        """
        wget -O {output.gff} '{params.url_base}{params.url_file}'
        """

rule filter_gencode_mRNA_ends:
    input:
        gff="data/gff/gencode.v{version}.annotation.gff3.gz"
    output:
        gff="data/gff/gencode.v{version}.mRNA_ends_found.gff3.gz",
        tbi="data/gff/gencode.v{version}.mRNA_ends_found.gff3.gz.tbi"
    wildcard_constraints:
        version='\d+'
    conda: "envs/bedtools.yaml"
    threads: 4
    shell:
        """
        bgzip -cd {input.gff} |\\
          awk '$0 !~ /^#/' |\\
          awk '$0 !~ /mRNA_end_NF/' |\\
          sort --parallel={threads} -S4G -k1,1 -k4,4n |\\
          bgzip -@ {threads} -c > {output.gff}

        tabix {output.gff}
        """
        
rule filter_validated_sites:
    input:
        bed_all="data/bed/cleavage-sites/utrome.cleavage.e{epsilon}.t{threshold}.bed.gz",
        gencode="data/gff/gencode.v{version}.mRNA_ends_found.gff3.gz"
    output:
        validated="data/bed/cleavage-sites/utrome.validated.e{epsilon}.t{threshold}.gc{version}.bed.gz",
        unvalidated="data/bed/cleavage-sites/utrome.unvalidated.e{epsilon}.t{threshold}.gc{version}.bed.gz"
    wildcard_constraints:
        epsilon="\d+",
        threshold="\d+",
        version="\d+"
    params:
        radius=config['radiusGENCODE']
    conda: "envs/bioc_3_14.yaml"
    script:
        "scripts/filter_validated_sites.R"

rule filter_supported_sites:
    input:
        unvalidated="data/bed/cleavage-sites/utrome.unvalidated.e{epsilon}.t{threshold}.gc{version}.bed.gz",
        atlas="data/cleavage-sites/polyAsite.atlas.tsv.gz"
    output:
        supported="data/bed/cleavage-sites/utrome.supported.e{epsilon}.t{threshold}.gc{version}.pas{tpm}.bed.gz",
        unsupported="data/bed/cleavage-sites/utrome.unsupported.e{epsilon}.t{threshold}.gc{version}.pas{tpm}.bed.gz"
    wildcard_constraints:
        epsilon="\d+",
        threshold="\d+",
        version="\d+",
        tpm="\d+"
    params:
        radius=config['radiusPAS']
    conda: "envs/bioc_3_14.yaml"
    script:
        "scripts/filter_supported_sites.R"

rule filter_likely_sites:
    input:
        unsupported="data/bed/cleavage-sites/utrome.unsupported.e{epsilon}.t{threshold}.gc{version}.pas{tpm}.bed.gz",
        passing="data/bed/cleavage-sites/utrome.cleavage.e{epsilon}.t{threshold}.f{likelihood}.bed.gz"
    output:
        likely="data/bed/cleavage-sites/utrome.likely.e{epsilon}.t{threshold}.gc{version}.pas{tpm}.f{likelihood}.bed.gz",
        unlikely="data/bed/cleavage-sites/utrome.unlikely.e{epsilon}.t{threshold}.gc{version}.pas{tpm}.f{likelihood}.bed.gz"
    wildcard_constraints:
        epsilon="\d+",
        threshold="\d+",
        version="\d+",
        tpm="\d+",
        likelihood="0.\d+"
    conda: "envs/bioc_3_14.yaml"
    script:
        "scripts/filter_likely_sites.R"

rule export_cleavage_sites:
    input:
        supported="data/bed/cleavage-sites/utrome.supported.e{epsilon}.t{threshold}.gc{version}.pas{tpm}.bed.gz",
        likely="data/bed/cleavage-sites/utrome.likely.e{epsilon}.t{threshold}.gc{version}.pas{tpm}.f{likelihood}.bed.gz",
        gencode="data/gff/gencode.v{version}.mRNA_ends_found.gff3.gz"
    output:
        utr3="data/bed/cleavage-sites/utrome.utr3.e{epsilon}.t{threshold}.gc{version}.pas{tpm}.f{likelihood}.bed.gz",
        extutr3="data/bed/cleavage-sites/utrome.extutr3.e{epsilon}.t{threshold}.gc{version}.pas{tpm}.f{likelihood}.bed.gz"
    wildcard_constraints:
        epsilon="\d+",
        threshold="\d+",
        version="\d+",
        tpm="\d+",
        likelihood="0.\d+"
    params:
        ext_utr3=config['extUTR3'],
        ext_utr5=config['extUTR5']
    conda: "envs/bioc_3_14.yaml"
    script:
        "scripts/export_cleavage_sites.R"

rule augment_transcriptome_utr3:
    input:
        utr3="data/bed/cleavage-sites/utrome.utr3.e{epsilon}.t{threshold}.gc{version}.pas{tpm}.f{likelihood}.bed.gz",
        gencode="data/gff/gencode.v{version}.mRNA_ends_found.gff3.gz"
    output:
        gff="data/gff/txs.utr3.e{epsilon}.t{threshold}.gc{version}.pas{tpm}.f{likelihood}.gff3.gz",
        gtf="data/gff/txs.utr3.e{epsilon}.t{threshold}.gc{version}.pas{tpm}.f{likelihood}.gtf.gz"
    wildcard_constraints:
        epsilon = "\d+",
        threshold = "\d+",
        version="\d+",
        tpm="\d+",
        likelihood = "0.\d+"
    conda: "envs/bioc_3_14.yaml"
    resources:
        mem_mb=16000
    script: "scripts/augment_transcriptome_utr3.R"

rule augment_transcriptome_extutr3:
    input:
        extutr3="data/bed/cleavage-sites/utrome.extutr3.e{epsilon}.t{threshold}.gc{version}.pas{tpm}.f{likelihood}.bed.gz",
        gencode="data/gff/gencode.v{version}.mRNA_ends_found.gff3.gz"
    output:
        gff="data/gff/txs.extutr3.e{epsilon}.t{threshold}.gc{version}.pas{tpm}.f{likelihood}.gff3.gz",
        gtf="data/gff/txs.extutr3.e{epsilon}.t{threshold}.gc{version}.pas{tpm}.f{likelihood}.gtf.gz"
    wildcard_constraints:
        epsilon = "\d+",
        threshold = "\d+",
        version="\d+",
        tpm="\d+",
        likelihood = "0.\d+"
    params:
        ext_utr3=config['extUTR3']
    conda: "envs/bioc_3_14.yaml"
    resources:
        mem_mb=16000
    script: "scripts/augment_transcriptome_extutr3.R"

rule create_chunked_granges:
    input:
        utr3="data/gff/txs.utr3.e{epsilon}.t{threshold}.gc{version}.pas{tpm}.f{likelihood}.gff3.gz",
        extutr3="data/gff/txs.extutr3.e{epsilon}.t{threshold}.gc{version}.pas{tpm}.f{likelihood}.gff3.gz",
        gencode="data/gff/gencode.v{version}.mRNA_ends_found.gff3.gz"
    output:
        plus="data/granges/augmented.plus.chunked.e{epsilon}.t{threshold}.gc{version}.pas{tpm}.f{likelihood}.Rds",
        minus="data/granges/augmented.minus.chunked.e{epsilon}.t{threshold}.gc{version}.pas{tpm}.f{likelihood}.Rds",
    wildcard_constraints:
        epsilon="\d+",
        threshold="\d+",
        version="\d+",
        tpm="\d+",
        likelihood="0.\d+"
    conda: "envs/bioc_3_14.yaml"
    resources:
        mem_mb=16000
    script: "scripts/create_chunked_granges.R"

rule truncate_plus_strand:
    input:
        granges="data/granges/augmented.plus.chunked.e{epsilon}.t{threshold}.gc{version}.pas{tpm}.f{likelihood}.Rds"
    output:
        granges="data/granges/utrome.raw.plus.e{epsilon}.t{threshold}.gc{version}.pas{tpm}.f{likelihood}.w{width}.Rds"
    wildcard_constraints:
        epsilon="\d+",
        threshold="\d+",
        version="\d+",
        tpm="\d+",
        likelihood="0.\d+",
        width="\d+"
    conda: "envs/bioc_3_14.yaml"
    threads: 24
    resources:
        mem_mb=3000
    script: "scripts/truncate_plus_strand.R"

rule truncate_minus_strand:
    input:
        granges="data/granges/augmented.minus.chunked.e{epsilon}.t{threshold}.gc{version}.pas{tpm}.f{likelihood}.Rds"
    output:
        granges="data/granges/utrome.raw.minus.e{epsilon}.t{threshold}.gc{version}.pas{tpm}.f{likelihood}.w{width}.Rds"
    wildcard_constraints:
        epsilon="\d+",
        threshold="\d+",
        version="\d+",
        tpm="\d+",
        likelihood="0.\d+",
        width="\d+"
    conda: "envs/bioc_3_14.yaml"
    threads: 24
    resources:
        mem_mb=3000
    script: "scripts/truncate_minus_strand.R"

rule export_unmerged_utrome:
    input:
        plus="data/granges/utrome.raw.plus.e{epsilon}.t{threshold}.gc{version}.pas{tpm}.f{likelihood}.w{width}.Rds",
        minus="data/granges/utrome.raw.minus.e{epsilon}.t{threshold}.gc{version}.pas{tpm}.f{likelihood}.w{width}.Rds"
    output:
        gtf="data/gff/utrome.e{epsilon}.t{threshold}.gc{version}.pas{tpm}.f{likelihood}.w{width}.gtf",
        fa="data/gff/utrome.e{epsilon}.t{threshold}.gc{version}.pas{tpm}.f{likelihood}.w{width}.fasta.gz"
    wildcard_constraints:
        epsilon="\d+",
        threshold="\d+",
        version="\d+",
        tpm="\d+",
        likelihood="0.\d+",
        width="\d+"
    conda: "envs/bioc_3_14.yaml"
    threads: 24
    resources:
        mem_mb=2000
    script: "scripts/export_unmerged_utrome.R"

rule sort_utrome:
    input:
        gtf="data/gff/utrome.{settings}.gtf"
    output:
        gtf="data/gff/utrome.{settings}.gtf.gz",
        idx="data/gff/utrome.{settings}.gtf.gz.tbi"
    conda: "envs/bedtools.yaml"
    threads: 4
    resources:
        mem_mb=6000
    shell:
        """
        cat {input.gtf} |\\
          awk '!( $0 ~ /^#/ )' |\\
          sort --parallel={threads} -S4G -k1,1 -k4,4n |\\
          bgzip -@ {threads} -c > {output.gtf}

        tabix {output.gtf}
        """

rule kallisto_index:
    input:
        fa="data/gff/utrome.{settings}.fasta.gz"
    output:
        kdx="data/kdx/utrome.{settings}.kdx"
    conda: "envs/kallisto.yaml"
    resources:
        mem_mb=16000
    shell:
        """
        kallisto index -i {output.kdx} {input.fa}
        """

rule export_merge_table:
    input:
        gtf="data/gff/utrome.{settings}.gtf.gz"
    output:
        tsv="data/gff/utrome.{settings}.m{merge}.tsv"
    params:
        genome="hg38"
    conda: "envs/bioc_3_14.yaml"
    resources:
        mem_mb=8000
    script: "scripts/export_merge_table.R"

################################################################################
## Reports
################################################################################

rule summarize_pear_results:
    input:
        script="scripts/parse_pear.awk"
    output:
        csv="qc/pear/pear_summary.csv"
    shell:
        """
        echo 'sample,cts_total,cts_assembled,cts_unassembled,cts_discarded,pct_assembled,pct_unassembled,pct_discarded' > {output.csv}
        for log in logs/*/pear_merge/*/*.out; do
            awk -f {input.script} $log >> {output.csv}
        done
        """

rule summarize_umitools_results:
    input:
        script="scripts/parse_umitools_extract.awk"
    output:
        csv="qc/barcodes/umitools_summary.csv"
    shell:
        """
        echo 'sample,total_in,total_out,raw_matched,raw_unmatched,excluded_bx,corrected_bx' > {output.csv}
        for log in logs/*/umitools_extract_assembled/*/*.out; do
            awk -f {input.script} $log >> {output.csv}
        done
        """

rule summarize_hisat_results:
    input:
        script="scripts/parse_hisat2.awk"
    output:
        csv="qc/hisat2/hisat2_summary.csv"
    shell:
        """
        echo 'sample,cts_total,cts_unaligned,cts_unique,cts_multi' > {output.csv}
        for log in qc/hisat2/*.log; do
            awk -f {input.script} $log >> {output.csv}
        done
        """
