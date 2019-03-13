#######################################################
# this script is to analysis chip-seq, before please carefor the same directory of config.yaml file, change the information
# write by zhanweimin, 630950832@qq.com
# 2019.3.13
#######################################################

configfile: "config.yaml"

samples = config["samples"]
contrals = config["contrals"]
sample = samples + contrals
reps = config["reps"]


#here will defination the output file
ALL_filter_fastq=[]
for SAMPLE in sample:
    for REP in reps:
        ALL_filter_fastq.append("1_clean_data/{sample}_{rep}_1.fastq.gz".format(sample=SAMPLE,rep=REP))
        ALL_filter_fastq.append("1_clean_data/{sample}_{rep}_1_unpaired.fastq.gz".format(sample=SAMPLE,rep=REP))
        ALL_filter_fastq.append("1_clean_data/{sample}_{rep}_2.fastq.gz".format(sample=SAMPLE,rep=REP))
        ALL_filter_fastq.append("1_clean_data/{sample}_{rep}_2_unpaired.fastq.gz".format(sample=SAMPLE,rep=REP))

ALL_bam=expand("3_filter/{sample}_{rep}_sorted.bam",sample=sample,rep=reps)
ALL_bam_bai=expand("3_filter/{sample}_{rep}_sorted.bam.bai",sample=sample,rep=reps)
ALL_peak=expand("4_call_peak/{samples}_{rep}_{contrals}",samples=samples,rep=reps,contrals=contrals)
ALL_bigwig=expand("5_transcrate/{sample}_{rep}.bw",sample=sample,rep=reps)

ALL_summary=[]
for SAMPLE in sample:
    for REP in reps:
        ALL_summary.append("6_summary/{sample}_{rep}/{sample}_{rep}_2k.matrix".format(sample=SAMPLE,rep=REP))
        ALL_summary.append("6_summary/{sample}_{rep}/{sample}_{rep}_2k_heatmap.pdf".format(sample=SAMPLE,rep=REP))
        ALL_summary.append("6_summary/{sample}_{rep}/{sample}_{rep}_2k_profile.pdf".format(sample=SAMPLE,rep=REP))

ALL_motif=expand("7_meme_analysis/{samples}_{rep}_{contrals}/",samples=samples,rep=reps,contrals=contrals)


TARGETS=[]
TARGETS.extend(ALL_filter_fastq)
TARGETS.extend(ALL_bam)
TARGETS.extend(ALL_bam_bai)
TARGETS.extend(ALL_peak)
TARGETS.extend(ALL_bigwig)
TARGETS.extend(ALL_motif)
TARGETS.extend(ALL_summary)


#all the rule begin
rule all:
    input:
        TARGETS

rule trimmomatic_filter:
    input:
        "raw_data/{sample}_{rep}.cleaned_1.fastq.gz",
        "raw_data/{sample}_{rep}.cleaned_2.fastq.gz"
    output:
        "1_clean_data/{sample}_{rep}_1.fastq.gz",
        "1_clean_data/{sample}_{rep}_1_unpaired.fastq.gz",
        "1_clean_data/{sample}_{rep}_2.fastq.gz",
        "1_clean_data/{sample}_{rep}_2_unpaired.fastq.gz"
    log:
        "1_clean_data/{sample}_{rep}_trimmomatic_filter.log"
    message: 
        "trimmomatic_filter"
    shell:
        "java -jar config[Trimmomatic] PE -phred33 \
        {input[0]} {input[1]} {output[0]} {output[1]} {output[2]} {output[3]}\
        ILLUMINACLIP:config[Truseq_fa]:2:30:10 \
        LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:50"


rule bowtie_alignment:
    input:
        "1_clean_data/{sample}_{rep}_1.fastq.gz",
        "1_clean_data/{sample}_{rep}_1_unpaired.fastq.gz",
        "1_clean_data/{sample}_{rep}_2.fastq.gz",
        "1_clean_data/{sample}_{rep}_2_unpaired.fastq.gz"
    output:
        "2_alignment/{sample}_{rep}.bam"
    log:
        "2_alignment/{sample}_{rep}_alignment.log"
    message:
        "bowtie2 alignment {input}"
    shell:
        "bowtie2 -p 4 config[Index] -1 {input[0]} -2 {input[2]} -U {input[1]} -U {input[3]} |samtools view -Sb -o {output} -"


rule samtools_sort:
    input:
        "2_alignment/{sample}_{rep}.bam"
    output:
        "3_filter/{sample}_{rep}_sorted.bam"
    log:
        "3_filter/{sample}_{rep}_sorted.log"
    message:
        "samtools sort {input} file"
    shell:
        "samtools view -@ 4 -Sb -q 30 {input} |samtools sort -o {output}"


rule samtools index:
    input:
        "3_filter/{sample}_{rep}_sorted.bam"
    output:
        "3_filter/{sample}_{rep}_sorted.bam.bai"
    message:
        "samtools build {input} index"
    log:
        "3_filter/{sample}_{rep}_sorted_bam_index.log"
    shell:
        "samtools index {input}"


rule GEM_call_peak:
    input:
        "3_filter/{samples}_{rep}_sorted.bam",
        "3_filter/{contrals}_{rep}_sorted.bam"
    output:
        directory("4_call_peak/{samples}_{rep}_{contrals}"),
        "4_call_peak/4_call_peak/{samples}_{rep}_{contrals}/{samples}_{rep}_{contrals}.GPS_events.narrowPeak"
    log:
        "4_call_peak/{samples}_{rep}_{contrals}_call_peak.log"
    shell:
        "java -Xmx6G -Xss4M -jar config[Gem]/gem.jar \
			--d config[Gem]/Read_Distribution_default.txt \
			--genome config[Reference_dir] \
			--g config[Chrom_size] \
			--expt {input[0]} \
			--ctrl {input[1]} \
			--q 5 --f SAM --k_min 6 --k_max 20 --outNP --sl --outHOMER \
			--out {output[0]}"


rule bamtools_bigwig:
    input:
        "3_filter/{sample}_{rep}_sorted.bam"
    output:
        "5_transcrate/{sample}_{rep}.bw"
    log:
        "5_transcrate/{sample}_{rep}_bigwig.log"
    shell:
        "bamCoverage -b {input} -o {output}"

rule get_matrix:
    input:
        "5_transcrate/{sample}_{rep}.bw"
    output:
        "6_summary/{sample}_{rep}/{sample}_{rep}_2k.matrix"
    log:
        "6_summary/{sample}_{rep}_get_matrix.log"
    shell:
        "computeMatrix scale-regions -S {input} -b 2000 -a 2000 -o {output} -R config[Gtf_dir]"


rule bamtools_summary:
    input:
        "6_summary/{sample}_{rep}/{sample}_{rep}_2k.matrix"
    output:
        "6_summary/{sample}_{rep}/{sample}_{rep}_2k_heatmap.pdf",
        "6_summary/{sample}_{rep}/{sample}_{rep}_2k_profile.pdf"
    log:
        "6_summary/{sample}_{rep}/{sample}_{rep}_bamtools_summary.log"
    shell:
        "plotHeatmap -m {input} --plotFileFormat pdf -o {output[0]} &&\
        plotProfile -m {input} --plotFileFormat pdf -o {output[1]}"


rule find_motif:
    input:
        "4_call_peak/4_call_peak/{samples}_{rep}_{contrals}/{samples}_{rep}_{contrals}.GPS_events.narrowPeak"
    output:
        directory("7_meme_analysis/{samples}_{rep}_{contrals}/")
    log:
        "7_meme_analysis/{samples}_{rep}_{contrals}/{samples}_{rep}_{contrals}_find_motif.log"
    shell:
        "cat {input} |sed 's/\s\+/\t/g' |sed 's/chr//g' | bedtools slop -i - -g config[Chrom_size] -l 50 -r 50 |findMotifsGenome.pl - config[reference] {output} -size 200"
