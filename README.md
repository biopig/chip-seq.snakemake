# chip-seq.snakemake
that is the script of chip-seq analysis use the snakemake workflow

使用这个脚本需要安装的软件有：
- [python3]()
- [trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic)
- [bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml)
- [samtools](http://samtools.sourceforge.net/)
- [deeptools](https://deeptools.readthedocs.io/en/develop/)
- [gem](http://groups.csail.mit.edu/cgs/gem/)
- [homer](http://homer.ucsd.edu/homer/)

安装完软件以后需要下载基因组数据库，并使用bowtie2构建索引
> bowtie2-build reference.fa reference

reference目录中，chrom.sizes是染色体和长度的信息，使用tab分割；还需要各个染色体单独的序列，应该是这样的信息
```
├── chr10.fa
├── chr1.fa
├── chr2.fa
├── chr3.fa
├── chr4.fa
├── chr5.fa
├── chr6.fa
├── chr7.fa
├── chr8.fa
├── chr9.fa
├── chrom.sizes
├── Zea_mays.B73_RefGen_v4.dna.chromosome.fa
├── Zea_mays.B73_RefGen_v4.dna_index.1.bt2
├── Zea_mays.B73_RefGen_v4.dna_index.2.bt2
├── Zea_mays.B73_RefGen_v4.dna_index.3.bt2
├── Zea_mays.B73_RefGen_v4.dna_index.4.bt2
├── Zea_mays.B73_RefGen_v4.dna_index.rev.1.bt2
└── Zea_mays.B73_RefGen_v4.dna_index.rev.2.bt2
```
最后的结果应该是这样的一种情况
```
├── 1_clean_data
├── 2_alignment
├── 3_filter
├── 4_call_peak
│   ├── AG_1_HALO
│   │   └── AG_1_HALO_outputs
│   └── AG_2_HALO
│       └── AG_2_HALO_outputs
├── 5_transcrate
├── 6_summary
│   ├── AG_1
│   ├── AG_2
│   ├── HALO_1
│   └── HALO_2
├── 7_meme_analysis
│   ├── AG_1_HALO
│   │   ├── homerResults
│   │   └── knownResults
│   ├── AG_2_HALO
│   │   ├── homerResults
│   │   └── knownResults
│   └── get_fasta
├── 8_get_gene
│   ├── AG_1_HALO
│   └── AG_2_HALO
├── raw_data
├── reference
│   └── preparsed
└── script
```
