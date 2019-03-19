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
