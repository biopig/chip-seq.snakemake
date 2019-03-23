这个目录存放的是**注释信息**
分别包括：
- **genome.fa**
- **genome.gtf**
- **indival.fa**
- **chrom.sizes**
- **genome.sqlit**
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
├── preparsed
│   ├── Zea_mays.B73_RefGen_v4.dna.chromosome.fa.200.cgbins
│   ├── Zea_mays.B73_RefGen_v4.dna.chromosome.fa.200.cgfreq
│   ├── Zea_mays.B73_RefGen_v4.dna.chromosome.fa.200.gcbins
│   ├── Zea_mays.B73_RefGen_v4.dna.chromosome.fa.200.pos
│   └── Zea_mays.B73_RefGen_v4.dna.chromosome.fa.200.seq
├── Zea_mays.B73_RefGen_v4.dna.chromosome.fa
├── Zea_mays.B73_RefGen_v4.dna_index.1.bt2
├── Zea_mays.B73_RefGen_v4.dna_index.2.bt2
├── Zea_mays.B73_RefGen_v4.dna_index.3.bt2
├── Zea_mays.B73_RefGen_v4.dna_index.4.bt2
├── Zea_mays.B73_RefGen_v4.dna_index.rev.1.bt2
└── Zea_mays.B73_RefGen_v4.dna_index.rev.2.bt2
```
注：
1. preparsed目录和里面的文件是中间分析时需要用到的文件，如果没有的话软件会根据genome.fa来生成的
2. \*.bt2是构建索引生成的文件
2. 因为之前我再分析数据的时候，gtf和sqlit不在这个文件夹中，所以这里没有显示，但是在使用时最好将注释文件都放在一个文件夹中，这样就比较清晰明了了
