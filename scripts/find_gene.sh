dir=/disks/wanjiong/dap


for i in AG_1 AG_2 GG_1 GG_2 AY_1 AY_2 GY_1 GY_2
do
        if [ ! -d "$dir"/8_get_gene/"$i" ]
        then mkdir "$dir"/8_get_gene/"$i"
        fi
        Rscript find_gene.R "$dir"/8_get_gene/maize_v4_2018_11.sqlite "$dir"/7_meme_analysis/get_fasta/"$i".GPS_peaks.bed "$i"
done
