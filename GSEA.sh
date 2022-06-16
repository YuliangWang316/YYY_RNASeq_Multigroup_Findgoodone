data_dir=/mnt/d/yyy/Merge_HUADA_YK_HUADA2/3vs3/gct1
for file_a in ${data_dir}/*;do
bash /mnt/d/GSEA_4.1.0/gsea-cli.sh GSEA -res $file_a -cls /mnt/d/yyy/sample3_vs_3.cls#KO_versus_WT -gmx  /mnt/d/h.all.v7.1.symbols.gmt  -collapse false -nperm 1000 -permute gene_set -rnd_type no_balance -scoring_scheme weighted -set_max 5000 -set_min 0 -rpt_label file -metric Signal2Noise -sort real -order descending -out /mnt/d/yyy/Merge_HUADA_YK_HUADA2/3vs3/result1
done

