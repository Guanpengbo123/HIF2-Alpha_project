source /etc/bashrc
#!/bin/bash
gmtlist='c2.cp.reactome.v7.5.1.symbols c5.go.v7.5.1.symbols Total_kegg'
for gmt in ${gmtlist}
do	
	gsea-cli.sh GSEA -res ./matrix.txt \
 -cls ./matrix.cls#type1_versus_type2 \
 -gmx ./${gmt}.gmt \
 -collapse Collapse -mode Max_probe -norm meandiv -nperm 1000 -permute gene_set -rnd_type no_balance -scoring_scheme weighted  -metric Diff_of_Classes -sort real -order descending -chip ./Mouse_Gene_Symbol_Remapping_Human_Orthologs_MSigDB.v7.4.chip \
 -create_gcts false -create_svgs false -include_only_symbols true -make_sets true -median false -num 100 -plot_top_x 20 -rnd_seed timestamp -save_rnd_lists false -set_max 500 -set_min 5 \
-zip_report false  -rpt_label ${gmt}
done
