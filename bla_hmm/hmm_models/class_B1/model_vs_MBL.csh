set bla=$1
set PFs=(class_B1 class_B2 class_B3 Group_00_MBL Group_01_MBL_UNIPROT_blasB Group_02_MBL_UNIPROT_GlyoxylaseII_FAMILY Group_03_MBL_UNIPROT_ZN_METALLOHYDROLASE_GROUP3 Group_04_MBL_UNIPROT_ATSA Group_05_MBL Group_06_MBL_UNIPROT_RNA_metabolizing_MBL Group_07_MBL_UNIPROT_DNA_REPAIR_MBL_FAMILY Group_08_MBL_IPR035681_ComA_like_MBL Group_09_MBL Group_10_MBL Group_10_MBL_Phnp_MBLs_IPR035682 Group_11_MBL_UNIPROT_CMP_NEU5_AC_HYDROXYLASE Group_12_MBL Group_13_MBL Group_14_MBL Group_15_MBL Group_16_MBL_UNIPROT_CAMP_PHOSPHODIESTERASE)
set cut_off=`sort -k2 -n scores.dat | head -1 | awk '{print int($2)}'`

foreach i ($PFs)
 hmmsearch class_"$bla"_model.hmm ../MBL_functionally/"$i".fasta > "$i"_vs_class_"$bla"_model.dat
 cat "$i"_vs_class_"$bla"_model.dat | sed -n '/Scores for complete sequences/,/Domain annotation for each/p' | sed '/------ inclusion threshold ------/d' | sed '1,4d' | head -n -3 > "$i"_vs_class_"$bla"_model.tab
 set seqs=`cat ../MBL_functionally/"$i".fasta | awk '$1 ~ ">" {print $1}' | wc -l`
 set hit_check=`cat "$i"_vs_class_"$bla"_model.tab | awk '$0 ~ "No hits detected that satisfy reporting thresholds" {print $0}' | wc -l`
 if ($hit_check == 0) then
 set hits=`cat "$i"_vs_class_"$bla"_model.tab | awk '{print $1}' | wc -l`
 set hits_co=`cat "$i"_vs_class_"$bla"_model.tab | awk -v c=$cut_off '$2 > c {print $2}' | wc -l`
 else
  set hits=0
  set hits_co=0
 endif
 echo "$i"" vs class_"$bla"_model cut-off=""$cut_off"" ""$seqs"" / ""$hits"" / ""$hits_co" >> class_"$bla"_model_vs_MBLs.dat
end

