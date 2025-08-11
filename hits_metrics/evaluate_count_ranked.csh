
  set input=$1
  set output=$2
  set class=$3
  set otu=$4
  set higher_otu=$5
  #count hits and frequencies
  set total_genomes=`cat $input | awk -v otu=$otu '$5 ~ otu"-" {print $1}' | sort -u | wc -l`
  if ($total_genomes == "0") then
   echo $otu no_genome_found >> $output/warnings.log
  endif

  set species=`cat $input | awk -v otu="$otu" '$5 ~ otu"-" {print $3}' | sort -u | wc -l`
  if ($species == "0") then
   echo $otu no_species_found >> $output/warnings.log
  endif

  set class_totals_hits=`cat $input | awk -v otu=$otu -v bla=$class '$5 ~ otu"-" && $8 == "class_"bla' | wc -l`

  set class_totals_genomes=`cat $input | awk -v otu=$otu -v bla=$class '$5 ~ otu"-" && $8 == "class_"bla {print $1}' | sort -u | wc -l`

  set class_totals_species=`cat $input | awk -v otu=$otu -v bla=$class '$5 ~ otu"-" && $8 == "class_"bla {print $3}' | sort -u | wc -l`

  set class_whole_hits=`cat $input | awk -v otu=$otu -v bla=$class '$5 ~ otu"-" && $8 == "class_"bla && $16 !~ "YYYYYY" && $16 !~ "PSEUDO" && $28 == 0 && $25 > 60 && $18 !~ "X" && $19 !~ "PSEUDO" ' | wc -l`

  set class_whole_genomes=`cat $input | awk -v otu=$otu -v bla=$class '$5 ~ otu"-" && $8 == "class_"bla && $16 !~ "YYYYYY" && $16 !~ "PSEUDO" && $28 == 0 && $25 > 60 && $18 !~ "X" && $19 !~ "PSEUDO" {print $1}' | sort -u | wc -l`

  set class_whole_species=`cat $input | awk -v otu=$otu -v bla=$class '$5 ~ otu"-" && $8 == "class_"bla && $16 !~ "YYYYYY" && $16 !~ "PSEUDO" && $28 == 0 && $25 > 60 && $18 !~ "X" && $19 !~ "PSEUDO" {print $3}' | sort -u | wc -l`
  
  set dif_clusters_cds=`cat $output/"$higher_otu""$otu"/"$otu"_"$class"_whole_proteins_positives_meta_ranked.dat | awk -v otu=$otu -v bla=$class '$6 ~ otu"-" && $9 == "class_"bla {print $31}' | sort -u | wc -l`

  set dif_clusters_proteins=`cat $output/"$higher_otu""$otu"/"$otu"_"$class"_whole_proteins_positives_meta_ranked.dat | awk -v otu=$otu -v bla=$class '$6 ~ otu"-" && $9 == "class_"bla {print $33}' | sort -u | wc -l`

  set class_totals_genomes_perc=`echo $class_totals_genomes $total_genomes | awk '{print $1/$2*100}'`

  set class_whole_genomes_perc=`echo $class_whole_genomes $total_genomes | awk '{print $1/$2*100}'`

  set class_totals_species_perc=`echo $class_totals_species $species | awk '{print $1/$2*100}'`

  set class_whole_species_perc=`echo $class_whole_species $species | awk '{print $1/$2*100}'`

  if ($class_whole_hits == "0") then
   set rel_richness_cds="NA"
   set rel_richness_proteins="NA"
  else
   set rel_richness_cds=`echo $dif_clusters_cds $class_whole_hits | awk '{print $1/$2}'`
   set rel_richness_proteins=`echo $dif_clusters_proteins $class_whole_hits | awk '{print $1/$2}'`
  endif
  
#Evaluation
  if ($class == "A") then
   set true_hits=`cat $output/"$higher_otu""$otu"/"$otu"_"$class"_whole_proteins_positives_meta_ranked.dat | awk -v otu=$otu -v bla=$class -v p1="._S_._._K" -v p2=E -v p3="[KR]_._._._." '$6 ~ otu"-" && $9 == "class_"bla && $38 ~ p1 && $40 == p2 && $41 ~ p3' | wc -l`
   set true_hits_genomes=`cat $output/"$higher_otu""$otu"/"$otu"_"$class"_whole_proteins_positives_meta_ranked.dat | awk -v otu=$otu -v bla=$class -v p1="._S_._._K" -v p2=E -v p3="[KR]_._._._." '$6 ~ otu"-" && $9 == "class_"bla && $38 ~ p1 && $40 == p2 && $41 ~ p3 {print $2}' | sort -u | wc -l`
   set true_hits_species=`cat $output/"$higher_otu""$otu"/"$otu"_"$class"_whole_proteins_positives_meta_ranked.dat | awk -v otu=$otu -v bla=$class -v p1="._S_._._K" -v p2=E -v p3="[KR]_._._._." '$6 ~ otu"-" && $9 == "class_"bla && $38 ~ p1 && $40 == p2 && $41 ~ p3 {print $4}' | sort -u | wc -l`
   set true_hits_cds_clu=`cat $output/"$higher_otu""$otu"/"$otu"_"$class"_whole_proteins_positives_meta_ranked.dat | awk -v otu=$otu -v bla=$class -v p1="._S_._._K" -v p2=E -v p3="[KR]_._._._." '$6 ~ otu"-" && $9 == "class_"bla && $38 ~ p1 && $40 == p2 && $41 ~ p3 {print $31}' | sort -u | wc -l`
   set true_hits_prot_clu=`cat $output/"$higher_otu""$otu"/"$otu"_"$class"_whole_proteins_positives_meta_ranked.dat | awk -v otu=$otu -v bla=$class -v p1="._S_._._K" -v p2=E -v p3="[KR]_._._._." '$6 ~ otu"-" && $9 == "class_"bla && $38 ~ p1 && $40 == p2 && $41 ~ p3 {print $33}' | sort -u | wc -l`
   set true_hits_secreted=`cat $output/"$higher_otu""$otu"/"$otu"_"$class"_whole_proteins_positives_meta_ranked.dat | awk -v otu=$otu -v bla=$class -v p1="._S_._._K" -v p2=E -v p3="[KR]_._._._." '$6 ~ otu"-" && $9 == "class_"bla && $38 ~ p1 && $40 == p2 && $41 ~ p3' | awk '$34 !~ "OTHER" || $36 !~ "OTHER"' | wc -l`
   set true_hits_secreted_genomes=`cat $output/"$higher_otu""$otu"/"$otu"_"$class"_whole_proteins_positives_meta_ranked.dat | awk -v otu=$otu -v bla=$class -v p1="._S_._._K" -v p2=E -v p3="[KR]_._._._." '$6 ~ otu"-" && $9 == "class_"bla && $38 ~ p1 && $40 == p2 && $41 ~ p3' | awk '$34 !~ "OTHER" || $36 !~ "OTHER" {print $2}' | sort -u | wc -l`
   set true_hits_secreted_species=`cat $output/"$higher_otu""$otu"/"$otu"_"$class"_whole_proteins_positives_meta_ranked.dat | awk -v otu=$otu -v bla=$class -v p1="._S_._._K" -v p2=E -v p3="[KR]_._._._." '$6 ~ otu"-" && $9 == "class_"bla && $38 ~ p1 && $40 == p2 && $41 ~ p3' | awk '$34 !~ "OTHER" || $36 !~ "OTHER" {print $4}' | sort -u | wc -l`
   set true_hits_secreted_cds_clu=`cat $output/"$higher_otu""$otu"/"$otu"_"$class"_whole_proteins_positives_meta_ranked.dat | awk -v otu=$otu -v bla=$class -v p1="._S_._._K" -v p2=E -v p3="[KR]_._._._." '$6 ~ otu"-" && $9 == "class_"bla && $38 ~ p1 && $40 == p2 && $41 ~ p3' | awk '$34 !~ "OTHER" || $36 !~ "OTHER" {print $31}' | sort -u | wc -l`
   set true_hits_secreted_prot_clu=`cat $output/"$higher_otu""$otu"/"$otu"_"$class"_whole_proteins_positives_meta_ranked.dat | awk -v otu=$otu -v bla=$class -v p1="._S_._._K" -v p2=E -v p3="[KR]_._._._." '$6 ~ otu"-" && $9 == "class_"bla && $38 ~ p1 && $40 == p2 && $41 ~ p3' | awk '$34 !~ "OTHER" || $36 !~ "OTHER" {print $33}' | sort -u | wc -l`
   set true_hits_genomes_perc=`echo $true_hits_genomes $total_genomes | awk '{print $1/$2*100}'`
   set true_hits_secreted_genomes_perc=`echo $true_hits_secreted_genomes $total_genomes | awk '{print $1/$2*100}'`
   set true_hits_species_perc=`echo $true_hits_species $species | awk '{print $1/$2*100}'`
   set true_hits_secreted_species_perc=`echo $true_hits_secreted_species $species | awk '{print $1/$2*100}'`
  else if ($class == "C") then
   set true_hits=`cat $output/"$higher_otu""$otu"/"$otu"_"$class"_whole_proteins_positives_meta_ranked.dat | awk -v otu=$otu -v bla=$class -v p1="S_._._K" -v p2="[KR]_._." '$6 ~ otu"-" && $9 == "class_"bla && $38 ~ p1 && $41 ~ p2' | wc -l`
   set true_hits_genomes=`cat $output/"$higher_otu""$otu"/"$otu"_"$class"_whole_proteins_positives_meta_ranked.dat | awk -v otu=$otu -v bla=$class -v p1="S_._._K" -v p2="[KR]_._." '$6 ~ otu"-" && $9 == "class_"bla && $38 ~ p1 && $41 ~ p2 {print $2}' | sort -u | wc -l`
   set true_hits_species=`cat $output/"$higher_otu""$otu"/"$otu"_"$class"_whole_proteins_positives_meta_ranked.dat | awk -v otu=$otu -v bla=$class -v p1="S_._._K" -v p2="[KR]_._." '$6 ~ otu"-" && $9 == "class_"bla && $38 ~ p1 && $41 ~ p2 {print $4}' | sort -u | wc -l`
   set true_hits_cds_clu=`cat $output/"$higher_otu""$otu"/"$otu"_"$class"_whole_proteins_positives_meta_ranked.dat | awk -v otu=$otu -v bla=$class -v p1="S_._._K" -v p2="[KR]_._." '$6 ~ otu"-" && $9 == "class_"bla && $38 ~ p1 && $41 ~ p2 {print $31}' | sort -u | wc -l`
   set true_hits_prot_clu=`cat $output/"$higher_otu""$otu"/"$otu"_"$class"_whole_proteins_positives_meta_ranked.dat | awk -v otu=$otu -v bla=$class -v p1="S_._._K" -v p2="[KR]_._." '$6 ~ otu"-" && $9 == "class_"bla && $38 ~ p1 && $41 ~ p2 {print $33}' | sort -u | wc -l`
   set true_hits_secreted=`cat $output/"$higher_otu""$otu"/"$otu"_"$class"_whole_proteins_positives_meta_ranked.dat | awk -v otu=$otu -v bla=$class -v p1="S_._._K" -v p2="[KR]_._." '$6 ~ otu"-" && $9 == "class_"bla && $38 ~ p1 && $41 ~ p2' | awk '$34 !~ "OTHER" || $36 !~ "OTHER"' | wc -l`
   set true_hits_secreted_genomes=`cat $output/"$higher_otu""$otu"/"$otu"_"$class"_whole_proteins_positives_meta_ranked.dat | awk -v otu=$otu -v bla=$class -v p1="S_._._K" -v p2="[KR]_._." '$6 ~ otu"-" && $9 == "class_"bla && $38 ~ p1 && $41 ~ p2' | awk '$34 !~ "OTHER" || $36 !~ "OTHER" {print $2}' | sort -u | wc -l`
   set true_hits_secreted_species=`cat $output/"$higher_otu""$otu"/"$otu"_"$class"_whole_proteins_positives_meta_ranked.dat | awk -v otu=$otu -v bla=$class -v p1="S_._._K" -v p2="[KR]_._." '$6 ~ otu"-" && $9 == "class_"bla && $38 ~ p1 && $41 ~ p2' | awk '$34 !~ "OTHER" || $36 !~ "OTHER" {print $4}' | sort -u | wc -l`
   set true_hits_secreted_cds_clu=`cat $output/"$higher_otu""$otu"/"$otu"_"$class"_whole_proteins_positives_meta_ranked.dat | awk -v otu=$otu -v bla=$class -v p1="S_._._K" -v p2="[KR]_._." '$6 ~ otu"-" && $9 == "class_"bla && $38 ~ p1 && $41 ~ p2' | awk '$34 !~ "OTHER" || $36 !~ "OTHER" {print $31}' | sort -u | wc -l`
   set true_hits_secreted_prot_clu=`cat $output/"$higher_otu""$otu"/"$otu"_"$class"_whole_proteins_positives_meta_ranked.dat | awk -v otu=$otu -v bla=$class -v p1="S_._._K" -v p2="[KR]_._." '$6 ~ otu"-" && $9 == "class_"bla && $38 ~ p1 && $41 ~ p2' | awk '$34 !~ "OTHER" || $36 !~ "OTHER" {print $33}' | sort -u | wc -l`
   set true_hits_genomes_perc=`echo $true_hits_genomes $total_genomes | awk '{print $1/$2*100}'`
   set true_hits_secreted_genomes_perc=`echo $true_hits_secreted_genomes $total_genomes | awk '{print $1/$2*100}'`
   set true_hits_species_perc=`echo $true_hits_species $species | awk '{print $1/$2*100}'`
   set true_hits_secreted_species_perc=`echo $true_hits_secreted_species $species | awk '{print $1/$2*100}'`
  else if ($class == "D") then
   set true_hits=`cat $output/"$higher_otu""$otu"/"$otu"_"$class"_whole_proteins_positives_meta_ranked.dat | awk -v otu=$otu -v bla=$class -v p1="S_._._K" -v p2="S_._[VILAGPFY]" -v p3="[KR]_._." '$6 ~ otu"-" && $9 == "class_"bla && $21 !~ "BlaR1" && $38 ~ p1 && $39 ~ p2 && $40 ~ p3' | wc -l`
   set true_hits_genomes=`cat $output/"$higher_otu""$otu"/"$otu"_"$class"_whole_proteins_positives_meta_ranked.dat | awk -v otu=$otu -v bla=$class -v p1="S_._._K" -v p2="S_._[VILAGPFY]" -v p3="[KR]_._." '$6 ~ otu"-" && $9 == "class_"bla && $21 !~ "BlaR1" && $38 ~ p1 && $39 ~ p2 && $40 ~ p3 {print $2}' | sort -u | wc -l`
   set true_hits_species=`cat $output/"$higher_otu""$otu"/"$otu"_"$class"_whole_proteins_positives_meta_ranked.dat | awk -v otu=$otu -v bla=$class -v p1="S_._._K" -v p2="S_._[VILAGPFY]" -v p3="[KR]_._." '$6 ~ otu"-" && $9 == "class_"bla && $21 !~ "BlaR1" && $38 ~ p1 && $39 ~ p2 && $40 ~ p3 {print $4}' | sort -u | wc -l`
   set true_hits_cds_clu=`cat $output/"$higher_otu""$otu"/"$otu"_"$class"_whole_proteins_positives_meta_ranked.dat | awk -v otu=$otu -v bla=$class -v p1="S_._._K" -v p2="S_._[VILAGPFY]" -v p3="[KR]_._." '$6 ~ otu"-" && $9 == "class_"bla && $21 !~ "BlaR1" && $38 ~ p1 && $39 ~ p2 && $40 ~ p3 {print $31}' | sort -u | wc -l`
   set true_hits_prot_clu=`cat $output/"$higher_otu""$otu"/"$otu"_"$class"_whole_proteins_positives_meta_ranked.dat | awk -v otu=$otu -v bla=$class -v p1="S_._._K" -v p2="S_._[VILAGPFY]" -v p3="[KR]_._." '$6 ~ otu"-" && $9 == "class_"bla && $21 !~ "BlaR1" && $38 ~ p1 && $39 ~ p2 && $40 ~ p3 {print $33}' | sort -u | wc -l`
   set true_hits_secreted=`cat $output/"$higher_otu""$otu"/"$otu"_"$class"_whole_proteins_positives_meta_ranked.dat | awk -v otu=$otu -v bla=$class -v p1="S_._._K" -v p2="S_._[VILAGPFY]" -v p3="[KR]_._." '$6 ~ otu"-" && $9 == "class_"bla && $21 !~ "BlaR1" && $38 ~ p1 && $39 ~ p2 && $40 ~ p3' | awk '$34 !~ "OTHER" || $36 !~ "OTHER"' | wc -l`
   set true_hits_secreted_genomes=`cat $output/"$higher_otu""$otu"/"$otu"_"$class"_whole_proteins_positives_meta_ranked.dat | awk -v otu=$otu -v bla=$class -v p1="S_._._K" -v p2="S_._[VILAGPFY]" -v p3="[KR]_._." '$6 ~ otu"-" && $9 == "class_"bla && $21 !~ "BlaR1" && $38 ~ p1 && $39 ~ p2 && $40 ~ p3' | awk '$34 !~ "OTHER" || $36 !~ "OTHER" {print $2}' | sort -u | wc -l`
   set true_hits_secreted_species=`cat $output/"$higher_otu""$otu"/"$otu"_"$class"_whole_proteins_positives_meta_ranked.dat | awk -v otu=$otu -v bla=$class -v p1="S_._._K" -v p2="S_._[VILAGPFY]" -v p3="[KR]_._." '$6 ~ otu"-" && $9 == "class_"bla && $21 !~ "BlaR1" && $38 ~ p1 && $39 ~ p2 && $40 ~ p3' | awk '$34 !~ "OTHER" || $36 !~ "OTHER" {print $4}' | sort -u | wc -l`
   set true_hits_secreted_cds_clu=`cat $output/"$higher_otu""$otu"/"$otu"_"$class"_whole_proteins_positives_meta_ranked.dat | awk -v otu=$otu -v bla=$class -v p1="S_._._K" -v p2="S_._[VILAGPFY]" -v p3="[KR]_._." '$6 ~ otu"-" && $9 == "class_"bla && $21 !~ "BlaR1" && $38 ~ p1 && $39 ~ p2 && $40 ~ p3' | awk '$34 !~ "OTHER" || $36 !~ "OTHER" {print $31}' | sort -u | wc -l`
   set true_hits_secreted_prot_clu=`cat $output/"$higher_otu""$otu"/"$otu"_"$class"_whole_proteins_positives_meta_ranked.dat | awk -v otu=$otu -v bla=$class -v p1="S_._._K" -v p2="S_._[VILAGPFY]" -v p3="[KR]_._." '$6 ~ otu"-" && $9 == "class_"bla && $21 !~ "BlaR1" && $38 ~ p1 && $39 ~ p2 && $40 ~ p3' | awk '$34 !~ "OTHER" || $36 !~ "OTHER" {print $33}' | sort -u | wc -l`
   set true_hits_genomes_perc=`echo $true_hits_genomes $total_genomes | awk '{print $1/$2*100}'`
   set true_hits_secreted_genomes_perc=`echo $true_hits_secreted_genomes $total_genomes | awk '{print $1/$2*100}'`
   set true_hits_species_perc=`echo $true_hits_species $species | awk '{print $1/$2*100}'`
   set true_hits_secreted_species_perc=`echo $true_hits_secreted_species $species | awk '{print $1/$2*100}'`
  else if ($class == "B1") then
   set true_hits=`cat $output/"$higher_otu""$otu"/"$otu"_"$class"_whole_proteins_positives_meta_ranked.dat | awk -v otu=$otu -v bla=$class -v p1="H_._H_._D" -v p2="H_._._._." -v p3=C -v p4=H '$6 ~ otu"-" && $9 == "class_"bla && $21 !~ "B2_POSITIVE" && $38 ~ p1 && $39 ~ p2 && $40 == p3 && $41 == p4' | wc -l`
   set true_hits_genomes=`cat $output/"$higher_otu""$otu"/"$otu"_"$class"_whole_proteins_positives_meta_ranked.dat | awk -v otu=$otu -v bla=$class -v p1="H_._H_._D" -v p2="H_._._._." -v p3=C -v p4=H '$6 ~ otu"-" && $9 == "class_"bla && $21 !~ "B2_POSITIVE" && $38 ~ p1 && $39 ~ p2 && $40 == p3 && $41 == p4 {print $2}'| sort -u | wc -l`
   set true_hits_species=`cat $output/"$higher_otu""$otu"/"$otu"_"$class"_whole_proteins_positives_meta_ranked.dat | awk -v otu=$otu -v bla=$class -v p1="H_._H_._D" -v p2="H_._._._." -v p3=C -v p4=H '$6 ~ otu"-" && $9 == "class_"bla && $21 !~ "B2_POSITIVE" && $38 ~ p1 && $39 ~ p2 && $40 == p3 && $41 == p4 {print $4}'| sort -u | wc -l`
   set true_hits_cds_clu=`cat $output/"$higher_otu""$otu"/"$otu"_"$class"_whole_proteins_positives_meta_ranked.dat | awk -v otu=$otu -v bla=$class -v p1="H_._H_._D" -v p2="H_._._._." -v p3=C -v p4=H '$6 ~ otu"-" && $9 == "class_"bla && $21 !~ "B2_POSITIVE" && $38 ~ p1 && $39 ~ p2 && $40 == p3 && $41 == p4 {print $31}'| sort -u | wc -l`
   set true_hits_prot_clu=`cat $output/"$higher_otu""$otu"/"$otu"_"$class"_whole_proteins_positives_meta_ranked.dat | awk -v otu=$otu -v bla=$class -v p1="H_._H_._D" -v p2="H_._._._." -v p3=C -v p4=H '$6 ~ otu"-" && $9 == "class_"bla && $21 !~ "B2_POSITIVE" && $38 ~ p1 && $39 ~ p2 && $40 == p3 && $41 == p4 {print $33}'| sort -u | wc -l`
   set true_hits_secreted=`cat $output/"$higher_otu""$otu"/"$otu"_"$class"_whole_proteins_positives_meta_ranked.dat | awk -v otu=$otu -v bla=$class -v p1="H_._H_._D" -v p2="H_._._._." -v p3=C -v p4=H '$6 ~ otu"-" && $9 == "class_"bla && $21 !~ "B2_POSITIVE" && $38 ~ p1 && $39 ~ p2 && $40 == p3 && $41 == p4' | awk '$34 !~ "OTHER" || $36 !~ "OTHER"' | wc -l`
   set true_hits_secreted_genomes=`cat $output/"$higher_otu""$otu"/"$otu"_"$class"_whole_proteins_positives_meta_ranked.dat | awk -v otu=$otu -v bla=$class -v p1="H_._H_._D" -v p2="H_._._._." -v p3=C -v p4=H '$6 ~ otu"-" && $9 == "class_"bla && $21 !~ "B2_POSITIVE" && $38 ~ p1 && $39 ~ p2 && $40 == p3 && $41 == p4' | awk '$34 !~ "OTHER" || $36 !~ "OTHER" {print $2}'| sort -u | wc -l`
   set true_hits_secreted_species=`cat $output/"$higher_otu""$otu"/"$otu"_"$class"_whole_proteins_positives_meta_ranked.dat | awk -v otu=$otu -v bla=$class -v p1="H_._H_._D" -v p2="H_._._._." -v p3=C -v p4=H '$6 ~ otu"-" && $9 == "class_"bla && $21 !~ "B2_POSITIVE" && $38 ~ p1 && $39 ~ p2 && $40 == p3 && $41 == p4' | awk '$34 !~ "OTHER" || $36 !~ "OTHER" {print $4}'| sort -u | wc -l`
   set true_hits_secreted_cds_clu=`cat $output/"$higher_otu""$otu"/"$otu"_"$class"_whole_proteins_positives_meta_ranked.dat | awk -v otu=$otu -v bla=$class -v p1="H_._H_._D" -v p2="H_._._._." -v p3=C -v p4=H '$6 ~ otu"-" && $9 == "class_"bla && $21 !~ "B2_POSITIVE" && $38 ~ p1 && $39 ~ p2 && $40 == p3 && $41 == p4' | awk '$34 !~ "OTHER" || $36 !~ "OTHER" {print $31}'| sort -u | wc -l`
   set true_hits_secreted_prot_clu=`cat $output/"$higher_otu""$otu"/"$otu"_"$class"_whole_proteins_positives_meta_ranked.dat | awk -v otu=$otu -v bla=$class -v p1="H_._H_._D" -v p2="H_._._._." -v p3=C -v p4=H '$6 ~ otu"-" && $9 == "class_"bla && $21 !~ "B2_POSITIVE" && $38 ~ p1 && $39 ~ p2 && $40 == p3 && $41 == p4' | awk '$34 !~ "OTHER" || $36 !~ "OTHER" {print $33}'| sort -u | wc -l`
   set true_hits_genomes_perc=`echo $true_hits_genomes $total_genomes | awk '{print $1/$2*100}'`
   set true_hits_secreted_genomes_perc=`echo $true_hits_secreted_genomes $total_genomes | awk '{print $1/$2*100}'`
   set true_hits_species_perc=`echo $true_hits_species $species | awk '{print $1/$2*100}'`
   set true_hits_secreted_species_perc=`echo $true_hits_secreted_species $species | awk '{print $1/$2*100}'`
  else if ($class == "B2") then
   set true_hits=`cat $output/"$higher_otu""$otu"/"$otu"_"$class"_whole_proteins_positives_meta_ranked.dat | awk -v otu=$otu -v bla=$class -v p1="._._._._D" -v p2=C -v p3=H '$6 ~ otu"-" && $9 == "class_"bla && $38 ~ p1 && $39 == p2 && $40 == p3' | wc -l`
   set true_hits_genomes=`cat $output/"$higher_otu""$otu"/"$otu"_"$class"_whole_proteins_positives_meta_ranked.dat | awk -v otu=$otu -v bla=$class -v p1="._._._._D" -v p2=C -v p3=H '$6 ~ otu"-" && $9 == "class_"bla && $38 ~ p1 && $39 == p2 && $40 == p3 {print $2}' | sort -u | wc -l`
   set true_hits_species=`cat $output/"$higher_otu""$otu"/"$otu"_"$class"_whole_proteins_positives_meta_ranked.dat | awk -v otu=$otu -v bla=$class -v p1="._._._._D" -v p2=C -v p3=H '$6 ~ otu"-" && $9 == "class_"bla && $38 ~ p1 && $39 == p2 && $40 == p3 {print $4}' | sort -u | wc -l`
   set true_hits_cds_clu=`cat $output/"$higher_otu""$otu"/"$otu"_"$class"_whole_proteins_positives_meta_ranked.dat | awk -v otu=$otu -v bla=$class -v p1="._._._._D" -v p2=C -v p3=H '$6 ~ otu"-" && $9 == "class_"bla && $38 ~ p1 && $39 == p2 && $40 == p3 {print $31}' | sort -u | wc -l`
   set true_hits_prot_clu=`cat $output/"$higher_otu""$otu"/"$otu"_"$class"_whole_proteins_positives_meta_ranked.dat | awk -v otu=$otu -v bla=$class -v p1="._._._._D" -v p2=C -v p3=H '$6 ~ otu"-" && $9 == "class_"bla && $38 ~ p1 && $39 == p2 && $40 == p3 {print $33}' | sort -u | wc -l`
   set true_hits_secreted=`cat $output/"$higher_otu""$otu"/"$otu"_"$class"_whole_proteins_positives_meta_ranked.dat | awk -v otu=$otu -v bla=$class -v p1="._._._._D" -v p2=C -v p3=H '$6 ~ otu"-" && $9 == "class_"bla && $38 ~ p1 && $39 == p2 && $40 == p3' | awk '$34 !~ "OTHER" || $36 !~ "OTHER"' | wc -l`
   set true_hits_secreted_genomes=`cat $output/"$higher_otu""$otu"/"$otu"_"$class"_whole_proteins_positives_meta_ranked.dat | awk -v otu=$otu -v bla=$class -v p1="._._._._D" -v p2=C -v p3=H '$6 ~ otu"-" && $9 == "class_"bla && $38 ~ p1 && $39 == p2 && $40 == p3' | awk '$34 !~ "OTHER" || $36 !~ "OTHER" {print $2}' | sort -u | wc -l`
   set true_hits_secreted_species=`cat $output/"$higher_otu""$otu"/"$otu"_"$class"_whole_proteins_positives_meta_ranked.dat | awk -v otu=$otu -v bla=$class -v p1="._._._._D" -v p2=C -v p3=H '$6 ~ otu"-" && $9 == "class_"bla && $38 ~ p1 && $39 == p2 && $40 == p3' | awk '$34 !~ "OTHER" || $36 !~ "OTHER" {print $4}' | sort -u | wc -l`
   set true_hits_secreted_cds_clu=`cat $output/"$higher_otu""$otu"/"$otu"_"$class"_whole_proteins_positives_meta_ranked.dat | awk -v otu=$otu -v bla=$class -v p1="._._._._D" -v p2=C -v p3=H '$6 ~ otu"-" && $9 == "class_"bla && $38 ~ p1 && $39 == p2 && $40 == p3' | awk '$34 !~ "OTHER" || $36 !~ "OTHER" {print $31}' | sort -u | wc -l`
   set true_hits_secreted_prot_clu=`cat $output/"$higher_otu""$otu"/"$otu"_"$class"_whole_proteins_positives_meta_ranked.dat | awk -v otu=$otu -v bla=$class -v p1="._._._._D" -v p2=C -v p3=H '$6 ~ otu"-" && $9 == "class_"bla && $38 ~ p1 && $39 == p2 && $40 == p3' | awk '$34 !~ "OTHER" || $36 !~ "OTHER" {print $33}' | sort -u | wc -l`
   set true_hits_genomes_perc=`echo $true_hits_genomes $total_genomes | awk '{print $1/$2*100}'`
   set true_hits_secreted_genomes_perc=`echo $true_hits_secreted_genomes $total_genomes | awk '{print $1/$2*100}'`
   set true_hits_species_perc=`echo $true_hits_species $species | awk '{print $1/$2*100}'`
   set true_hits_secreted_species_perc=`echo $true_hits_secreted_species $species | awk '{print $1/$2*100}'`
  else if ($class == "B3") then
   set true_hits=`cat $output/"$higher_otu""$otu"/"$otu"_"$class"_whole_proteins_positives_meta_ranked.dat | awk -v otu=$otu -v bla=$class -v p1="[HQ]_._H_._D_H" -v p2=H -v p3=H '$6 ~ otu"-" && $9 == "class_"bla && $21 !~ "HARLDQ_MOTIF" && $38 ~ p1 && $39 == p2 && $40 == p3' | wc -l`
   set true_hits_genomes=`cat $output/"$higher_otu""$otu"/"$otu"_"$class"_whole_proteins_positives_meta_ranked.dat | awk -v otu=$otu -v bla=$class -v p1="[HQ]_._H_._D_H" -v p2=H -v p3=H '$6 ~ otu"-" && $9 == "class_"bla && $21 !~ "HARLDQ_MOTIF" && $38 ~ p1 && $39 == p2 && $40 == p3 {print $2}' | sort -u | wc -l`
   set true_hits_species=`cat $output/"$higher_otu""$otu"/"$otu"_"$class"_whole_proteins_positives_meta_ranked.dat | awk -v otu=$otu -v bla=$class -v p1="[HQ]_._H_._D_H" -v p2=H -v p3=H '$6 ~ otu"-" && $9 == "class_"bla && $21 !~ "HARLDQ_MOTIF" && $38 ~ p1 && $39 == p2 && $40 == p3 {print $4}' | sort -u | wc -l`
   set true_hits_cds_clu=`cat $output/"$higher_otu""$otu"/"$otu"_"$class"_whole_proteins_positives_meta_ranked.dat | awk -v otu=$otu -v bla=$class -v p1="[HQ]_._H_._D_H" -v p2=H -v p3=H '$6 ~ otu"-" && $9 == "class_"bla && $21 !~ "HARLDQ_MOTIF" && $38 ~ p1 && $39 == p2 && $40 == p3 {print $31}' | sort -u | wc -l`
   set true_hits_prot_clu=`cat $output/"$higher_otu""$otu"/"$otu"_"$class"_whole_proteins_positives_meta_ranked.dat | awk -v otu=$otu -v bla=$class -v p1="[HQ]_._H_._D_H" -v p2=H -v p3=H '$6 ~ otu"-" && $9 == "class_"bla && $21 !~ "HARLDQ_MOTIF" && $38 ~ p1 && $39 == p2 && $40 == p3 {print $33}' | sort -u | wc -l`
   set true_hits_secreted=`cat $output/"$higher_otu""$otu"/"$otu"_"$class"_whole_proteins_positives_meta_ranked.dat | awk -v otu=$otu -v bla=$class -v p1="[HQ]_._H_._D_H" -v p2=H -v p3=H '$6 ~ otu"-" && $9 == "class_"bla && $21 !~ "HARLDQ_MOTIF" && $38 ~ p1 && $39 == p2 && $40 == p3' | awk '$34 !~ "OTHER" || $36 !~ "OTHER"' | wc -l`
   set true_hits_secreted_genomes=`cat $output/"$higher_otu""$otu"/"$otu"_"$class"_whole_proteins_positives_meta_ranked.dat | awk -v otu=$otu -v bla=$class -v p1="[HQ]_._H_._D_H" -v p2=H -v p3=H '$6 ~ otu"-" && $9 == "class_"bla && $21 !~ "HARLDQ_MOTIF" && $38 ~ p1 && $39 == p2 && $40 == p3' | awk '$34 !~ "OTHER" || $36 !~ "OTHER" {print $2}' | sort -u | wc -l`
   set true_hits_secreted_species=`cat $output/"$higher_otu""$otu"/"$otu"_"$class"_whole_proteins_positives_meta_ranked.dat | awk -v otu=$otu -v bla=$class -v p1="[HQ]_._H_._D_H" -v p2=H -v p3=H '$6 ~ otu"-" && $9 == "class_"bla && $21 !~ "HARLDQ_MOTIF" && $38 ~ p1 && $39 == p2 && $40 == p3' | awk '$34 !~ "OTHER" || $36 !~ "OTHER" {print $4}' | sort -u | wc -l`
   set true_hits_secreted_cds_clu=`cat $output/"$higher_otu""$otu"/"$otu"_"$class"_whole_proteins_positives_meta_ranked.dat | awk -v otu=$otu -v bla=$class -v p1="[HQ]_._H_._D_H" -v p2=H -v p3=H '$6 ~ otu"-" && $9 == "class_"bla && $21 !~ "HARLDQ_MOTIF" && $38 ~ p1 && $39 == p2 && $40 == p3' | awk '$34 !~ "OTHER" || $36 !~ "OTHER" {print $31}' | sort -u | wc -l`
   set true_hits_secreted_prot_clu=`cat $output/"$higher_otu""$otu"/"$otu"_"$class"_whole_proteins_positives_meta_ranked.dat | awk -v otu=$otu -v bla=$class -v p1="[HQ]_._H_._D_H" -v p2=H -v p3=H '$6 ~ otu"-" && $9 == "class_"bla && $21 !~ "HARLDQ_MOTIF" && $38 ~ p1 && $39 == p2 && $40 == p3' | awk '$34 !~ "OTHER" || $36 !~ "OTHER" {print $33}' | sort -u | wc -l`
   set true_hits_genomes_perc=`echo $true_hits_genomes $total_genomes | awk '{print $1/$2*100}'`
   set true_hits_secreted_genomes_perc=`echo $true_hits_secreted_genomes $total_genomes | awk '{print $1/$2*100}'`
   set true_hits_species_perc=`echo $true_hits_species $species | awk '{print $1/$2*100}'`
   set true_hits_secreted_species_perc=`echo $true_hits_secreted_species $species | awk '{print $1/$2*100}'`
  endif
  # evaluation end

   echo $total_genomes $species $class_totals_hits $class_totals_genomes $class_totals_genomes_perc $class_totals_species $class_totals_species_perc $class_whole_hits $class_whole_genomes $class_whole_genomes_perc $class_whole_species $class_whole_species_perc $dif_clusters_cds $rel_richness_cds $dif_clusters_proteins $rel_richness_proteins $true_hits $true_hits_genomes $true_hits_genomes_perc $true_hits_species $true_hits_species_perc $true_hits_cds_clu $true_hits_prot_clu $true_hits_secreted $true_hits_secreted_genomes $true_hits_secreted_genomes_perc $true_hits_secreted_species $true_hits_secreted_species_perc $true_hits_secreted_cds_clu $true_hits_secreted_prot_clu >> $output/"$higher_otu""$otu"/"$otu"_class_"$class"_statistics.out

