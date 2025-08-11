#script for generating cluster annotation on EXABAYES tree

set cluster_data_file=$1
set phylip_name=$2
set clusters_exa=`cat $cluster_data_file | awk '{print $1}'`

foreach i ($clusters_exa)
 set clu_exabayes=`echo $i | sed 's/Cluster_/'$phylip_name'_/g'`
 set members_transf=`cat $cluster_data_file | awk -v cl=$i '$1 == cl {print log($2)/log(10)+1}'`
 set num_genera=`cat $cluster_data_file | awk -v cl=$i '$1 == cl {print log($6)/log(10)+1}'`
 set phyla=`cat $cluster_data_file | awk -v cl=$i '$1 == cl {print $15}'`
 set phyla_num=`cat $cluster_data_file | awk -v cl=$i '$1 == cl {print $14}'`
 set classes=`cat $cluster_data_file | awk -v cl=$i '$1 == cl {print $13}'`
 set classes_num=`cat $cluster_data_file | awk -v cl=$i '$1 == cl {print $12}'`
 
 set motifs_check=`cat $cluster_data_file | sed 's/gnc_region_disrupted_incoclusive/-/g' | sed 's/dhy_region_disrupted_incoclusive/-/g' | awk -v cl=$i -v p1="._._._._D" -v p2=C -v p3=H '$1 == cl && $19 ~ p1 && $20 ~ p2 && $21 ~ p3 {print $19,$20,$21}'`
 set motifs_check_num=`cat $cluster_data_file | sed 's/gnc_region_disrupted_incoclusive/-/g' | sed 's/dhy_region_disrupted_incoclusive/-/g' | awk -v cl=$i -v p1="._._._._D" -v p2=C -v p3=H '$1 == cl && $19 ~ p1 && $20 ~ p2 && $21 ~ p3 {print $19,$20,$21}' | wc -l`

 if ($motifs_check_num != "0") then
  # co-existence of atypical motifs 
  set motif_1=`echo $motifs_check | awk '{print $1}' | sed 's/,/ /g' |  awk -v p1="._._._._D" '{for(i=1;i<=NF;i++){if($i !~ p1) print $i;}}' | wc -l`
  set motif_2=`echo $motifs_check | awk '{print $2}' | sed 's/,/ /g' |  awk -v p2=C '{for(i=1;i<=NF;i++){if($i != p2) print $i;}}' | wc -l`
  set motif_3=`echo $motifs_check | awk '{print $3}' | sed 's/,/ /g' |  awk -v p3=H '{for(i=1;i<=NF;i++){if($i != p3) print $i;}}' | wc -l`
  if ($motif_1 != "0" || $motif_2 != "0" || $motif_3 != "0") then
   set triangle="0"
  else
   set triangle="1"
  endif
 else
  set triangle="-1"
 endif
 set secretion_check=`cat $cluster_data_file | awk -v cl=$i '$1 == cl && $23 == "OTHER" && $24 == "OTHER"' | wc -l`
 if ($secretion_check != "0") then
  set circle="-1"
 else
 set circle="1"
 endif
 set characterization=`cat $cluster_data_file | awk -v cl=$i '$1 == cl && $25 == "YES"' | wc -l`
 if ($characterization != "0") then
  set star="1"
 else
  set star="-1"
 endif

 # remember to put space as SEPERATOR
 #1 for multibar chart
 echo $clu_exabayes $members_transf $num_genera >> for_multibar.dat
 #2 for phyla and classes colored strips
 if ($phyla_num > 1) then
  echo $clu_exabayes "#000000" Multiple_Phyla >> for_phyla_color_strip.dat
 else
  set color=`cat color_codes_phyla.dat | awk -v cla=$phyla '$1 == cla {print $2}'`
  echo $clu_exabayes "$color" $phyla >> for_phyla_color_strip.dat # auto 8a alla3ei. giati?
 endif
 if ($classes_num > 1) then
  echo $clu_exabayes "#000000" Multiple_Classes >> for_classes_color_strip.dat
 else
  set color=`cat B2color_codes_classes.dat | awk -v cla=$classes '$1 == cla {print $2}'`
  echo $clu_exabayes "$color" $classes >> for_classes_color_strip.dat # auto 8a alla3ei
 endif
 #3 for binary
 echo $clu_exabayes $triangle $circle $star >> for_binary.dat
end

