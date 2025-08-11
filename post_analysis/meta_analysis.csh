#meta analysis scirpt of whole genes/proteins. Catalytic motifs and signal peptide prediction
set class=$1

set results=$2 #file used as input for cluster analysis

cat "$results"_class_"$class"/"$results"_class_"$class"_wcodes.dat | awk '{print $0,length($19)}' > "$results"_class_"$class"/"$results"_class_"$class"_wcodes_plen.dat

set entries_whole=`cat "$results"_class_"$class"/"$class"_whole_proteins_code.dat | awk '{print $1}' | sort -u`

# add clustering, signalP, motifs information
foreach j ($entries_whole)
 echo $j
 set cluster_cds=`cat "$results"_class_"$class"/class_"$class"_clusters_cds.dat | awk -v ent=$j '$3 == ent {print $1}'`
 set ref_db_compar=`cat "$results"_class_"$class"/class_"$class"_clusters_cds_wrefs.dat | awk -v ent=$j '$3 == ent {print $4}'`
 set cluster_proteins=`cat "$results"_class_"$class"/class_"$class"_clusters_proteins.dat | awk -v ent=$j '$3 == ent {print $1}'`
 set gr_neg_spt=`cat "$results"_class_"$class"/gramneg_summary.signalp5 | awk -v ent=$j '$1 == ent {print $2}' | sed 's/(/_/g' | sed 's/)/_/g'`
 if ($gr_neg_spt == "OTHER") then
  set gr_neg_site=0
 else
  set gr_neg_site=`cat "$results"_class_"$class"/gramneg_summary.signalp5 | awk -v ent=$j '$1 == ent {print $9}' | sed 's/\.//g'`
 endif
 set gr_pos_spt=`cat "$results"_class_"$class"/grampos_summary.signalp5 | awk -v ent=$j '$1 == ent {print $2}' | sed 's/(/_/g' | sed 's/)/_/g'`
 if ($gr_pos_spt == "OTHER") then
  set gr_pos_site=0
 else
 set gr_pos_site=`cat "$results"_class_"$class"/grampos_summary.signalp5 | awk -v ent=$j '$1 == ent {print $9}' | sed 's/\.//g'`
 endif
 set aligned_profile=`cat "$results"_class_"$class"/"$results"_class_"$class"_wcodes_plen.dat | awk -v ent=$j '$1 == ent {print $22}'`
 set aligned_query=`cat "$results"_class_"$class"/"$results"_class_"$class"_wcodes_plen.dat | awk -v ent=$j '$1 == ent {print $23}'`
# CLASS B1 MOTIFS CHECK
 if ($class == "B1") then
  #motif 1
  set m1="HfHeD"
  #motif 2
  set m2="HtkDN"
  #motif 3
  set m3="GGC"
  #motif 4
  set m4="GH"
  set m1_check=`echo $aligned_profile | awk -v m1=$m1 '$1 ~ m1 {print $1}' | wc -l`
   if ($m1_check == "0") then
    set comment1=$m1"_region_disrupted_incoclusive"
   else
    set H1=`echo $aligned_profile | awk -v m1=$m1 '{print match($1,m1)}'`
    set H1hit=`echo $aligned_query | cut -c$H1-$H1`
    set H2=`echo $aligned_profile | awk -v m1=$m1 '{print match($1,m1)+2}'`
    set H2hit=`echo $aligned_query | cut -c$H2-$H2`
    set D3=`echo $aligned_profile | awk -v m1=$m1 '{print match($1,m1)+4}'`
    set D3hit=`echo $aligned_query | cut -c$D3-$D3`
    set comment1=$H1hit"_x_"$H2hit"_x_"$D3hit
   endif
  set m2_check=`echo $aligned_profile | awk -v m2=$m2 '$1 ~ m2 {print $1}' | wc -l`
   if ($m2_check == "0") then
    set comment2=$m2"_region_disrupted_incoclusive" 
   else
   set Hm2=`echo $aligned_profile | awk -v m2=$m2 '{print match($1,m2)}'`
   set Hm2hit=`echo $aligned_query | cut -c$Hm2-$Hm2`
   set Dm2=`echo $aligned_profile | awk -v m2=$m2 '{print match($1,m2)+3}'`
   set Dm2hit=`echo $aligned_query | cut -c$Dm2-$Dm2`
   set Nm2=`echo $aligned_profile | awk -v m2=$m2 '{print match($1,m2)+4}'`
   set Nm2hit=`echo $aligned_query | cut -c$Nm2-$Nm2`      
   set comment2=$Hm2hit"_x_x_""$Dm2hit""_""$Nm2hit"
   endif
  set m3_check=`echo $aligned_profile | awk -v m3=$m3 '$1 ~ m3 {print $1}' | wc -l`
   if ($m3_check == "0") then
    set comment3=$m3"_region_disrupted_incoclusive" 
   else
    set Cm3=`echo $aligned_profile | awk -v m3=$m3 '{print match($1,m3)+2}'`
    set Cm3hit=`echo $aligned_query | cut -c$Cm3-$Cm3`
    set comment3=$Cm3hit
   endif
  set m4_check=`echo $aligned_profile | awk -v m4=$m4 '$1 ~ m4 {print $1}' | wc -l`
   if ($m4_check == "0") then
    set comment4=$m4"_region_disrupted_incoclusive" 
   else
    set Hm4=`echo $aligned_profile | awk -v m4=$m4 '{print match($1,m4)+1}'`
    set Hm4hit=`echo $aligned_query | cut -c$Hm4-$Hm4`
    set comment4=$Hm4hit
   endif 
  set comment5="."     
# CLASS A MOTIFS CHECK
 else if ($class == "A") then
  #motif 1
  set m1="StfK"
  #motif 2
  set m2="SDN"
  #motif 3
  set m3="Epe"
  #motif 4
  set m4="KTG"
   set m1_check=`echo $aligned_profile | awk -v m1=$m1 '$1 ~ m1 {print $1}' | wc -l`
   if ($m1_check == "0") then
    set comment1=$m1"_region_disrupted_incoclusive"
   else
    set Sminus=`echo $aligned_profile | awk -v m1=$m1 '{print match($1,m1)-1}'`
    set Sminushit=`echo $aligned_query | cut -c$Sminus-$Sminus`
    set S=`echo $aligned_profile | awk -v m1=$m1 '{print match($1,m1)}'`
    set Shit=`echo $aligned_query | cut -c$S-$S`
    set K=`echo $aligned_profile | awk -v m1=$m1 '{print match($1,m1)+3}'`
    set Khit=`echo $aligned_query | cut -c$K-$K`
    set comment1="$Sminushit""_""$Shit""_x_x_"$Khit
   endif 
  set m2_check=`echo $aligned_profile | awk -v m2=$m2 '$1 ~ m2 {print $1}' | wc -l`
   if ($m2_check == "0") then
    set comment2=$m2"_region_disrupted_incoclusive" 
   else
   set Sm2=`echo $aligned_profile | awk -v m2=$m2 '{print match($1,m2)}'`
   set Sm2hit=`echo $aligned_query | cut -c$Sm2-$Sm2`
   set Dm2=`echo $aligned_profile | awk -v m2=$m2 '{print match($1,m2)+1}'`
   set Dm2hit=`echo $aligned_query | cut -c$Dm2-$Dm2`
   set Nm2=`echo $aligned_profile | awk -v m2=$m2 '{print match($1,m2)+2}'`
   set Nm2hit=`echo $aligned_query | cut -c$Nm2-$Nm2`      
   set comment2="$Sm2hit""_""$Dm2hit""_""$Nm2hit"
   endif
  set m3_check=`echo $aligned_profile | awk -v m3=$m3 '$1 ~ m3 {print $1}' | wc -l`
   if ($m3_check == "0") then
    set comment3=$m3"_region_disrupted_incoclusive"
    set comment5="."
   else
    set Em3=`echo $aligned_profile | awk -v m3=$m3 '{print match($1,m3)}'`
    set Em3hit=`echo $aligned_query | cut -c$Em3-$Em3`
    set comment3=$Em3hit
    if ($comment3 != "E") then
     set B1=`expr $Em3 - 5`
     set B2=`expr $Em3 + 5`
     set E_vic=`echo $aligned_query | cut -c$B1-$B2 | grep E | wc -l`
     if ($E_vic != "0") then
      set comment5="E_PRESENT_IN_166_AREA(+_-5_AA)"
     else
      set comment5="E_NOT_PRESENT_IN_166_AREA(+_-5_AA)"
     endif
    else
     set comment5="."
    endif  
   endif
  set m4_check=`echo $aligned_profile | awk -v m4=$m4 '$1 ~ m4 {print $1}' | wc -l`
   if ($m4_check == "0") then
    set comment4=$m4"_region_disrupted_incoclusive" 
   else
    set Km4=`echo $aligned_profile | awk -v m4=$m4 '{print match($1,m4)}'`
    set Km4hit=`echo $aligned_query | cut -c$Km4-$Km4`
    set Tm4=`echo $aligned_profile | awk -v m4=$m4 '{print match($1,m4)+1}'`
    set Tm4hit=`echo $aligned_query | cut -c$Tm4-$Tm4`
    set Gm4=`echo $aligned_profile | awk -v m4=$m4 '{print match($1,m4)+2}'`
    set Gm4hit=`echo $aligned_query | cut -c$Gm4-$Gm4`
    set Xm4=`echo $aligned_profile | awk -v m4=$m4 '{print match($1,m4)+3}'`
    set Xm4hit=`echo $aligned_query | cut -c$Xm4-$Xm4`
    set Zm4=`echo $aligned_profile | awk -v m4=$m4 '{print match($1,m4)+4}'`
    set Zm4hit=`echo $aligned_query | cut -c$Zm4-$Zm4`            
    set comment4="$Km4hit""_""$Tm4hit""_""$Gm4hit""_""$Xm4hit""_""$Zm4hit"
   endif
 
# CLASS C MOTIFS CHECK
 else if ($class == "C") then
  #motif 1
  set m1="SvSK"
  #motif 2
  set m2="YsN"
  #motif 3
  set m3="YG"
  #motif 4
  set m4="KtG"
   set m1_check=`echo $aligned_profile | awk -v m1=$m1 '$1 ~ m1 {print $1}' | wc -l`
   if ($m1_check == "0") then
    set comment1=$m1"_region_disrupted_incoclusive"
   else
    set S=`echo $aligned_profile | awk -v m1=$m1 '{print match($1,m1)}'`
    set Shit=`echo $aligned_query | cut -c$S-$S`
    set S2=`echo $aligned_profile | awk -v m1=$m1 '{print match($1,m1)+2}'`
    set S2hit=`echo $aligned_query | cut -c$S2-$S2`
    set K=`echo $aligned_profile | awk -v m1=$m1 '{print match($1,m1)+3}'`
    set Khit=`echo $aligned_query | cut -c$K-$K`
    set comment1="$Shit""_x_""$S2hit""_""$Khit"
   endif 
  set m2_check=`echo $aligned_profile | awk -v m2=$m2 '$1 ~ m2 {print $1}' | wc -l`
   if ($m2_check == "0") then
    set comment2=$m2"_region_disrupted_incoclusive" 
   else
   set Ym2=`echo $aligned_profile | awk -v m2=$m2 '{print match($1,m2)}'`
   set Ym2hit=`echo $aligned_query | cut -c$Ym2-$Ym2`
   set sm2=`echo $aligned_profile | awk -v m2=$m2 '{print match($1,m2)+1}'`
   set sm2hit=`echo $aligned_query | cut -c$sm2-$sm2`
   set Nm2=`echo $aligned_profile | awk -v m2=$m2 '{print match($1,m2)+2}'`
   set Nm2hit=`echo $aligned_query | cut -c$Nm2-$Nm2`      
   set comment2="$Ym2hit""_""$sm2hit""_""$Nm2hit"
   endif  
  set m3_check=`echo $aligned_profile | awk -v m3=$m3 '$1 ~ m3 {print $1}' | wc -l`
   if ($m3_check == "0") then
    set comment3=$m3"_region_disrupted_incoclusive" 
   else
    set Ym3=`echo $aligned_profile | awk -v m3=$m3 '{print match($1,m3)}'`
    set Ym3hit=`echo $aligned_query | cut -c$Ym3-$Ym3`
    set comment3=$Ym3hit
   endif
  set m4_check=`echo $aligned_profile | awk -v m4=$m4 '$1 ~ m4 {print $1}' | wc -l`
   if ($m4_check == "0") then
    set comment4=$m4"_region_disrupted_incoclusive" 
   else
    set Km4=`echo $aligned_profile | awk -v m4=$m4 '{print match($1,m4)}'`
    set Km4hit=`echo $aligned_query | cut -c$Km4-$Km4`
    set Tm4=`echo $aligned_profile | awk -v m4=$m4 '{print match($1,m4)+1}'`
    set Tm4hit=`echo $aligned_query | cut -c$Tm4-$Tm4`
    set Gm4=`echo $aligned_profile | awk -v m4=$m4 '{print match($1,m4)+2}'`
    set Gm4hit=`echo $aligned_query | cut -c$Gm4-$Gm4`        
    set comment4="$Km4hit""_""$Tm4hit""_""$Gm4hit"
   endif
 # check class C/class D bifunctional enzymes
 set protein_seq=`cat "$results"_class_"$class"/"$results"_class_"$class"_wcodes.dat | awk -v ent=$j '$1 == ent {print $19}'`
 set class_d_check=`cat "$results"_class_D/"$results"_class_D_wcodes.dat | grep $protein_seq | wc -l`
 if ($class_d_check != "0") then
  set comment5="Class_C/Class_D"
 else
  set comment5="."
 endif
# CLASS D MOTIFS CHECK
 else if ($class == "D") then
  #motif 1
  set m1="STFK"
  #motif 2
  set m2="SvV"
  #motif 3
  set m3="KtG"
  #motif 4
  set m4="GW"
   set m1_check=`echo $aligned_profile | awk -v m1=$m1 '$1 ~ m1 {print $1}' | wc -l`
   if ($m1_check == "0") then
    set comment1=$m1"_region_disrupted_incoclusive"
   else
    set S=`echo $aligned_profile | awk -v m1=$m1 '{print match($1,m1)}'`
    set Shit=`echo $aligned_query | cut -c$S-$S`
    set T=`echo $aligned_profile | awk -v m1=$m1 '{print match($1,m1)+1}'`
    set Thit=`echo $aligned_query | cut -c$T-$T`
    set F=`echo $aligned_profile | awk -v m1=$m1 '{print match($1,m1)+2}'`
    set Fhit=`echo $aligned_query | cut -c$F-$F`
    set K=`echo $aligned_profile | awk -v m1=$m1 '{print match($1,m1)+3}'`
    set Khit=`echo $aligned_query | cut -c$K-$K`
    set comment1="$Shit""_""$Thit""_""$Fhit""_""$Khit"
   endif
  set m2_check=`echo $aligned_profile | awk -v m2=$m2 '$1 ~ m2 {print $1}' | wc -l`
   if ($m2_check == "0") then
    set comment2=$m2"_region_disrupted_incoclusive" 
   else
   set Sm2=`echo $aligned_profile | awk -v m2=$m2 '{print match($1,m2)}'`
   set Sm2hit=`echo $aligned_query | cut -c$Sm2-$Sm2`
   set vm2=`echo $aligned_profile | awk -v m2=$m2 '{print match($1,m2)+1}'`
   set vm2hit=`echo $aligned_query | cut -c$vm2-$vm2`
   set Vm2=`echo $aligned_profile | awk -v m2=$m2 '{print match($1,m2)+2}'`
   set Vm2hit=`echo $aligned_query | cut -c$Vm2-$Vm2`      
   set comment2="$Sm2hit""_""$vm2hit""_""$Vm2hit"
   endif  
  set m3_check=`echo $aligned_profile | awk -v m3=$m3 '$1 ~ m3 {print $1}' | wc -l`
   if ($m3_check == "0") then
    set comment3=$m3"_region_disrupted_incoclusive" 
   else
    set Km3=`echo $aligned_profile | awk -v m3=$m3 '{print match($1,m3)}'`
    set Km3hit=`echo $aligned_query | cut -c$Km3-$Km3`
    set Tm3=`echo $aligned_profile | awk -v m3=$m3 '{print match($1,m3)+1}'`
    set Tm3hit=`echo $aligned_query | cut -c$Tm3-$Tm3`
    set Gm3=`echo $aligned_profile | awk -v m3=$m3 '{print match($1,m3)+2}'`
    set Gm3hit=`echo $aligned_query | cut -c$Gm3-$Gm3`        
    set comment3="$Km3hit""_""$Tm3hit""_""$Gm3hit"
   endif 
  set m4_check=`echo $aligned_profile | awk -v m4=$m4 '$1 ~ m4 {print $1}' | wc -l`
   if ($m4_check == "0") then
    set comment4=$m4"_region_disrupted_incoclusive" 
   else
    set Gm4=`echo $aligned_profile | awk -v m4=$m4 '{print match($1,m4)}'`
    set Gm4hit=`echo $aligned_query | cut -c$Gm4-$Gm4`
    set Wm4=`echo $aligned_profile | awk -v m4=$m4 '{print match($1,m4)+1}'`
    set Wm4hit=`echo $aligned_query | cut -c$Wm4-$Wm4`    
    set comment4="$Gm4hit""_""$Wm4hit"
   endif
 # check class C/class D bifunctional enzymes
 set protein_seq=`cat "$results"_class_"$class"/"$results"_class_"$class"_wcodes.dat | awk -v ent=$j '$1 == ent {print $19}'`
 set class_c_check=`cat "$results"_class_C/"$results"_class_C_wcodes.dat | grep $protein_seq | wc -l`
  if ($class_c_check != "0") then
   set comment5="Class_C/Class_D"
  else
   set comment5="."
  endif 
# CLASS B3 MOTIFS CHECK
 else if ($class == "B3") then  
  #motif 1
  set m1="HfDH"
  #motif 2
  set m2="GHT"
  #motif 3
  set m3="aHp"
   set m1_check=`echo $aligned_profile | awk -v m1=$m1 '$1 ~ m1 {print $1}' | wc -l`
   if ($m1_check == "0") then
    set comment1=$m1"_region_disrupted_incoclusive"
   else
    set Z=`echo $aligned_profile | awk -v m1=$m1 '{print match($1,m1)-2}'`
    set Zhit=`echo $aligned_query | cut -c$Z-$Z`
    set X=`echo $aligned_profile | awk -v m1=$m1 '{print match($1,m1)-1}'`
    set Xhit=`echo $aligned_query | cut -c$X-$X`
    set H=`echo $aligned_profile | awk -v m1=$m1 '{print match($1,m1)}'`
    set Hhit=`echo $aligned_query | cut -c$H-$H`
    set Y=`echo $aligned_profile | awk -v m1=$m1 '{print match($1,m1)+1}'`
    set Yhit=`echo $aligned_query | cut -c$Y-$Y`
    set D=`echo $aligned_profile | awk -v m1=$m1 '{print match($1,m1)+2}'`
    set Dhit=`echo $aligned_query | cut -c$D-$D`
    set H2=`echo $aligned_profile | awk -v m1=$m1 '{print match($1,m1)+3}'`
    set H2hit=`echo $aligned_query | cut -c$H2-$H2`
    set comment1="$Zhit""_""$Xhit""_""$Hhit""_""$Yhit""_""$Dhit""_""$H2hit"
   endif
  set m2_check=`echo $aligned_profile | awk -v m2=$m2 '$1 ~ m2 {print $1}' | wc -l`
   if ($m2_check == "0") then
    set comment2=$m2"_region_disrupted_incoclusive" 
   else
   set Hm2=`echo $aligned_profile | awk -v m2=$m2 '{print match($1,m2)+1}'`
   set Hm2hit=`echo $aligned_query | cut -c$Hm2-$Hm2`    
   set comment2="$Hm2hit"
   endif  
  set m3_check=`echo $aligned_profile | awk -v m3=$m3 '$1 ~ m3 {print $1}' | wc -l`
   if ($m3_check == "0") then
    set comment3=$m3"_region_disrupted_incoclusive" 
   else
   set Hm3=`echo $aligned_profile | awk -v m3=$m3 '{print match($1,m3)+1}'`
   set Hm3hit=`echo $aligned_query | cut -c$Hm3-$Hm3`    
   set comment3="$Hm3hit"
   endif 
  set comment4="."
  set comment5="."
# CLASS B2 MOTIFS CHECK
 else if ($class == "B2") then  
  #motif 1
  set m1="nyhtd"
  #motif 2
  set m2="gnc"
  #motif 3
  set m3="dhy"
  set m1_check=`echo $aligned_profile | awk -v m1=$m1 '$1 ~ m1 {print $1}' | wc -l`
   if ($m1_check == "0") then
    set comment1=$m1"_region_disrupted_incoclusive"
   else
    set N1=`echo $aligned_profile | awk -v m1=$m1 '{print match($1,m1)}'`
    set N1hit=`echo $aligned_query | cut -c$N1-$N1`
    set Y1=`echo $aligned_profile | awk -v m1=$m1 '{print match($1,m1)+1}'`
    set Y1hit=`echo $aligned_query | cut -c$Y1-$Y1`
    set H1=`echo $aligned_profile | awk -v m1=$m1 '{print match($1,m1)+2}'`
    set H1hit=`echo $aligned_query | cut -c$H1-$H1`
    set T1=`echo $aligned_profile | awk -v m1=$m1 '{print match($1,m1)+3}'`
    set T1hit=`echo $aligned_query | cut -c$T1-$T1`
    set D1=`echo $aligned_profile | awk -v m1=$m1 '{print match($1,m1)+4}'`
    set D1hit=`echo $aligned_query | cut -c$D1-$D1`
    set comment1=$N1hit"_"$Y1hit"_""$H1hit""_""$T1hit""_""$D1hit"
   endif
  set m2_check=`echo $aligned_profile | awk -v m2=$m2 '$1 ~ m2 {print $1}' | wc -l`
   if ($m2_check == "0") then
    set comment2=$m2"_region_disrupted_incoclusive" 
   else
   set Cm2=`echo $aligned_profile | awk -v m2=$m2 '{print match($1,m2)+2}'`
   set Cm2hit=`echo $aligned_query | cut -c$Cm2-$Cm2`
   set comment2="$Cm2hit"
  endif
  set m3_check=`echo $aligned_profile | awk -v m3=$m3 '$1 ~ m3 {print $1}' | wc -l`
   if ($m3_check == "0") then
    set comment3=$m3"_region_disrupted_incoclusive" 
   else
   set Hm3=`echo $aligned_profile | awk -v m3=$m3 '{print match($1,m3)+1}'`
   set Hm3hit=`echo $aligned_query | cut -c$Hm3-$Hm3`    
   set comment3="$Hm3hit"
   endif
  set comment4="."
  set comment5="."   
endif
cat "$results"_class_"$class"/"$results"_class_"$class"_wcodes_plen.dat | awk -v ent=$j -v clstrc=$cluster_cds -v arb=$ref_db_compar -v clstrp=$cluster_proteins -v nspt=$gr_neg_spt  -v nsps=$gr_neg_site -v pspt=$gr_pos_spt -v psps=$gr_pos_site -v cm1=$comment1 -v cm2=$comment2 -v cm3=$comment3 -v cm4=$comment4 -v cm5="$comment5" '$1 == ent {print $0,clstrc,arb,clstrp,nspt,nsps,pspt,psps,cm1,cm2,cm3,cm4,cm5}' >> "$results"_class_"$class"/"$class"_whole_proteins_positives_meta.dat
end


