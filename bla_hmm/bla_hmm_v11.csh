#!/usr/bin/csh
# SCRIPT FOR HMM SEARCHES OF BETA-LACTAMASES IN BACTERIAL GENOMES. REQUIRED SOFTWARE: PYTHON, HMMER V3, PROKKA, BLAST+ 
# STATHIS KOTSAKIS. LABORATORY OF BACTERIOLOGY. HPI. 2019

set BLA_HMM=$HOME/bla_hmm
set s=$1
set classes=$2
set strain=$3
set outdir=$4
set frames=$5
# Usage: csh bla_hmm.csh <sequence in fasta> <bla class>[all]  <strain>[strain] <output_directory>[blas_hmm_results] <translation frames>[6] <print alignemnts>[no]
# Example: csh bla_hmm.csh contigs.fasta all strain_name blas_hmm_example 6 yes
# Example w defaults: csh bla_hmm.csh contigs.fasta all

# set defaults

if ($classes == "all") then
 set classes=(A C D B1 B2 B3)
endif
set nonomatch

if ($outdir == "") then
 set outdir=blas_hmm_results
endif

rm -rf $outdir
mkdir $outdir

#annotation check
set pattern_gff = ".gff"
set filetype_gff =  ( *$pattern_gff* ) 
set pattern_faa = ".faa"
set filetype_faa =  ( *$pattern_faa* ) 

if ( -f $filetype_gff ) then
 set gff=1
else
 set gff=0
endif

if ( -f $filetype_faa ) then
 set faa=1
else
 set faa=0
endif

set annotation_c=`echo $gff $faa | awk '{print $1+$2}'`
if ($annotation_c == "2") then
 set annotation="yes"
 set gff_file=$filetype_gff
 set faa_file=$filetype_faa
 mkdir $outdir/annotation
 cp $gff_file $outdir/annotation/
 cp $faa_file $outdir/annotation/
else
 set annotation="no"
endif

if ($strain == "") then
 set strain=strain
endif



if ($frames == "") then
 set frames=6
endif

# extract cut-off scores and profile lengths
foreach j ($classes)
 set cut_off_$j=`sort -k2 -n $BLA_HMM/hmm_models/class_"$j"/scores.dat | head -1 | awk '{print int($2)}'`
 set neg_cut_off_$j=`cat $BLA_HMM/hmm_models/class_"$j"/negative.dat | awk '{print int($2)}'`
 set prof_len_$j=`cat $BLA_HMM/hmm_models/class_"$j"/class_"$j"_model.hmm | awk 'NR==3 {print $2}'`
end

mkdir $outdir/translations
mkdir $outdir/positives
mkdir $outdir/greys
mkdir $outdir/negatives

cat $s | awk '{ sub("\r$", ""); print }' | sed 's/)/_/g' | sed 's/(/_/g' > $outdir/sequence_temp_fas
set genomesize=`cat $outdir/sequence_temp_fas | sed '/>/d' | tr -d '\n' | wc -c`

echo "-------------------------AGCT SEQUENCE PROCESSOR HAS STARTED------------------------------"

timeout 1800 $BLA_HMM/nuc_translator/./*.py $outdir/sequence_temp_fas $outdir/translations/trans DNA $frames 1 NOBIN 24 > -
if ($? == 124) then
echo "----------------------------Translation Not Finished--------------------------------------"
echo "-------------------------------DNA Input Problem------------------------------------------"
echo $strain "DNA_Input_Problem" > $outdir/greys/"$strain"_greys_total.dat
echo $strain "DNA_Input_Problem" > $outdir/positives/"$strain"_positives_total.dat
echo $strain "DNA_Input_Problem" > $outdir/positives/"$strain"_pos_nums.dat
echo $strain "DNA_Input_Problem" > $outdir/greys/"$strain"_greys_nums.dat
rm -
rm -rf $outdir/translations
exit 0
else
rm -
echo "----------------------------Translation Has Finished--------------------------------------"
mkdir $outdir/hmm_hits
echo "---------------HMMER SCORE CUT-OFFS FOR EACH BETA-LACTAMASE CLASS-------------------------"
echo "-------------SEE hmm_models/model_validation.pdf FOR JUSTIFICATION------------------------"
echo "++++++++++++++CLASS-A: POSITIVE>=""$cut_off_A"", NEGATIVE<""$neg_cut_off_A""+++++++++"
echo "++++++++++++++CLASS-C: POSITIVE>=""$cut_off_C"", NEGATIVE<""$neg_cut_off_C""+++++++++"
echo "++++++++++++++CLASS-D: POSITIVE>=""$cut_off_D"", NEGATIVE<""$neg_cut_off_D""+++++++++"
echo "++++++++++++++CLASS-B1: POSITIVE>=""$cut_off_B1"", NEGATIVE<""$neg_cut_off_B1""++++++++"
echo "++++++++++++++CLASS-B2: POSITIVE>=""$cut_off_B2"", NEGATIVE<""$neg_cut_off_B2""++++++++"
echo "++++++++++++++CLASS-B3: POSITIVE>=""$cut_off_B3"", NEGATIVE<""$neg_cut_off_B3""++++++++"
echo "------------------------------------------------------------------------------------------"
set i=1
while ( $i <= $frames)
 mkdir $outdir/hmm_hits/frame_"$i"
 cp $outdir/translations/trans.tr_frame"$i" $outdir/hmm_hits/frame_"$i"/trans.tr_frame"$i"
 foreach k ($classes)
  hmmsearch $BLA_HMM/hmm_models/class_"$k"/class_"$k"_model.hmm $outdir/translations/trans.tr_frame"$i" > $outdir/hmm_hits/frame_"$i"/class_"$k"_frame_"$i".dat
  cat $outdir/hmm_hits/frame_"$i"/class_"$k"_frame_"$i".dat | sed -n '/Scores for complete sequences/,/Domain annotation for each/p' | sed '/------ inclusion threshold ------/d' | sed '1,4d' | head -n -3 > $outdir/hmm_hits/frame_"$i"/class_"$k"_frame_"$i".tab
  set hit_check=`cat $outdir/hmm_hits/frame_"$i"/class_"$k"_frame_"$i".tab | awk '$0 ~ "No hits detected that satisfy reporting thresholds" {print $0}' | wc -l`
  if ($hit_check == 0) then
   echo ">" >> $outdir/sequence_temp_fas
   set comment="TRUE_HIT"
   set cut_off=`eval echo \$cut_off_$k`
   set neg_cut_off=`eval echo \$neg_cut_off_$k`
   set prof_len=`eval echo \$prof_len_$k`
   set hits_nodes=`cat $outdir/hmm_hits/frame_"$i"/class_"$k"_frame_"$i".tab | awk '{print $9}'`
   cat $outdir/hmm_hits/frame_"$i"/class_"$k"_frame_"$i".dat | sed -n '/Domain annotation for each/,/Internal pipeline statistics summary/p' | sed '1,1d' | head -n -5 > $outdir/hmm_hits/frame_"$i"/class_"$k"_frame_"$i".ali.temp
   echo ">>" >> $outdir/hmm_hits/frame_"$i"/class_"$k"_frame_"$i".ali.temp
   foreach m ($hits_nodes)
    set contig_hit=`echo $m |  sed 's/_fr*.//g'`
    cat $outdir/hmm_hits/frame_"$i"/class_"$k"_frame_"$i".ali.temp | sed -n '/>> '$m'/,/>>/p' | sed '$ d' | sed -n '/>> '$m'/,/Alignments for each/p' | sed '1,3d' | head -n -2 | sort -k 3,3nr > $outdir/hmm_hits/frame_"$i"/"$m"_class_"$k"_frame_"$i"_hits.dat
    set hits_start=`cat $outdir/hmm_hits/frame_"$i"/"$m"_class_"$k"_frame_"$i"_hits.dat | awk '{print $10}' | sort -u`
    foreach v ($hits_start)
     set hit_score=`cat $outdir/hmm_hits/frame_"$i"/"$m"_class_"$k"_frame_"$i"_hits.dat | awk -v hit=$v '$10 == hit {print $3}'`
     set hit_check=`cat $outdir/hmm_hits/frame_"$i"/"$m"_class_"$k"_frame_"$i"_hits.dat | awk -v hit=$v '$10 == hit {print int($3)}'`
     set hit_start=$v
     set hit_end=`cat $outdir/hmm_hits/frame_"$i"/"$m"_class_"$k"_frame_"$i"_hits.dat | awk -v hit=$v '$10 == hit {print $11}'`
     set hit_domain=`cat $outdir/hmm_hits/frame_"$i"/"$m"_class_"$k"_frame_"$i"_hits.dat | awk -v hit=$v '$10 == hit {print $1}'`
     set hit_hmmfrom=`cat $outdir/hmm_hits/frame_"$i"/"$m"_class_"$k"_frame_"$i"_hits.dat | awk -v hit=$v '$10 == hit {print $7}'`
     set hit_hmmto=`cat $outdir/hmm_hits/frame_"$i"/"$m"_class_"$k"_frame_"$i"_hits.dat | awk -v hit=$v '$10 == hit {print $8}'`
     if ($hit_check < $neg_cut_off) then
      echo $m class_$k negative hit start $hit_start end $hit_end score $hit_score >> $outdir/negatives/"$strain"_negatives.dat
     else if ($hit_check < $cut_off) then # greys
      echo $m class_$k grey_zone hit start $hit_start end $hit_end score $hit_score >> $outdir/greys/"$strain"_greys.dat
      sed -n '/>'$contig_hit'/,/>/p' $outdir/sequence_temp_fas > $outdir/greys/contig_"$contig_hit"_class_"$k"_grey.seq
      sed -i '/>/d' $outdir/greys/contig_"$contig_hit"_class_"$k"_grey.seq
      set locus_grey=`echo $contig_hit | sed 's/\..*//g' | sed 's/.*_//g'`
      # SEQUENCE EXTRACTION
      # cds       
      if ($i < 4) then
       set cds_start_g=`echo $hit_start | awk -v i=$i '{print 3*$1+i-1}'`
       set cds_end_g=`echo $hit_end | awk -v i=$i '{print 3*$1+i-1}'`
       echo $strain $genomesize $contig_hit class_"$k" $hit_score + $cds_start_g $cds_end_g > $outdir/greys/"$contig_hit"_class_"$k"_"$cds_start_g"_"$cds_end_g"_grey_hits.dat
       set comp="plus"
      else
       set cds_start_c_g=`echo $hit_start | awk -v i=$i '{print 3*$1+i-4}'`
       set cds_end_c_g=`echo $hit_end | awk -v i=$i '{print 3*$1+i-4}'`
       set contig_len_g=`cat $outdir/greys/contig_"$contig_hit"_class_"$k"_grey.seq | tr -d '\n' | wc -c`
       set cds_start_g=`echo $contig_len_g $cds_end_c_g | awk '{print $1-$2}'`
       set cds_end_g=`echo $contig_len_g $cds_start_c_g | awk '{print $1-$2}'`
       echo $strain $genomesize $contig_hit class_"$k" $hit_score - $cds_start_g $cds_end_g > $outdir/greys/"$contig_hit"_class_"$k"_"$cds_start_g"_"$cds_end_g"_grey_hits.dat
       set comp="minus"
      endif
       echo ">" >> $outdir/translations/trans.tr_frame"$i"
       sed -n '/>'$m'/,/>/p' $outdir/translations/trans.tr_frame"$i" > $outdir/translations/trans.tr_frame"$i"_"$m"_"$k".temp
       sed -i '/>/d' $outdir/translations/trans.tr_frame"$i"_"$m"_"$k".temp
       set pfrag_seq_trimmed_g=`cat $outdir/translations/trans.tr_frame"$i"_"$m"_"$k".temp | tr -d '\n' | cut -c$hit_start-$hit_end | sed 's/[*]/_/g' | cut -c15- | rev | cut -c15- | rev`
       sed -i '$ d' $outdir/translations/trans.tr_frame"$i"
      csh $BLA_HMM/cds_extractor_v2.csh $contig_hit $outdir/sequence_temp_fas $cds_start_g $cds_end_g $comp $outdir/greys/contig_"$contig_hit"_class_"$k"_grey.seq $pfrag_seq_trimmed_g $annotation 900 "$contig_hit"_class_"$k"_"$cds_start_g"_"$cds_end_g"_genbank $outdir
      cp -r $outdir/cds_extraction/* $outdir/greys/

       # knowledge based evaluation of hits
       # check if a B3 hit is a HARLDQ motif MBL
       if ( $k == "B3" ) then
        cat $outdir/greys/"$contig_hit"_class_"$k"_"$cds_start_g"_"$cds_end_g"_genbank_wseqs.dat | awk '{print ">"$4"\n" $6}' | sed 's/_STOP_/\*/g' > $outdir/greys/B3_"$contig_hit"_"$cds_start_g"_grey.temp
        hmmsearch $BLA_HMM/hmm_models/class_HARLDQ/class_HARLDQ_model.hmm $outdir/greys/B3_"$contig_hit"_"$cds_start_g"_grey.temp > $outdir/greys/B3_"$contig_hit"_"$cds_start_g"_grey_vs_HARLDQ.dat.temp
        set h_greycheck=`cat $outdir/greys/B3_"$contig_hit"_"$cds_start_g"_grey_vs_HARLDQ.dat.temp | sed -n '/Scores for complete sequences/,/Domain annotation for each/p' | sed '/------ inclusion threshold ------/d' | sed '1,4d' | head -n -3 | awk '$0 ~ "No hits detected that satisfy reporting thresholds" {print $0}' | wc -l`
        if ($h_greycheck == 0) then
         cat $outdir/greys/B3_"$contig_hit"_"$cds_start_g"_grey_vs_HARLDQ.dat.temp | sed -n '/Domain annotation for each/,/Internal pipeline statistics summary/p' | sed '1,1d' | head -n -5 > $outdir/greys/B3_"$contig_hit"_"$cds_start_g"_grey_vs_HARLDQ.ali.temp
         cat $outdir/greys/B3_"$contig_hit"_"$cds_start_g"_grey_vs_HARLDQ.ali.temp | sed -n '/>>/,/Alignments for each/p' | sed '1,3d' | head -n -2 | sort -k 3,3nr > $outdir/greys/B3_"$contig_hit"_"$cds_start_g"_grey_vs_HARLDQ.hits.temp
         set grey_co_h=`cat $outdir/greys/B3_"$contig_hit"_"$cds_start_g"_grey_vs_HARLDQ.hits.temp | awk '$3 >= 415 {print $3}' | wc -l`
         if ( $grey_co_h != 0 ) then
          set comment="HARLDQ_MOTIF"
         else
          set comment="TRUE_HIT"
         endif
        endif
        rm $outdir/greys/*.temp
       endif
       # check if a D hit is a BlaR1
       if ( $k == "D" ) then
        cat $outdir/greys/"$contig_hit"_class_"$k"_"$cds_start_g"_"$cds_end_g"_genbank_wseqs.dat | awk '{print ">"$4"\n" $6}' | sed 's/_STOP_/\*/g' > $outdir/greys/D_"$contig_hit"_"$cds_start_g"_grey.temp
        hmmsearch $BLA_HMM/hmm_models/class_BlaR1/class_BlaR1_model.hmm $outdir/greys/D_"$contig_hit"_"$cds_start_g"_grey.temp > $outdir/greys/D_"$contig_hit"_"$cds_start_g"_grey_vs_BlaR1.dat.temp
        set b_greycheck=`cat $outdir/greys/D_"$contig_hit"_"$cds_start_g"_grey_vs_BlaR1.dat.temp | sed -n '/Scores for complete sequences/,/Domain annotation for each/p' | sed '/------ inclusion threshold ------/d' | sed '1,4d' | head -n -3 | awk '$0 ~ "No hits detected that satisfy reporting thresholds" {print $0}' | wc -l`
        if ($b_greycheck == 0) then
         cat $outdir/greys/D_"$contig_hit"_"$cds_start_g"_grey_vs_BlaR1.dat.temp | sed -n '/Domain annotation for each/,/Internal pipeline statistics summary/p' | sed '1,1d' | head -n -5 > $outdir/greys/D_"$contig_hit"_"$cds_start_g"_grey_vs_BlaR1.ali.temp
         cat $outdir/greys/D_"$contig_hit"_"$cds_start_g"_grey_vs_BlaR1.ali.temp | sed -n '/>>/,/Alignments for each/p' | sed '1,3d' | head -n -2 | sort -k 3,3nr > $outdir/greys/D_"$contig_hit"_"$cds_start_g"_grey_vs_BlaR1.hits.temp
         set grey_co_b=`cat $outdir/greys/D_"$contig_hit"_"$cds_start_g"_grey_vs_BlaR1.hits.temp | awk '$3 >= 260 && $7 < 300 {print $3}' | wc -l`
         if ( $grey_co_b != 0 ) then
          set comment="BlaR1"
         else
          set comment="TRUE_HIT"
         endif
        endif
        rm $outdir/greys/*.temp
       endif                
       # run B2 model for B1 greys
       if ($k == "B1") then
        cat $outdir/greys/"$contig_hit"_class_"$k"_"$cds_start_g"_"$cds_end_g"_genbank_wseqs.dat | awk '{print ">"$4"\n" $6}' > $outdir/greys/B1_"$contig_hit"_"$cds_start_g"_grey.temp
        hmmsearch $BLA_HMM/hmm_models/class_B2/class_B2_model.hmm $outdir/greys/B1_"$contig_hit"_"$cds_start_g"_grey.temp > $outdir/greys/B1_"$contig_hit"_"$cds_start_g"_grey_vs_B2.dat.temp
        set b2_hitcheck=`cat $outdir/greys/B1_"$contig_hit"_"$cds_start_g"_grey_vs_B2.dat.temp | sed -n '/Scores for complete sequences/,/Domain annotation for each/p' | sed '/------ inclusion threshold ------/d' | sed '1,4d' | head -n -3 | awk '$0 ~ "No hits detected that satisfy reporting thresholds" {print $0}' | wc -l`
        if ($b2_hitcheck == 0) then
         cat $outdir/greys/B1_"$contig_hit"_"$cds_start_g"_grey_vs_B2.dat.temp | sed -n '/Domain annotation for each/,/Internal pipeline statistics summary/p' | sed '1,1d' | head -n -5 > $outdir/greys/B1_"$contig_hit"_"$cds_start_g"_grey_vs_B2.ali.temp
         cat $outdir/greys/B1_"$contig_hit"_"$cds_start_g"_grey_vs_B2.ali.temp | sed -n '/>>/,/Alignments for each/p' | sed '1,3d' | head -n -2 | sort -k 3,3nr > $outdir/greys/B1_"$contig_hit"_"$cds_start_g"_grey_vs_B2.hits.temp
         set hits_co_b2=`cat $outdir/greys/B1_"$contig_hit"_"$cds_start_g"_grey_vs_B2.hits.temp | awk -v c=$cut_off_B2 '$3 >= c {print $3}' | wc -l`
         if ( $hits_co_b2 != 0 ) then
          set comment="B2_POSITIVE"
         else
          set comment="TRUE_HIT"
         endif
        endif
        rm $outdir/greys/*.temp
       endif
       # alignment extraction
       cat $outdir/hmm_hits/frame_"$i"/class_"$k"_frame_"$i".ali.temp | sed -n '/>> '$m'/,/>>/p' > $outdir/greys/"$m"_class_"$k"_frame_"$i"_grey.ali.temp
       echo "==" >> $outdir/greys/"$m"_class_"$k"_frame_"$i"_grey.ali.temp
       cat $outdir/greys/"$m"_class_"$k"_frame_"$i"_grey.ali.temp | sed -n '/== domain '$hit_domain' /,/==/p' | sed '/==/d' | sed '/>>/d' > $outdir/greys/"$m"_class_"$k"_frame_"$i"_"$cds_start_g"_grey.ali
       set aligned_profile=`cat $outdir/greys/"$m"_class_"$k"_frame_"$i"_"$cds_start_g"_grey.ali | sed '/+/d' | awk '{print $3}' | awk 'NF > 0' | awk 'ORS=NR%2?FS:RS' | awk '{print $1}' | tr -d '\n'`
       set aligned_query=`cat $outdir/greys/"$m"_class_"$k"_frame_"$i"_"$cds_start_g"_grey.ali | sed '/+/d' | awk '{print $3}' | awk 'NF > 0' | awk 'ORS=NR%2?FS:RS' | awk '{print $2}' | tr -d '\n' | sed 's/[*]/_/g'`
       set aligned_perc=`echo $hit_hmmfrom $hit_hmmto | awk -v pf=$prof_len '{print ($2-$1)/pf*100}'`
       set insertions=`echo $aligned_profile | tr -cd '.' | wc -c`
       set deletions=`echo $aligned_query | tr -cd '-' | wc -c`
       set stops=`echo $aligned_query | tr -cd '_' | wc -c`
       # check if an A hit contains E166
       if ( $k == "A" ) then
        set Epel=`echo $aligned_profile | awk '$1 ~ "EpeL" {print $1}' | wc -l`
        if ($Epel == "0") then
         set comment="E166_region_missing_incoclusive"
        else
         set Echar=`echo $aligned_profile | awk '{print match($1,"EpeL")}'`
         set Ehit=`echo $aligned_query | cut -c$Echar-$Echar`
          if ($Ehit != "E") then
           set comment="E166_is_"$Ehit"_Likely_PBP"
         else
          set comment="TRUE_HIT"
          endif
         endif
        endif  
       # check if a B1 true grey hit contains C221
       if ( $k == "B1" ) then
        if ( $comment != "B2_POSITIVE" ) then
        set GGC=`echo $aligned_profile | awk '$1 ~ "GGC" {print $1}' | wc -l`
        if ($GGC == "0") then
         set comment="C221_region_missing_or_disrupted"
        else
         set Cchar=`echo $aligned_profile | awk '{print match($1,"GGC")+2}'`
         set Chit=`echo $aligned_query | cut -c$Cchar-$Cchar`
         if ($Chit != "C") then
          if ($Chit == "E" || $Chit == "D") then
           set comment="C221_is_"$Chit"_Likely_non_Bla_MBL"
          else
           set comment="C221_is_"$Chit"_Intersting_MBL"
          endif
          else
          set comment="TRUE_HIT"
         endif
        endif
       endif
       endif     
        echo `cat $outdir/greys/"$contig_hit"_class_"$k"_"$cds_start_g"_"$cds_end_g"_grey_hits.dat` `cat $outdir/greys/"$contig_hit"_class_"$k"_"$cds_start_g"_"$cds_end_g"_genbank_wseqs.dat` | awk -v com=$comment -v d=$outdir -v ap=$aligned_profile -v aq=$aligned_query -v ape=$aligned_perc -v ins=$insertions -v del=$deletions -v sto=$stops -v from=$hit_hmmfrom -v to=$hit_hmmto '{print $0,com,d,ap,aq,from,to,ape,ins,del,sto}' >> $outdir/greys/"$strain"_greys_total.dat
        rm -rf $outdir/cds_extraction
        rm $outdir/greys/*.temp
        # end of sequence extraction
      else #positives
       sed -n '/>'$contig_hit'/,/>/p' $outdir/sequence_temp_fas > $outdir/positives/contig_"$contig_hit"_class_"$k"_pos.seq
       sed -i '/>/d' $outdir/positives/contig_"$contig_hit"_class_"$k"_pos.seq
       set locus_pos=`echo $contig_hit | sed 's/\..*//g' | sed 's/.*_//g'`
       # SEQUENCE EXTRACTION
       # cds       
       if ($i < 4) then
        set cds_start=`echo $hit_start | awk -v i=$i '{print 3*$1+i-1}'`
        set cds_end=`echo $hit_end | awk -v i=$i '{print 3*$1+i-1}'`
        echo $strain $genomesize $contig_hit class_"$k" $hit_score + $cds_start $cds_end >> $outdir/positives/"$contig_hit"_class_"$k"_"$cds_start"_"$cds_end"_positive_hits.dat
        set comp="plus"
       else
        set cds_start_c=`echo $hit_start | awk -v i=$i '{print 3*$1+i-4}'`
        set cds_end_c=`echo $hit_end | awk -v i=$i '{print 3*$1+i-4}'` 
        set contig_len=`cat $outdir/positives/contig_"$contig_hit"_class_"$k"_pos.seq | tr -d '\n' | wc -c`
        set cds_start=`echo $contig_len $cds_end_c | awk '{print $1-$2}'`
        set cds_end=`echo $contig_len $cds_start_c | awk '{print $1-$2}'`
        echo $strain $genomesize $contig_hit class_"$k" $hit_score - $cds_start $cds_end > $outdir/positives/"$contig_hit"_class_"$k"_"$cds_start"_"$cds_end"_positive_hits.dat
        set comp="minus"
       endif
        echo ">" >> $outdir/translations/trans.tr_frame"$i"
        sed -n '/>'$m'/,/>/p' $outdir/translations/trans.tr_frame"$i" > $outdir/translations/trans.tr_frame"$i"_"$m"_"$k".temp
        sed -i '/>/d' $outdir/translations/trans.tr_frame"$i"_"$m"_"$k".temp
        set pfrag_seq_trimmed=`cat $outdir/translations/trans.tr_frame"$i"_"$m"_"$k".temp | tr -d '\n' | cut -c$hit_start-$hit_end | sed 's/[*]/_/g' | cut -c15- | rev | cut -c15- | rev`
        sed -i '$ d' $outdir/translations/trans.tr_frame"$i"
      csh $BLA_HMM/cds_extractor_v2.csh $contig_hit $outdir/sequence_temp_fas $cds_start $cds_end $comp $outdir/positives/contig_"$contig_hit"_class_"$k"_pos.seq $pfrag_seq_trimmed $annotation 900 "$contig_hit"_class_"$k"_"$cds_start"_"$cds_end"_genbank $outdir
       cp -r $outdir/cds_extraction/* $outdir/positives/
       # check if a B3 hit is a HARLDQ motif MBL
       if ( $k == "B3" ) then
        cat $outdir/positives/"$contig_hit"_class_"$k"_"$cds_start"_"$cds_end"_genbank_wseqs.dat | awk '{print ">"$4"\n" $6}' | sed 's/_STOP_/\*/g' > $outdir/positives/B3_"$contig_hit"_"$cds_start"_positive.temp
        hmmsearch $BLA_HMM/hmm_models/class_HARLDQ/class_HARLDQ_model.hmm $outdir/positives/B3_"$contig_hit"_"$cds_start"_positive.temp > $outdir/positives/B3_"$contig_hit"_"$cds_start"_positive_vs_HARLDQ.dat.temp
        set h_hitcheck=`cat $outdir/positives/B3_"$contig_hit"_"$cds_start"_positive_vs_HARLDQ.dat.temp | sed -n '/Scores for complete sequences/,/Domain annotation for each/p' | sed '/------ inclusion threshold ------/d' | sed '1,4d' | head -n -3 | awk '$0 ~ "No hits detected that satisfy reporting thresholds" {print $0}' | wc -l`
        if ($h_hitcheck == 0) then
         cat $outdir/positives/B3_"$contig_hit"_"$cds_start"_positive_vs_HARLDQ.dat.temp | sed -n '/Domain annotation for each/,/Internal pipeline statistics summary/p' | sed '1,1d' | head -n -5 > $outdir/positives/B3_"$contig_hit"_"$cds_start"_positive_vs_HARLDQ.ali.temp
         cat $outdir/positives/B3_"$contig_hit"_"$cds_start"_positive_vs_HARLDQ.ali.temp | sed -n '/>>/,/Alignments for each/p' | sed '1,3d' | head -n -2 | sort -k 3,3nr > $outdir/positives/B3_"$contig_hit"_"$cds_start"_positive_vs_HARLDQ.hits.temp
         set hits_co_h=`cat $outdir/positives/B3_"$contig_hit"_"$cds_start"_positive_vs_HARLDQ.hits.temp | awk '$3 >= 415 {print $3}' | wc -l`
         if ( $hits_co_h != 0 ) then
          set comment="HARLDQ_MOTIF"
         else
          set comment="TRUE_HIT"
         endif
        endif
        rm $outdir/positives/*.temp
       endif
       # check if a D hit is a BlaR1
       if ( $k == "D" ) then
        cat $outdir/positives/"$contig_hit"_class_"$k"_"$cds_start"_"$cds_end"_genbank_wseqs.dat | awk '{print ">"$4"\n" $6}' | sed 's/_STOP_/\*/g' > $outdir/positives/D_"$contig_hit"_"$cds_start"_positive.temp
        hmmsearch $BLA_HMM/hmm_models/class_BlaR1/class_BlaR1_model.hmm $outdir/positives/D_"$contig_hit"_"$cds_start"_positive.temp > $outdir/positives/D_"$contig_hit"_"$cds_start"_positive_vs_BlaR1.dat.temp
        set r1_hitcheck=`cat $outdir/positives/D_"$contig_hit"_"$cds_start"_positive_vs_BlaR1.dat.temp | sed -n '/Scores for complete sequences/,/Domain annotation for each/p' | sed '/------ inclusion threshold ------/d' | sed '1,4d' | head -n -3 | awk '$0 ~ "No hits detected that satisfy reporting thresholds" {print $0}' | wc -l`
        if ($r1_hitcheck == 0) then
         cat $outdir/positives/D_"$contig_hit"_"$cds_start"_positive_vs_BlaR1.dat.temp | sed -n '/Domain annotation for each/,/Internal pipeline statistics summary/p' | sed '1,1d' | head -n -5 > $outdir/positives/D_"$contig_hit"_"$cds_start"_positive_vs_BlaR1.ali.temp
         cat $outdir/positives/D_"$contig_hit"_"$cds_start"_positive_vs_BlaR1.ali.temp | sed -n '/>>/,/Alignments for each/p' | sed '1,3d' | head -n -2 | sort -k 3,3nr > $outdir/positives/D_"$contig_hit"_"$cds_start"_positive_vs_BlaR1.hits.temp
         set hits_co_r1=`cat $outdir/positives/D_"$contig_hit"_"$cds_start"_positive_vs_BlaR1.hits.temp | awk '$3 >= 260 && $7 < 300 {print $3}' | wc -l`
         if ( $hits_co_r1 != 0 ) then
          set comment="BlaR1"
         else
          set comment="TRUE_HIT"         
         endif
        endif
        rm $outdir/positives/*.temp
       endif
       # alignment extraction
       cat $outdir/hmm_hits/frame_"$i"/class_"$k"_frame_"$i".ali.temp | sed -n '/>> '$m'/,/>>/p' > $outdir/positives/"$m"_class_"$k"_frame_"$i"_positive.ali.temp
       echo "==" >> $outdir/positives/"$m"_class_"$k"_frame_"$i"_positive.ali.temp
       cat $outdir/positives/"$m"_class_"$k"_frame_"$i"_positive.ali.temp | sed -n '/== domain '$hit_domain' /,/==/p' | sed '/==/d'  | sed '/>>/d' > $outdir/positives/"$m"_class_"$k"_frame_"$i"_"$cds_start"_positive.ali
       set aligned_profile=`cat $outdir/positives/"$m"_class_"$k"_frame_"$i"_"$cds_start"_positive.ali | sed '/+/d' | awk '{print $3}' | awk 'NF > 0' | awk 'ORS=NR%2?FS:RS' | awk '{print $1}' | tr -d '\n'`
       set aligned_query=`cat $outdir/positives/"$m"_class_"$k"_frame_"$i"_"$cds_start"_positive.ali | sed '/+/d' | awk '{print $3}' | awk 'NF > 0' | awk 'ORS=NR%2?FS:RS' | awk '{print $2}' | tr -d '\n' | sed 's/[*]/_/g'`
       set aligned_perc=`echo $hit_hmmfrom $hit_hmmto | awk -v pf=$prof_len '{print ($2-$1)/pf*100}'`
       set insertions=`echo $aligned_profile | tr -cd '.' | wc -c`
       set deletions=`echo $aligned_query | tr -cd '-' | wc -c`
       set stops=`echo $aligned_query | tr -cd '_' | wc -c`
       # check if an A hit contains E166
       if ( $k == "A" ) then
        set Epel=`echo $aligned_profile | awk '$1 ~ "EpeL" {print $1}' | wc -l`
        if ($Epel == "0") then
         set comment="E166_region_missing_or_disrupted"
        else
         set Echar=`echo $aligned_profile | awk '{print match($1,"EpeL")}'`
         set Ehit=`echo $aligned_query | cut -c$Echar-$Echar`
         if ($Ehit != "E") then
          set comment="E166_is_"$Ehit"_Likely_PBP"
         else
          set comment="TRUE_HIT"
         endif
        endif
       endif
       # check if a B1 hit contains C221
       if ( $k == "B1" ) then
        set GGC=`echo $aligned_profile | awk '$1 ~ "GGC" {print $1}' | wc -l`
        if ($GGC == "0") then
         set comment="C221_region_missing_or_disrupted"
        else
         set Cchar=`echo $aligned_profile | awk '{print match($1,"GGC")+2}'`
         set Chit=`echo $aligned_query | cut -c$Cchar-$Cchar`
         if ($Chit != "C") then
          if ($Chit == "E" || $Chit == "D") then
           set comment="C221_is_"$Chit"_Likely_non_Bla_MBL"
          else
           set comment="C221_is_"$Chit"_Intersting_MBL"
          endif
          else
          set comment="TRUE_HIT"
         endif
        endif
       endif
       #rm $outdir/positives/*.temp       
       echo `cat $outdir/positives/"$contig_hit"_class_"$k"_"$cds_start"_"$cds_end"_positive_hits.dat` `cat $outdir/positives/"$contig_hit"_class_"$k"_"$cds_start"_"$cds_end"_genbank_wseqs.dat` | awk -v com=$comment -v d=$outdir -v ap=$aligned_profile -v aq=$aligned_query -v ape=$aligned_perc -v ins=$insertions -v del=$deletions -v sto=$stops -v from=$hit_hmmfrom -v to=$hit_hmmto '{print $0,com,d,ap,aq,from,to,ape,ins,del,sto}' >> $outdir/positives/"$strain"_positives_total.dat
       rm -rf $outdir/cds_extraction
       # end of sequence extraction      
     endif
    end  
   end
  sed -i '$ d' $outdir/sequence_temp_fas  
  endif
 end

@ i++
end

if ( -f $outdir/positives/"$strain"_positives_total.dat ) then
 echo "POSITIVE HITS FOUND"
else
 echo $strain $genomesize "no_positives" > $outdir/positives/"$strain"_positives_total.dat
endif

if ( -f $outdir/greys/"$strain"_greys_total.dat ) then
 echo "GREY HITS FOUND"
else
 echo $strain $genomesize "no_greys" > $outdir/greys/"$strain"_greys_total.dat
endif


foreach q ($classes)
 set class=class_$q
 set num_greys_$q=`cat $outdir/greys/"$strain"_greys_total.dat | awk -v q=$class '$4 ~ q {print $4}' | wc -l`
 set num_pos_$q=`cat $outdir/positives/"$strain"_positives_total.dat | awk -v q=$class '$4 ~ q {print $4}' | wc -l`
end

echo $strain $num_pos_A $num_pos_C $num_pos_D $num_pos_B1 $num_pos_B2 $num_pos_B3 > $outdir/positives/"$strain"_pos_nums.dat

echo $strain $num_greys_A $num_greys_C $num_greys_D $num_greys_B1 $num_greys_B2 $num_greys_B3 > $outdir/greys/"$strain"_greys_nums.dat

rm $outdir/sequence_temp_fas
rm -rf $outdir/translations
if ( -f null.gbff ) then
rm null.gbff
endif

echo "-------------------------------HMM SEARCH FINISHED---------------------------------------------"
endif


