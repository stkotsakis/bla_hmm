set bla=$1
set PFs=(PF00753 PF12706 PF13483 PF13691 PF14597 PF16661 PF17030)
set cut_off=`sort -k2 -n scores.dat | head -1 | awk '{print int($2)}'`

foreach i ($PFs)
 hmmsearch class_"$bla"_model.hmm ../PFAM_MbL_clan/"$i"_clustered > "$i"_vs_class_"$bla"_model.dat
 cat "$i"_vs_class_"$bla"_model.dat | sed -n '/Scores for complete sequences/,/Domain annotation for each/p' | sed '/------ inclusion threshold ------/d' | sed '1,4d' | head -n -3 > "$i"_vs_class_"$bla"_model.tab
 set seqs=`cat ../PFAM_MbL_clan/"$i"_clustered | awk '$1 ~ ">" {print $1}' | wc -l`
 set hit_check=`cat "$i"_vs_class_"$bla"_model.tab | awk '$0 ~ "No hits detected that satisfy reporting thresholds" {print $0}' | wc -l`
 if ($hit_check == 0) then
 set hits=`cat "$i"_vs_class_"$bla"_model.tab | awk '{print $1}' | wc -l`
 set hits_co=`cat "$i"_vs_class_"$bla"_model.tab | awk -v c=$cut_off '$2 > c {print $2}' | wc -l`
 else
  set hits=0
  set hits_co=0
 endif
 echo "$i"" vs class_"$bla"_model cut-off=""$cut_off"" ""$seqs"" / ""$hits"" / ""$hits_co" >> class_"$bla"_model_vs_PFs.dat
end

