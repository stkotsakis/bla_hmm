
export protein_ids=`cat class_B2.fasta | awk '$1 ~ ">" {print $1}' | sed 's/>//g' | sed 's/ /\n/g'`
sed -i "\$a>" class_B2.fasta
for i in $protein_ids
 do
 mkdir $i
 cat class_B2.fasta | sed '/^>'$i'/,/^>/{/^>'$i'/!{/^>/!d}}' | sed '/>'$i'/d' > "$i"/class_B2_w_o_"$i".fasta
 cat class_B2.fasta | sed -n '/>'$i'/,/>/p' | sed '$ d' > "$i"/"$i".fasta
 sed -i '$ d' "$i"/class_B2_w_o_"$i".fasta
 clustalo -i "$i"/class_B2_w_o_"$i".fasta -o "$i"/class_B2_w_o_"$i"_aligned.fasta
 hmmbuild "$i"/class_B2_w_o_"$i".hmm "$i"/class_B2_w_o_"$i"_aligned.fasta
 hmmsearch "$i"/class_B2_w_o_"$i".hmm "$i"/"$i".fasta > "$i"/"$i"_vs_model.out
 export score=`cat "$i"/"$i"_vs_model.out | sed -n 15p | awk '{print $2}'`
 echo $i $score >> scores.dat
done
sed -i '$ d' class_B2.fasta
clustalo -i class_B2.fasta -o class_B2_aligned.fasta
hmmbuild class_B2_model.hmm class_B2_aligned.fasta


