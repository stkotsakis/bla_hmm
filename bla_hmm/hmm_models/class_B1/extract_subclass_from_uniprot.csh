# extract MBL subclasses from uniprot class-b beta-lactamase family based on HMMER results for cross-over validation
set i=$1 # data set
set bla=$2
set cut_off=`sort -k2 -n scores.dat | head -1 | awk '{print int($2)}'`
set hits_acc=`cat "$i"_vs_class_"$bla"_model.tab | awk -v c=$cut_off '$5 > c {print $9}'`
echo ">" >> ../MBL_functionally/"$i"
mkdir uniprot_class_$bla

foreach j ($hits_acc)
 sed -n '/>'$j'/,/>/p' ../MBL_functionally/"$i".fasta > uniprot_class_$bla/"$j"_temp
 sed -i '$ d' uniprot_class_$bla/"$j"_temp
end

cat uniprot_class_$bla/*_temp > uniprot_class_$bla/class_"$bla"_uniprot.fasta
rm uniprot_class_$bla/*_temp

sed -i '$ d' ../MBL_functionally/"$i"

