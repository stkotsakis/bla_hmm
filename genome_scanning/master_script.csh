# SCRIPT FOR DOWNLOADING AND BLA_HMM ANALYSIS OF PROKARYOTIC GENOMES
set domain=$1
set assembly_summary=$2
set cpus=$3

set curr_date=`date +%D | sed 's/\//_/g'`

mkdir $domain
cp $assembly_summary $domain/
mkdir $domain/"$curr_date"_assemblies
mkdir $domain/"$curr_date"_bla_hmm_results

#1 get assembly codes

cat $domain/$assembly_summary | awk -F'\t' '{print $1}' | sort -u | sed '/#/d' > $domain/"$curr_date"_assemblies/assembly_codes.dat

set assemblies_number=`cat $domain/$assembly_summary | awk -F'\t' '{print $1}' | sort -u | sed '/#/d' | wc -l`

if ($assemblies_number >= $cpus) then # split workload onto given number of cpus
 cd $domain/"$curr_date"_assemblies
 set per_cpu=`echo $assemblies_number | awk -v cpus=$cpus '{print int($1/cpus)}'`
 awk -v a=$per_cpu 'NR%a==1{x="ASS"++i;}{print > x}'  assembly_codes.dat
 set file_no=`ls ASS* | wc -l`
 # call parallel script
 cd ../..
 csh parallel_run.csh 1 $file_no $domain $assembly_summary $curr_date
else
cp $domain/"$curr_date"_assemblies/assembly_codes.dat "$curr_date"_assemblies/ASS1
# run genome analysis
mkdir "$curr_date"_logs

csh download_and_bla_hmm.csh 1 $domain $assembly_summary > "$curr_date"_logs/out.1

endif

cat $domain/"$curr_date"_bla_hmm_results/*_positives_total.tab > "$domain"_positive_total.dat
cat $domain/"$curr_date"_bla_hmm_results/*_greys_total.tab > "$domain"_greys_total.dat




 

