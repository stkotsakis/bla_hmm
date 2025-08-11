# Nucleotide translator. 1=Input file (plain sequence or fasta format), 2= Frame +
# Stathis D. Kotsakis. Laboratory of Bacteriology. Hellenic Pasteur Institute. 2019
set inp=$1
set fr=$2
set RESBAR=$HOME/resbar
set len=`cat $inp |  sed '/>/d' | tr -d '\n' | wc -c`
set len=`expr $len - $fr + 1`

set end=`echo $len | awk '{ v=$1/3; print int(v)*3}'`

set end=`expr $end + $fr - 1`

set triplets=`cat $inp | sed '/>/d' | tr -d '\n' | cut -c"$fr"-"$end" | sed -e 's/.\{3\}/& /g'`
rm -rf $RESBAR/nuc_translator/temp
mkdir $RESBAR/nuc_translator/temp

# Universal Genetic Code
set Ter= ( TAA TAG TGA )
set Met = (ATG)
set Leu = (TTA TTG CTT CTC CTA CTG CTN)
set Phe = (TTT TTC)
set Ser = (TCT TCC TCA TCG TCN AGT AGC)
set Tyr = (TAT TAC)
set Cys	= (TGT TGC)
set Trp = (TGG)
set Pro = (CCT CCC CCA CCG CCN)
set His = (CAT CAC)
set Gln = (CAA CAG)
set Arg = (CGT CGC CGA CGG CGN AGA AGG)
set Ile = (ATT ATC ATA)
set Thr = (ACT ACC ACA ACG ACN)
set Asn = (AAT AAC)
set Lys = (AAA AAG)
set Val = (GTT GTC GTA GTG GTN)
set Ala = (GCT GCC GCA GCG GCN)
set Asp = (GAT GAC)
set Glu = (GAA GAG)
set Gly = (GGT GGC GGA GGG GGN)

foreach i ($triplets)
 foreach t ($Ter:q)
  if ($i == $t) then
  set triplets=`echo $triplets | sed 's/'$i'/Ter/g'`
  endif
 end
 foreach m ($Met:q)
  if ($i == $m) then
  set triplets=`echo $triplets | sed 's/'$i'/Met/g'`
  endif
  end
 foreach l ($Leu:q)
  if ($i == $l) then
  set triplets=`echo $triplets | sed 's/'$i'/Leu/g'`
  endif
 end
  foreach fe ($Phe:q)
  if ($i == $fe) then
  set triplets=`echo $triplets | sed 's/'$i'/Phe/g'`
  endif
 end
 foreach s ($Ser:q)
  if ($i == $s) then
  set triplets=`echo $triplets | sed 's/'$i'/Ser/g'`
  endif
 end
  foreach y ($Tyr:q)
  if ($i == $y) then
  set triplets=`echo $triplets | sed 's/'$i'/Tyr/g'`
  endif
 end
  foreach c ($Cys:q)
  if ($i == $c) then
  set triplets=`echo $triplets | sed 's/'$i'/Cys/g'`
  endif
 end
  foreach tp ($Trp:q)
  if ($i == $tp) then
  set triplets=`echo $triplets | sed 's/'$i'/Trp/g'`
  endif
 end
  foreach p ($Pro:q)
  if ($i == $p) then
  set triplets=`echo $triplets | sed 's/'$i'/Pro/g'`
  endif
 end
  foreach h ($His:q)
  if ($i == $h) then
  set triplets=`echo $triplets | sed 's/'$i'/His/g'`
  endif
 end
  foreach q ($Gln:q)
  if ($i == $q) then
  set triplets=`echo $triplets | sed 's/'$i'/Gln/g'`
  endif
 end
  foreach r ($Arg:q)
  if ($i == $r) then
  set triplets=`echo $triplets | sed 's/'$i'/Arg/g'`
  endif
 end
  foreach il ($Ile:q)
  if ($i == $il) then
  set triplets=`echo $triplets | sed 's/'$i'/Ile/g'`
  endif
 end
  foreach th ($Thr:q)
  if ($i == $th) then
  set triplets=`echo $triplets | sed 's/'$i'/Thr/g'`
  endif
 end
  foreach n ($Asn:q)
  if ($i == $n) then
  set triplets=`echo $triplets | sed 's/'$i'/Asn/g'`
  endif
 end
  foreach k ($Lys:q)
  if ($i == $k) then
  set triplets=`echo $triplets | sed 's/'$i'/Lys/g'`
  endif
 end
  foreach v ($Val:q)
  if ($i == $v) then
  set triplets=`echo $triplets | sed 's/'$i'/Val/g'`
  endif
 end
  foreach a ($Ala:q)
  if ($i == $a) then
  set triplets=`echo $triplets | sed 's/'$i'/Ala/g'`
  endif
 end
  foreach d ($Asp:q)
  if ($i == $d) then
  set triplets=`echo $triplets | sed 's/'$i'/Asp/g'`
  endif
 end
  foreach e ($Glu:q)
  if ($i == $e) then
  set triplets=`echo $triplets | sed 's/'$i'/Glu/g'`
  endif
 end
  foreach g ($Gly:q)
  if ($i == $g) then
  set triplets=`echo $triplets | sed 's/'$i'/Gly/g'`
  endif
 end
end
echo $triplets > $RESBAR/nuc_translator/temp/tempseq_3letters.temp
cat $RESBAR/nuc_translator/temp/tempseq_3letters.temp | sed -e 's/Met/M/g' | sed -e 's/Ter/\*/g' | sed -e 's/Leu/L/g' | sed -e 's/Phe/F/g' | sed -e 's/Ser/S/g' | sed -e 's/Tyr/Y/g' | sed -e 's/Cys/C/g' | sed -e 's/Trp/W/g' | sed -e 's/Pro/P/g' | sed -e 's/His/H/g' | sed -e 's/Gln/Q/g' | sed -e 's/Arg/R/g' | sed -e 's/Ile/I/g' | sed -e 's/Thr/T/g' | sed -e 's/Asn/N/g' | sed -e 's/Lys/K/g' | sed -e 's/Val/V/g' | sed -e 's/Ala/A/g' | sed -e 's/Asp/D/g' | sed -e 's/Glu/E/g' | sed -e 's/Gly/G/g' | sed -e 's/ //g' | sed -e 's/.\{100\}/&\n/g' > "$inp"_translation





