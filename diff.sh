#!/bin/bash

#parameter setting
samplistA="SRR493366 SRR493367 SRR493368"
samplistB="SRR493369 SRR493370 SRR493371"
numA=3
numB=3
bs=100
fil=5
diff=0.05

echo "Start m5c level differential analysis..."
date

#summary estimated counts
let bsi=bs-1
i=0
for name in ${samplistA}; do
sed '1d' ../results/$name/bootstrap/abundance.tsv >> abundance_A.tsv
while [ $i -le $bsi ]; do
  sed '1d' ../results/$name/bootstrap/bs_abundance_$i.tsv >> bootstrap_A.tsv
  let i=i+1
done
i=0
done
sort -k 1 abundance_A.tsv > abundance_A_order.tsv
rm abundance_A.tsv
sort -k 1 bootstrap_A.tsv > bootstrap_A_order.tsv
rm bootstrap_A.tsv
./sum_counts abundance_A_order.tsv bootstrap_A_order.tsv $numA $bs $fil
mv temp.tsv summary_A.tsv
rm abundance_A_order.tsv
rm bootstrap_A_order.tsv
for name in ${samplistB}; do
sed '1d' ../results/$name/bootstrap/abundance.tsv >> abundance_B.tsv
while [ $i -le $bsi ]; do
  sed '1d' ../results/$name/bootstrap/bs_abundance_$i.tsv >> bootstrap_B.tsv
  let i=i+1
done
i=0
done
sort -k 1 abundance_B.tsv > abundance_B_order.tsv
rm abundance_B.tsv
sort -k 1 bootstrap_B.tsv > bootstrap_B_order.tsv
rm bootstrap_B.tsv
./sum_counts abundance_B_order.tsv bootstrap_B_order.tsv $numB $bs $fil
mv temp.tsv summary_B.tsv
rm abundance_B_order.tsv
rm bootstrap_B_order.tsv

i=0
for name in ${samplistA}; do
sed '1d' ../results_methy/$name/bootstrap/abundance.tsv >> abundance_A.tsv
while [ $i -le $bsi ]; do
  sed '1d' ../results_methy/$name/bootstrap/bs_abundance_$i.tsv >> bootstrap_A.tsv
  let i=i+1
done
i=0
done
sort -k 1 abundance_A.tsv > abundance_A_order.tsv
rm abundance_A.tsv
sort -k 1 bootstrap_A.tsv > bootstrap_A_order.tsv
rm bootstrap_A.tsv
./sum_counts abundance_A_order.tsv bootstrap_A_order.tsv $numA $bs $fil
mv temp.tsv summary_methy_A.tsv
rm abundance_A_order.tsv
rm bootstrap_A_order.tsv
for name in ${samplistB}; do
sed '1d' ../results_methy/$name/bootstrap/abundance.tsv >> abundance_B.tsv
while [ $i -le $bsi ]; do
  sed '1d' ../results_methy/$name/bootstrap/bs_abundance_$i.tsv >> bootstrap_B.tsv
  let i=i+1
done
i=0
done
sort -k 1 abundance_B.tsv > abundance_B_order.tsv
rm abundance_B.tsv
sort -k 1 bootstrap_B.tsv > bootstrap_B_order.tsv
rm bootstrap_B.tsv
./sum_counts abundance_B_order.tsv bootstrap_B_order.tsv $numB $bs $fil
mv temp.tsv summary_methy_B.tsv
rm abundance_B_order.tsv
rm bootstrap_B_order.tsv

echo "complete counting counts."
date

#computing m5C level
./sel_compare summary_A.tsv summary_methy_A.tsv
rm summary_A.tsv
rm summary_methy_A.tsv
mv compare_out summary_A.tsv
mv compare_out_1 summary_methy_A.tsv
./compute_m5c summary_A.tsv summary_methy_A.tsv
mv temp.tsv m5c_A.tsv
./sel_compare summary_B.tsv summary_methy_B.tsv
rm summary_B.tsv
rm summary_methy_B.tsv
mv compare_out summary_B.tsv
mv compare_out_1 summary_methy_B.tsv
./compute_m5c summary_B.tsv summary_methy_B.tsv
mv temp.tsv m5c_B.tsv

echo "complete computing m5c level."
date

#differential analysis for m5C level
./sel_compare m5c_A.tsv m5c_B.tsv
rm m5c_A.tsv
rm m5c_B.tsv
mv compare_out m5c_A.tsv
mv compare_out_1 m5c_B.tsv
./diff_m5c m5c_A.tsv m5c_B.tsv $diff

echo "End m5c level differential analysis."
date
