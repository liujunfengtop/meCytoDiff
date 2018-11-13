#!/bin/bash

#create directoty
mkdir ../diff_command

#compile C code
echo "Start compile ..."
date

gcc -o ../diff_command/anti_bisulfite anti_bisulfite.c
gcc -o ../diff_command/anti_bisulfite_single_batch anti_bisulfite_single_batch.c -lm
gcc -o ../diff_command/combination_diff combination_diff.c
gcc -o ../diff_command/compute_m5c compute_m5c.c -lm
gcc -o ../diff_command/diff_m5c diff_m5c.c -lm
gcc -o ../diff_command/sel_compare sel_compare.c
gcc -o ../diff_command/m5c_filter m5c_filter.c
gcc -o ../diff_command/selmethy selmethy.c -lm
gcc -o ../diff_command/sum_counts sum_counts.c -lm
cp anti_bisulfite.ctl ../diff_command
cp anti_bisulfite_single_batch.ctl ../diff_command
cp diff.sh ../diff_command/
cp diff_whole.sh ../diff_command/
cp diff_single.sh ../diff_command/
cp kallisto ../diff_command/
chmod u=rwx,g=rx,o=x ../diff_command/diff.sh
chmod u=rwx,g=rx,o=x ../diff_command/diff_whole.sh
chmod u=rwx,g=rx,o=x ../diff_command/diff_single.sh
chmod u=rwx,g=rx,o=x ../diff_command/kallisto

echo "End compile."
date
