#!/usr/bi/env bash

for na_bbone in 20; do
   for na_sc in 8; do
       if [[ $na_bbone -eq 1 ]] && [[ $na_sc -eq 0 ]] ; then
           continue
        fi
        for fex in '0.001' '0.005' '0.01' '0.025' '0.05' '0.1' '0.2' '0.4' '0.8' \
            '1.6' '3.2' '4.8' '6.4' '8.0' '9.0' '10.0'; do
            #dir_loc=chn-$na_bbone/lk-0/fex-$fex
            dir_loc=bb-$na_bbone-$na_sc/fex-$fex
            mkdir -p $dir_loc
            echo $dir_loc
            for icfg in {1..24}; do
                #python create-chn-ub-teth.py $na_bbone $fex
                python create-chn-br-teth.py $na_bbone $na_sc $fex
                #mv chn-$na_bbone.cfg $dir_loc/chn-$na_bbone.cfg.$icfg
                mv bb-$na_bbone-$na_sc.cfg $dir_loc/bb-$na_bbone-$na_sc.cfg.$icfg
            done
        done
   done
done
