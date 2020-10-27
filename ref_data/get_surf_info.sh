#!/bin/bash
Z=$(grep -n -m 1 "SAS" msmslog|sed  's/\([0-9]*\).*/\1/')
W=$(grep -n -m 1 "SES_area" msmslog|sed  's/\([0-9]*\).*/\1/')
l_z=$(head -n "$(( Z+1 ))" msmslog| tail -n 1)
l_w=$(head -n "$(( W+1 ))" msmslog| tail -n 1)
l_z_n=$(head -n "$(( Z+1 ))" msmslog| tail -n 1|wc -w)
l_w_n=$(head -n "$(( W+1 ))" msmslog| tail -n 1|wc -w)
echo $l_z | cut -d ' ' -f "$(( l_z_n-1 ))"- > surf_data
echo $l_w | cut -d ' ' -f "$(( l_w_n-1 ))"- >> surf_data
