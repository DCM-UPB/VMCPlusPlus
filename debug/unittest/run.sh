#!/bin/bash

for folder in ut*
do
   cp run_single_unittest.sh ${folder}/
   cd ${folder}
      echo "-> running unittest ${folder}"
      echo ""
      ./run_single_unittest.sh
   cd ..
done