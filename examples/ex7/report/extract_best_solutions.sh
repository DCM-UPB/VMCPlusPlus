#!/bin/bash


numCompare() {
   awk -v n1="$1" -v n2="$2" 'BEGIN {exit (n1 > n2)}' /dev/null
}


for folder in "1l2u" "1l3u" "1l5u" "1l10u" "2l2u" "2l3u" "2l5u" "2l10u"
do
    cd ${folder}

    echo "==================================================="
    echo "-> ${folder}"

    best_energy=1000.
    best_energy_error=1000.
    best_subfolder=1
    for subfolder in {1..8}
    do
        cd ${subfolder}
        echo "    -> ${subfolder}"
        sed '37q;d' output.txt | sed "s/Total/    Total/g"

        energy=$(sed '37q;d' output.txt | sed "s/       Total Energy        = //g" | sed "s/ +- [0-9.]*//g")
        energy_error=$(sed '37q;d' output.txt | sed "s/       Total Energy        = [0-9.]* +- //g")

        if numCompare $energy $best_energy
        then
            # echo "better energy"
            if numCompare $energy_error $best_energy_error
            then
                # echo "better energy error"
                best_energy=$energy
                best_energy_error=$energy_error
                best_subfolder=$subfolder
            # else
            #     echo "worse energy error"
            fi
        # else
        #     echo "worse energy"
        fi

        cd ..
    done

    echo "==================================================="
    echo "best energy = $best_energy +- $best_energy_error"
    echo "==================================================="

    rm -f -r best
    mkdir best
    cp $best_subfolder/*.txt best/

    cd best
        echo "set term png" >> gnuplot.in

        echo "set output \"foo.png\"" >> gnuplot.in
        echo "plot \"plot_opt_wf.txt\"" >> gnuplot.in

        echo "set output \"plot.png\"" >> gnuplot.in
        echo "plot \"plot_opt_wf.txt\" title \"$folder\", GPVAL_DATA_Y_MAX*exp(-0.5*x*x) title \"exact solution\"" >> gnuplot.in

        gnuplot gnuplot.in
        rm -f foo.png
    cd ..

    cd ..
done
