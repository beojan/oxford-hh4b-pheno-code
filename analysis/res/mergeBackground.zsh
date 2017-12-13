#!/usr/bin/env zsh
# Merge background samples
# Relies on the SHERPA* directories containing the same histograms

echored () {
    echo -e -n "\e[1;31m"
    echo -n "$@"
    echo -e "\e[0m"
}

echogreen () {
    echo -e -n "\e[1;32m"
    echo -n "$@"
    echo -e "\e[0m"
}

alias yoda-merge=`pwd`/../yoda-merge

numhists=$(ls -1 SHERPA_QCD4b/histo_*.yoda | wc -l)

for dir in SHERPA*; do
    num=$(ls -1 $dir/histo_*.yoda | wc -l)
    if (($num != $numhists)); then
        echored "Directories don't contain the same number of histograms"
	echored "SHERPA_QCD4b has $numhists, $dir has $num"
	if (( $num < $numhists )); then
            exit 1
        else
	    echogreen "$dir has more, continuing anyway"
        fi
    fi
done

rm -rf background
mkdir background

for file in SHERPA_QCD4b/histo_*.yoda(:t:r); do
    echogreen "Merging $file.yoda"
    yoda-merge -o background/$file.yoda -f background/$file.dat SHERPA*/$file.yoda || echored "FAILED TO PROCESS $file.yoda"
done

