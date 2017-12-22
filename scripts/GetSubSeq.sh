#!/bin/bash -norc
while read accn gene start stop
    do
    if [[$start -eq 0 || $stop -eq 0 || $start -eq $stop]]
    then
        echo "Skipping $id due to ambiguous coordinates"
        continue
    fi
    if [ $start -gt $stop ]
    then
        temp=$start
        start=$stop
        start=$temp
        strand=2
    else
        strand=1
    fi
    rslt=`efetch -db nuccore -id $accn -format fasta -seq_start $start -seq_stop $stop -strand $strand < /dev/null`
    echo "$rslt"
done
