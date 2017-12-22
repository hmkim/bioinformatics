 ThreePrimeUTRs() {
    xtract -pattern INSDSeq -ACC INSDSeq_accession-version -SEQ INSDSeq_sequence \
      -group INSDFeature -if INSDFeature_key -equals CDS -PRD "(-)" \
        -block INSDQualifier -if INSDQualifier_name \
          -equals product -PRD INSDQualifier_value \
        -block INSDFeature -pfc "\n" -element "&ACC" -rst \
          -last INSDInterval_to -element "&SEQ" "&PRD" |
    while read acc pos seq prd
    do
      if [ $pos -lt ${#seq} ]
      then
        echo -e ">$acc 3'UTR: $((pos+1))..${#seq} $prd"
        echo "${seq:$pos}" | fold -w 50
      elif [ $pos -ge ${#seq} ]
      then
        echo -e ">$acc NO 3'UTR"
      fi
    done
  }


cat z2 | ThreePrimeUTRs
