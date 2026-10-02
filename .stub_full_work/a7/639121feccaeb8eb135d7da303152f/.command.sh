#!/usr/bin/env bash -C -e -u -o pipefail
mkdir -p auto_loci_prediction
printf ">stub_locus\nACGTACGTACGT\n" > auto_loci_prediction/1_stub_1-12.fa
printf "locus_id\tsubject\taln_start\taln_end\textract_start\textract_end\tn_members\tmembers\n1\tstub\t1\t12\t1\t12\t1\tstub:1-12\n" > auto_loci_prediction/loci.tsv

# capture process environment
set +u
set +e
cd "$NXF_TASK_WORKDIR"

nxf_eval_cmd() {
    {
        IFS=$'\n' read -r -d '' "${1}";
        IFS=$'\n' read -r -d '' "${2}";
        (IFS=$'\n' read -r -d '' _ERRNO_; return ${_ERRNO_});
    } < <((printf '\0%s\0%d\0' "$(((({ shift 2; "${@}"; echo "${?}" 1>&3-; } | tr -d '\0' 1>&4-) 4>&2- 2>&1- | tr -d '\0' 1>&4-) 3>&1- | exit "$(cat)") 4>&1-)" "${?}" 1>&2) 2>&1)
}

echo '' > .command.env
#
nxf_eval_cmd STDOUT STDERR bash -c "blastdbcmd -version | sed '1!d; s/blastdbcmd: //'"
status=$?
if [ $status -eq 0 ]; then
  echo nxf_out_eval_19="$STDOUT" >> .command.env
  echo /nxf_out_eval_19/=exit:0 >> .command.env
else
  echo nxf_out_eval_19="$STDERR" >> .command.env
  echo /nxf_out_eval_19/=exit:$status >> .command.env
fi
