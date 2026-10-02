#!/usr/bin/env bash -C -e -u -o pipefail
mkdir -p all_hits
printf ">stub_hit\nACGTACGTACGT\n" > all_hits/1_stub.fa
printf ">stub_hit\nACGTACGTACGT\n" > cmc-087-1-35_C.merolae_unpurged_combined.extracted.fa

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
  echo nxf_out_eval_18="$STDOUT" >> .command.env
  echo /nxf_out_eval_18/=exit:0 >> .command.env
else
  echo nxf_out_eval_18="$STDERR" >> .command.env
  echo /nxf_out_eval_18/=exit:$status >> .command.env
fi
