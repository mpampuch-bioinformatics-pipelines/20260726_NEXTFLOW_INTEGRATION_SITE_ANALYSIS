#!/usr/bin/env bash -C -e -u -o pipefail
touch cmc-087-1-26_C.merolae_unpurged_raw_primary-busco.batch_summary.txt
mkdir -p cmc-087-1-26_C.merolae_unpurged_raw_primary-busco/cmc-087-1-26_C.merolae_unpurged.p_ctg.fa/run_eukaryota_odb10/busco_sequences

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
nxf_eval_cmd STDOUT STDERR bash -c "busco --version 2> /dev/null | sed 's/BUSCO //g'"
status=$?
if [ $status -eq 0 ]; then
  echo nxf_out_eval_11="$STDOUT" >> .command.env
  echo /nxf_out_eval_11/=exit:0 >> .command.env
else
  echo nxf_out_eval_11="$STDERR" >> .command.env
  echo /nxf_out_eval_11/=exit:$status >> .command.env
fi
