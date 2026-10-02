#!/usr/bin/env bash -C -e -u -o pipefail
touch cmc-087-1-26_C.merolae_unpurged_primary.ndb cmc-087-1-26_C.merolae_unpurged_primary.nhr cmc-087-1-26_C.merolae_unpurged_primary.nin cmc-087-1-26_C.merolae_unpurged_primary.not cmc-087-1-26_C.merolae_unpurged_primary.nsq cmc-087-1-26_C.merolae_unpurged_primary.ntf cmc-087-1-26_C.merolae_unpurged_primary.nto

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
nxf_eval_cmd STDOUT STDERR bash -c "blastn -version | sed '1!d; s/blastn: //'"
status=$?
if [ $status -eq 0 ]; then
  echo nxf_out_eval_16="$STDOUT" >> .command.env
  echo /nxf_out_eval_16/=exit:0 >> .command.env
else
  echo nxf_out_eval_16="$STDERR" >> .command.env
  echo /nxf_out_eval_16/=exit:$status >> .command.env
fi
