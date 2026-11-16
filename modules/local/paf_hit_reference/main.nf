process PAF_HIT_REFERENCE {
    tag "${meta.id}_${meta.db_type}_${meta.hit_id}"
    label 'process_single'

    // Same image as nf-core samtools/faidx (samtools 1.25). That image also
    // provides mawk at /usr/bin/awk, which runs paf_hit_reference.awk.
    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/87/87660348a4d032cc08d6ac8101ff9d4d27164300d39c9472043bb5001f9ba801/data'
        : 'community.wave.seqera.io/library/htslib_samtools:1.25--704fb80d75745bd4'}"

    input:
    tuple val(meta), path(paf), path(reference)

    output:
    tuple val(meta), path("parsed.tsv")          , emit: parsed
    tuple val(meta), path("insertion_gaps.tsv")  , emit: gaps
    tuple val(meta), path("best_chromosome.txt") , emit: best_chromosome
    tuple val(meta), path("best_chromosome.fa")  , emit: best_fasta
    tuple val(meta), path("moddotplot.pairs")    , emit: pairs
    tuple val("${task.process}"), val('samtools'), eval("samtools version | sed '1!d;s/.* //'"), topic: versions, emit: versions_samtools

    when:
    task.ext.when == null || task.ext.when

    script:
    def min_insertion = 50
    def max_target_gap = 50
    """
    # The script's shebang is `env awk -f`, which this image's env cannot parse.
    awk -f "\$(command -v paf_hit_reference.awk)" \\
        -v min_insertion=${min_insertion} \\
        -v max_target_gap=${max_target_gap} \\
        ${paf} > parsed.tsv

    printf 'query\\ttarget\\tstrand\\tquery_start\\tquery_end\\ttarget_start\\ttarget_end\\tquery_gap\\ttarget_gap\\n' > insertion_gaps.tsv
    awk -F '\\t' '\$1 == "gap" { print \$2,\$3,\$4,\$5,\$6,\$7,\$8,\$9,\$10 }' OFS='\\t' parsed.tsv >> insertion_gaps.tsv

    # Longest merged syntenic interval from the awk "interval" rows.
    # Summing PAF block length (column 11) prefers a fragmented off-target
    # whose pieces add up to more bases than the chromosome the alignment
    # actually sits on. Ties break on matched bases, then the reference name.
    chrom=\$(awk -F '\\t' -v paf="${paf}" '
        \$1 == "interval" {
            len = \$5 - \$4 + 1
            if (len > span[\$3]) span[\$3] = len
        }
        END {
            while ((getline line < paf) > 0) {
                n = split(line, f, "\\t")
                if (n < 12 || f[1] ~ /^#/) continue
                primary = 1
                for (i = 13; i <= n; i++) if (f[i] ~ /^tp:A:/ && f[i] != "tp:A:P") primary = 0
                if (!primary) continue
                matches[f[6]] += f[10]
            }
            best = ""
            for (target in span) {
                if (best == "" \\
                    || span[target] > span[best] \\
                    || (span[target] == span[best] && matches[target] > matches[best]) \\
                    || (span[target] == span[best] && matches[target] == matches[best] && target < best))
                    best = target
            }
            if (best == "") exit 1
            print best
        }
    ' parsed.tsv)
    printf '%s\\n' "\$chrom" > best_chromosome.txt
    samtools faidx ${reference} "\$chrom" > best_chromosome.fa

    awk -F '\\t' '\$1 == "interval" { print \$3 "\\t" \$2 }' parsed.tsv | sort -u > moddotplot.pairs
    if [[ ! -s moddotplot.pairs ]]; then
        echo "No reference interval for modDotPlot" >&2
        exit 1
    fi
    """

    stub:
    """
    printf 'alignment\\thit\\tchrom_long\\t+\\t1\\t10\\t1\\t10\\n' > parsed.tsv
    printf 'query\\ttarget\\tstrand\\tquery_start\\tquery_end\\ttarget_start\\ttarget_end\\tquery_gap\\ttarget_gap\\n' > insertion_gaps.tsv
    printf 'chrom_long\\n' > best_chromosome.txt
    printf '>chrom_long\\nACGT\\n' > best_chromosome.fa
    printf 'chrom_long\\thit\\n' > moddotplot.pairs
    """
}
