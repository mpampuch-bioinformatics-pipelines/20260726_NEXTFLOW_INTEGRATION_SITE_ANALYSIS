#!/usr/bin/env bash
# Map every sequencing read to extracted hit contigs and write a directory
# IGV and JBrowse 2 can open.
#
# Read inspection stays off until a later Nextflow option passes the original
# FASTQ here. This script never downsamples: every read in the FASTQ is aligned.
# Unmapped records are dropped after alignment so the BAM is the reads that
# touch the contigs, which is what both viewers display.
#
# Modules: module load minimap2/2.31 samtools/1.16.1 htslib/1.16
set -euo pipefail

preset=map-hifi
threads="${THREADS:-4}"
outdir=""

usage() {
    echo "usage: read_inspector.sh --outdir DIR [--preset map-hifi] [--threads N] reads.fastq[.gz] hit.fa [hit.fa...]" >&2
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        --outdir)
            outdir="$2"
            shift 2
            ;;
        --preset)
            preset="$2"
            shift 2
            ;;
        --threads)
            threads="$2"
            shift 2
            ;;
        -h|--help)
            usage
            exit 0
            ;;
        --)
            shift
            break
            ;;
        -*)
            echo "Unknown option: $1" >&2
            usage
            exit 1
            ;;
        *)
            break
            ;;
    esac
done

reads="${1:?FASTQ required}"
shift
if [[ $# -lt 1 ]]; then
    echo "At least one hit FASTA is required" >&2
    usage
    exit 1
fi
if [[ -z "$outdir" ]]; then
    echo "--outdir is required" >&2
    usage
    exit 1
fi
if [[ ! -s "$reads" ]]; then
    echo "Reads not found or empty: $reads" >&2
    exit 1
fi

mkdir -p "$outdir"
seen="$(mktemp)"
trap 'rm -f "$seen"' EXIT
: > "$outdir/contigs.fa"

# A colon in a sequence name is an IGV locus, so the coordinate suffix that
# blastdbcmd writes (ptg000011l:441832-517080) is rewritten to an underscore.
# A repeated name is prefixed with the file basename so two hits can share
# a contig id and still be separate sequences.
for hit in "$@"; do
    if [[ ! -s "$hit" ]]; then
        echo "Hit FASTA not found or empty: $hit" >&2
        exit 1
    fi
    base="$(basename "$hit")"
    base="${base%.gz}"
    base="${base%.fa}"
    base="${base%.fasta}"
    base="${base%.fna}"
    awk -v base="$base" -v seen="$seen" '
        BEGIN {
            while ((getline name < seen) > 0) used[name] = 1
            close(seen)
        }
        /^>/ {
            id = $1
            sub(/^>/, "", id)
            gsub(/:/, "_", id)
            orig = id
            if (id in used) id = base "_" id
            n = 2
            while (id in used) {
                id = orig "_" n
                n++
            }
            used[id] = 1
            rest = substr($0, length($1) + 1)
            print ">" id rest
            next
        }
        { print }
        END {
            for (id in used) print id >> seen
        }
    ' "$hit" >> "$outdir/contigs.fa"
done

if [[ ! -s "$outdir/contigs.fa" ]]; then
    echo "No sequences were written from the hit FASTAs" >&2
    exit 1
fi

bgzip -f "$outdir/contigs.fa"
samtools faidx "$outdir/contigs.fa.gz"

# -F 4 drops unmapped reads. They were aligned; they do not belong in the
# viewer track, and keeping them would copy the whole FASTQ into the BAM.
minimap2 -x "$preset" -a -t "$threads" "$outdir/contigs.fa.gz" "$reads" \
    | samtools view -F 4 -u - \
    | samtools sort -@ "$threads" -o "$outdir/reads.bam" -
samtools index "$outdir/reads.bam"

cat > "$outdir/igv_session.xml" << 'EOF'
<?xml version="1.0" encoding="UTF-8" standalone="no"?>
<Session genome="contigs.fa.gz" hasGeneTrack="false" hasSequenceTrack="true" version="8">
    <Resources>
        <Resource path="reads.bam"/>
    </Resources>
</Session>
EOF

cat > "$outdir/config.json" << 'EOF'
{
  "assemblies": [
    {
      "name": "hits",
      "sequence": {
        "type": "ReferenceSequenceTrack",
        "trackId": "hits-ReferenceSequenceTrack",
        "adapter": {
          "type": "BgzipFastaAdapter",
          "fastaLocation": {
            "uri": "contigs.fa.gz",
            "locationType": "UriLocation"
          },
          "faiLocation": {
            "uri": "contigs.fa.gz.fai",
            "locationType": "UriLocation"
          },
          "gziLocation": {
            "uri": "contigs.fa.gz.gzi",
            "locationType": "UriLocation"
          }
        }
      }
    }
  ],
  "tracks": [
    {
      "type": "AlignmentsTrack",
      "trackId": "reads",
      "name": "reads",
      "assemblyNames": ["hits"],
      "adapter": {
        "type": "BamAdapter",
        "bamLocation": {
          "uri": "reads.bam",
          "locationType": "UriLocation"
        },
        "index": {
          "location": {
            "uri": "reads.bam.bai",
            "locationType": "UriLocation"
          },
          "indexType": "BAI"
        }
      }
    }
  ]
}
EOF
