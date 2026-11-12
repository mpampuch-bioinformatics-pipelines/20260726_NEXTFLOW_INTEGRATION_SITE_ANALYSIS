#!/usr/bin/env awk -f
# Syntenic intervals and insertion gaps from a minimap2 PAF.
#
# Produce the PAF with: minimap2 -x asm5 -c --secondary=no reference.fa hit.fa
# Query is the extracted hit. Target is the reference genome. Coordinates
# written here are 1-based and inclusive. PAF coordinates are 0-based
# half-open; target start is converted by adding 1 and target end is kept.
#
# minimap2 -x asm5 usually keeps both host flanks and the construct in one
# record. The construct is then an I operation in the cg tag (700M400I700M),
# not a gap between two records. A query gap between two consecutive
# same-target alignments is the other shape, used when the chain Z-drops.
# Both are reported. An insertion is a query gap of at least min_insertion
# bases whose reference gap is no larger than max_target_gap, on either strand.
#
# Lines:
#   alignment  query target strand q_start q_end t_start t_end
#   interval   query target ref_start ref_end
#   gap        query target strand q_start q_end t_start t_end query_gap target_gap
# interval rows are alignments on one target merged across reference gaps
# of at most max_target_gap. gap t_start is the reference base on the
# low side of the breakpoint (a gene spans the insertion when it covers
# that base and the next one).

BEGIN {
    FS = "\t"
    OFS = "\t"
    if (min_insertion == "") min_insertion = 50
    if (max_target_gap == "") max_target_gap = 50
    min_insertion += 0
    max_target_gap += 0
}

NF < 12 || $1 ~ /^#/ { next }

{
    primary = 1
    cig = ""
    for (i = 13; i <= NF; i++) {
        if ($i ~ /^tp:A:/ && $i != "tp:A:P") primary = 0
        if ($i ~ /^cg:Z:/) cig = $i
    }
    if (!primary) next
    n++
    qname[n] = $1
    qstart[n] = $3 + 0
    qend[n] = $4 + 0
    strand[n] = $5
    tname[n] = $6
    tstart[n] = $8 + 0
    tend[n] = $9 + 0
    cigar[n] = cig
}

function swap(a, b,    t) {
    t = qname[a]; qname[a] = qname[b]; qname[b] = t
    t = qstart[a]; qstart[a] = qstart[b]; qstart[b] = t
    t = qend[a]; qend[a] = qend[b]; qend[b] = t
    t = strand[a]; strand[a] = strand[b]; strand[b] = t
    t = tname[a]; tname[a] = tname[b]; tname[b] = t
    t = tstart[a]; tstart[a] = tstart[b]; tstart[b] = t
    t = tend[a]; tend[a] = tend[b]; tend[b] = t
    t = cigar[a]; cigar[a] = cigar[b]; cigar[b] = t
}

function before(a, b) {
    if (qname[a] < qname[b]) return 1
    if (qname[a] > qname[b]) return 0
    return qstart[a] < qstart[b]
}

# t_low is the 1-based reference base on the low side of the breakpoint.
function emit_gap(query, target, strand, q0, gap_len, t_low) {
    if (t_low < 1) t_low = 1
    print "gap", query, target, strand, q0 + 1, q0 + gap_len, t_low, t_low, gap_len, 0
}

# Walk the CIGAR in query order. For a minus alignment the reference
# cursor starts at the high end and moves down, so the cursor is always
# the low-side base of the next junction. For a plus alignment the
# 0-based cursor is the next reference base, which is the 1-based
# coordinate of the last base already consumed.
function walk_cigar(query, target, strand, q0, t0, tend, cig,    s, tok, len, op, q, t) {
    s = cig
    sub(/^cg:Z:/, "", s)
    q = q0 + 0
    if (strand == "+") t = t0 + 0
    else t = tend + 0
    while (match(s, /[0-9]+[MIDNSHPX=]/)) {
        tok = substr(s, RSTART, RLENGTH)
        len = tok + 0
        op = substr(tok, length(tok), 1)
        if (op == "M" || op == "=" || op == "X") {
            if (strand == "+") t += len
            else t -= len
            q += len
        } else if (op == "I" || op == "S") {
            if (op == "I" && len >= min_insertion) emit_gap(query, target, strand, q, len, t)
            q += len
        } else if (op == "D" || op == "N") {
            if (strand == "+") t += len
            else t -= len
        }
        s = substr(s, RSTART + RLENGTH)
    }
}

END {
    for (a = 2; a <= n; a++) {
        b = a
        while (b > 1 && before(b, b - 1)) {
            swap(b, b - 1)
            b--
        }
    }

    for (i = 1; i <= n; i++) {
        print "alignment", qname[i], tname[i], strand[i], qstart[i] + 1, qend[i], tstart[i] + 1, tend[i]
        if (cigar[i] != "") walk_cigar(qname[i], tname[i], strand[i], qstart[i], tstart[i], tend[i], cigar[i])
        if (i > 1 && qname[i] == qname[i - 1] && tname[i] == tname[i - 1] && strand[i] == strand[i - 1]) {
            qgap = qstart[i] - qend[i - 1]
            if (strand[i] == "+") tgap = tstart[i] - tend[i - 1]
            else tgap = tstart[i - 1] - tend[i]
            if (qgap >= min_insertion && tgap <= max_target_gap && tgap >= -max_target_gap) {
                if (strand[i] == "+") t_low = tend[i - 1]
                else t_low = tstart[i - 1] + 1
                print "gap", qname[i], tname[i], strand[i], qend[i - 1] + 1, qstart[i], t_low, t_low, qgap, tgap
            }
        }
    }

    # Merge reference spans after the per-alignment output. Group by
    # query and target; input is already ordered by query and query start,
    # so re-order each target by reference start while merging.
    for (i = 1; i <= n; i++) {
        key = qname[i] SUBSEP tname[i]
        m = ++nspan[key]
        s_q[key, m] = qname[i]
        s_t[key, m] = tname[i]
        s_a[key, m] = tstart[i] + 1
        s_b[key, m] = tend[i]
    }
    for (key in nspan) {
        m = nspan[key]
        for (a = 2; a <= m; a++) {
            b = a
            while (b > 1 && s_a[key, b] < s_a[key, b - 1]) {
                tmp = s_a[key, b]; s_a[key, b] = s_a[key, b - 1]; s_a[key, b - 1] = tmp
                tmp = s_b[key, b]; s_b[key, b] = s_b[key, b - 1]; s_b[key, b - 1] = tmp
                b--
            }
        }
        cur_a = s_a[key, 1]
        cur_b = s_b[key, 1]
        for (a = 2; a <= m; a++) {
            if (s_a[key, a] <= cur_b + max_target_gap + 1) {
                if (s_b[key, a] > cur_b) cur_b = s_b[key, a]
            } else {
                print "interval", s_q[key, 1], s_t[key, 1], cur_a, cur_b
                cur_a = s_a[key, a]
                cur_b = s_b[key, a]
            }
        }
        print "interval", s_q[key, 1], s_t[key, 1], cur_a, cur_b
    }
}
