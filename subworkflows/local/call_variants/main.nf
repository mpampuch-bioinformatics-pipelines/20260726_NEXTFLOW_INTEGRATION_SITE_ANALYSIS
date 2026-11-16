//
// Optional HiFi variant calling.
//
// Small variants: DeepVariant PACBIO on the existing per-hit BAM. The
// reference is that hit FASTA, so coordinates are the extracted window.
// Structural variants: one minimap2 map-hifi alignment per sample to the
// host reference with the plasmid appended as plasmid_ contigs, then Sniffles2.
// An SV call is not an integration call.
//

include { SAMTOOLS_FAIDX                                      } from '../../../modules/nf-core/samtools/faidx/main'
include { DEEPVARIANT_RUNDEEPVARIANT                          } from '../../../modules/nf-core/deepvariant/rundeepvariant/main'
include { BCFTOOLS_STATS as BCFTOOLS_STATS_SMALL              } from '../../../modules/nf-core/bcftools/stats/main'
include { BCFTOOLS_STATS as BCFTOOLS_STATS_SV                 } from '../../../modules/nf-core/bcftools/stats/main'
include { MINIMAP2_ALIGN as MINIMAP2_SV                       } from '../../../modules/nf-core/minimap2/align/main'
include { SNIFFLES                                            } from '../../../modules/nf-core/sniffles/main'
include { PREPARE_SV_REFERENCE                                } from '../../../modules/local/prepare_sv_reference/main'
include { SV_IGV_REPORT                                       } from '../../../modules/local/sv_igv_report/main'

workflow CALL_VARIANTS {
    take:
    ch_alignment // channel: [ val(meta), path(bam), path(bai), path(hit_fasta) ]
    ch_reads     // channel: [ val(meta), path(reads) ] per sample
    ch_host      // channel: [ val(meta), path(host_fasta) ] per sample
    ch_plasmid   // channel: [ val(meta), path(plasmid) ] or an empty list when unset

    main:
    ch_small_vcf = channel.empty()
    ch_small_vcf_index = channel.empty()
    ch_small_gvcf = channel.empty()
    ch_small_gvcf_index = channel.empty()
    ch_small_stats = channel.empty()
    ch_sv_vcf = channel.empty()
    ch_sv_tbi = channel.empty()
    ch_sv_bam = channel.empty()
    ch_sv_bai = channel.empty()
    ch_sv_stats = channel.empty()
    ch_sv_reference = channel.empty()
    ch_sv_igv = channel.empty()

    if (params.call_small_variants) {
        ch_aln = ch_alignment.multiMap { meta, bam, bai, ref ->
            index: [meta, ref, []]
            reads: [meta, bam, bai, []]
            ref: [meta, ref]
        }
        SAMTOOLS_FAIDX(ch_aln.index, false)
        // faidx emits *.fa and *.fasta only. Hit names can be .fna, so keep
        // the original FASTA and attach the .fai it just wrote.
        // DeepVariant names its tuple inputs meta, meta2, meta3, so Nextflow
        // pairs those channels in arrival order. faidx completion reorders
        // the FASTA channel. Join on meta before the split, or one hit's BAM
        // is called against another hit's FASTA.
        ch_indexed = ch_aln.reads
            .join(ch_aln.ref.join(SAMTOOLS_FAIDX.out.fai))
            .multiMap { meta, bam, bai, intervals, fasta, fai ->
                reads: [meta, bam, bai, intervals]
                fasta: [meta, fasta]
                fai: [meta, fai]
                gzi: [meta, []]
                par: [meta, []]
                stats_fasta: [meta, fasta]
            }
        DEEPVARIANT_RUNDEEPVARIANT(
            ch_indexed.reads,
            ch_indexed.fasta,
            ch_indexed.fai,
            ch_indexed.gzi,
            ch_indexed.par,
        )
        // Each process output is one channel. Split it before publishing and
        // before bcftools, or the second consumer never receives the file.
        // bcftools pairs inputs by arrival too, and DeepVariant finishes out
        // of order, so the stats FASTA has to be joined on meta as well.
        ch_small_called = DEEPVARIANT_RUNDEEPVARIANT.out.vcf
            .join(DEEPVARIANT_RUNDEEPVARIANT.out.vcf_index)
            .join(ch_indexed.stats_fasta)
            .multiMap { meta, vcf, index, fasta ->
                vcf: [meta, vcf]
                index: [meta, index]
                stats: [meta, vcf, index]
                regions: [meta, []]
                targets: [meta, []]
                samples: [meta, []]
                exons: [meta, []]
                fasta: [meta, fasta]
            }
        ch_small_vcf = ch_small_called.vcf
        ch_small_vcf_index = ch_small_called.index
        ch_small_gvcf = DEEPVARIANT_RUNDEEPVARIANT.out.gvcf
        ch_small_gvcf_index = DEEPVARIANT_RUNDEEPVARIANT.out.gvcf_index
        BCFTOOLS_STATS_SMALL(
            ch_small_called.stats,
            ch_small_called.regions,
            ch_small_called.targets,
            ch_small_called.samples,
            ch_small_called.exons,
            ch_small_called.fasta,
        )
        ch_small_stats = BCFTOOLS_STATS_SMALL.out.stats
    }

    if (params.call_structural_variants) {
        ch_sv_input = ch_reads
            .map { meta, reads -> [meta.id, meta, reads] }
            .join(ch_host.map { meta, fasta -> [meta.id, fasta] })
            .join(ch_plasmid.map { meta, plasmid -> [meta.id, plasmid] })
            .multiMap { _id, meta, reads, host, plasmid ->
                prepare: [meta, host, plasmid]
                reads: [meta, reads]
            }
        PREPARE_SV_REFERENCE(ch_sv_input.prepare)
        ch_prepared = PREPARE_SV_REFERENCE.out.fasta.multiMap { meta, fasta ->
            reference: [meta, fasta]
            for_reads: [meta.id, fasta]
            for_bam: [meta.id, fasta]
        }
        ch_sv_reference = ch_prepared.reference
        ch_sv_ready = ch_sv_input.reads
            .map { meta, reads -> [meta.id, meta, reads] }
            .join(ch_prepared.for_reads)
            .multiMap { _id, meta, reads, fasta ->
                reads: [meta, reads]
                ref: [meta, fasta]
            }
        MINIMAP2_SV(ch_sv_ready.reads, ch_sv_ready.ref, true, 'bai', false, false)
        ch_mapped = MINIMAP2_SV.out.bam
            .join(MINIMAP2_SV.out.index)
            .multiMap { meta, bam, bai ->
                bam: [meta, bam]
                bai: [meta, bai]
                keyed: [meta.id, meta, bam, bai]
            }
        ch_sv_bam = ch_mapped.bam
        ch_sv_bai = ch_mapped.bai
        ch_sv_bam_ref = ch_mapped.keyed
            .join(ch_prepared.for_bam)
            .multiMap { _id, meta, bam, bai, fasta ->
                sniffles_bam: [meta, bam, bai, [], []]
                sniffles_ref: [meta, fasta]
                tandem: [meta, []]
                stats_ref: [meta, fasta]
                igv: [meta.id, meta, bam, bai, fasta]
            }
        SNIFFLES(
            ch_sv_bam_ref.sniffles_bam,
            ch_sv_bam_ref.sniffles_ref,
            ch_sv_bam_ref.tandem,
            true,
            false,
        )
        // Same arrival-order pairing as DeepVariant. Sniffles can finish in
        // a different order from the reference channel, so join the FASTA
        // on meta before splitting the bcftools inputs.
        ch_sv_called = SNIFFLES.out.vcf
            .join(SNIFFLES.out.tbi)
            .join(ch_sv_bam_ref.stats_ref)
            .multiMap { meta, vcf, tbi, fasta ->
                vcf: [meta, vcf]
                tbi: [meta, tbi]
                stats: [meta, vcf, tbi]
                regions: [meta, []]
                targets: [meta, []]
                samples: [meta, []]
                exons: [meta, []]
                fasta: [meta, fasta]
                igv: [meta.id, vcf]
            }
        ch_sv_vcf = ch_sv_called.vcf
        ch_sv_tbi = ch_sv_called.tbi
        BCFTOOLS_STATS_SV(
            ch_sv_called.stats,
            ch_sv_called.regions,
            ch_sv_called.targets,
            ch_sv_called.samples,
            ch_sv_called.exons,
            ch_sv_called.fasta,
        )
        ch_sv_stats = BCFTOOLS_STATS_SV.out.stats
        SV_IGV_REPORT(
            ch_sv_bam_ref.igv
                .join(ch_sv_called.igv)
                .map { _id, meta, bam, bai, fasta, vcf -> [meta, bam, bai, fasta, vcf] }
        )
        ch_sv_igv = SV_IGV_REPORT.out.report
    }

    emit:
    small_vcf           = ch_small_vcf
    small_vcf_index     = ch_small_vcf_index
    small_gvcf          = ch_small_gvcf
    small_gvcf_index    = ch_small_gvcf_index
    small_stats         = ch_small_stats
    sv_reference        = ch_sv_reference
    sv_bam              = ch_sv_bam
    sv_bai              = ch_sv_bai
    sv_vcf              = ch_sv_vcf
    sv_tbi              = ch_sv_tbi
    sv_stats            = ch_sv_stats
    sv_igv              = ch_sv_igv
}
