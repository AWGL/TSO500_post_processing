include { FUSION_CHECK } from '../modules/database/fusion_check.nf'
include { WRITE_SAMPLE_DB_LINE ; MERGE_SAMPLES_DATABASE } from '../modules/database/create_samples_database.nf'
include { GET_NTC_READS } from '../modules/database/get_ntc_reads.nf'
include { CONTAMINATION } from '../modules/qc/contamination.nf'

workflow DATABASE {
    take:
    results
    ntc_bams
    sample_qc
    run_id
    rna_samples_list

    main:
    // Get NTC results for each worklist in their own channel so we can pass it to database file generation later
    ntc_channel = results.filter { sample_id, _worksheet, _referral, _results ->
        {
            sample_id.contains("NTC")
        }
    }

    database_creation_ch = results.combine(ntc_channel, by: 1)

    FUSION_CHECK(database_creation_ch)

    GET_NTC_READS(ntc_bams)

    // merge sample_qc, results and ntc bam together
    sample_results_ntc_bam_ch = sample_qc.filter { sample_id, _ws, _referral, _qc -> { !sample_id.contains('NTC') } }.join(GET_NTC_READS.out.ntc_reads, by: 1)

    WRITE_SAMPLE_DB_LINE(run_id, sample_results_ntc_bam_ch)

    MERGE_SAMPLES_DATABASE(
        WRITE_SAMPLE_DB_LINE.out.sample_db_entry.groupTuple().map { worksheet, lines ->
            {
                return tuple(worksheet, lines.join("\n"))
            }
        }
    )

    fusion_check_worksheet_ch = FUSION_CHECK.out.fusion_check.groupTuple(by: 1)

    CONTAMINATION(
        rna_samples_list,
        fusion_check_worksheet_ch,
        file("${projectDir}/RNA_referrals"),
    )
}
