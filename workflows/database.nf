include { CREATE_DATABASE_FILES } from '../modules/database/create_database_files.nf'

workflow DATABASE {
    take:
    results

    main:
    // Get NTC results for each worklist in their own channel so we can pass it to database file generation later
    split_results_channel = results.branch { sample_id, _worksheet, _referral, _results ->
        ntcs: sample_id.contains("NTC")
        samples: true
    }

    database_creation_ch = split_results_channel.samples.combine(split_results_channel.ntcs, by: 1)

    CREATE_DATABASE_FILES(database_creation_ch) 
}
