include { APP_DEMULTIPLEX } from '../modules/tso500_app/demultiplex.nf'

workflow DEMULTIPLEX {
    take:
    sample_sheet
    resources
    run_folder

    main:
    APP_DEMULTIPLEX(
        sample_sheet,
        resources,
        run_folder
    )
}

