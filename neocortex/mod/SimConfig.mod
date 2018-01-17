COMMENT
/**
 * @file SimConfig.mod
 * @brief Interface to write simulation configuration for CoreNEURON
 */
ENDCOMMENT

NEURON {
    THREADSAFE
    ARTIFICIAL_CELL SimConfig
}

VERBATIM
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>

extern double* vector_vec();
extern int vector_capacity();
extern void* vector_arg();
extern int nrnmpi_myid;

// name of config files
const char* SIM_CONFIG_FILE = "sim.conf";
const char* REPORT_CONFIG_FILE = "report.conf";

// helper function to open file and error checking
FILE* open_file(const char *filename, const char *mode) {
    FILE *fp = fopen(filename, mode);
    if(!fp) {
        printf("Error while writing simulation configuration in %s\n", REPORT_CONFIG_FILE);
        abort();
    }
    return fp;
}
ENDVERBATIM

: write report defined in BlueConfig
PROCEDURE write_report_config() {
    VERBATIM
        // gids to be reported is double vector
        double *gid_vec = vector_vec(vector_arg(10));
        int num_gids = vector_capacity(vector_arg(10));

        // copy doible gids to int array
        int *gids = (int*) calloc(num_gids, sizeof(int));
        int i;
        for(i = 0; i < num_gids; i++) {
            gids[i] = (int)gid_vec[i];
        }

        // total number of gids across all ranks in this report
        int tot_num_gids = 0;
        MPI_Allreduce(&num_gids, &tot_num_gids, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

        // write report information
        if(nrnmpi_myid == 0) {
            FILE *fp = open_file(REPORT_CONFIG_FILE, "a");
            fprintf(fp, "\n%s %s %s %s %s %s %lf %lf %lf %d\n",
                        hoc_gargstr(1),
                        hoc_gargstr(2),
                        hoc_gargstr(3),
                        hoc_gargstr(4),
                        hoc_gargstr(5),
                        hoc_gargstr(6),
                        *getarg(7),
                        *getarg(8),
                        *getarg(9),
                        tot_num_gids);
            fclose(fp);
        }

        // write gids using parallel mpi i/o
        MPI_File fh;
        MPI_File_open(MPI_COMM_WORLD, REPORT_CONFIG_FILE, MPI_MODE_APPEND|MPI_MODE_WRONLY, MPI_INFO_NULL, &fh);

        // offset for the current end of file
        MPI_Offset mpi_offset = 0;
        MPI_File_get_position(fh, &mpi_offset);

        // calculate offset considering all ranks
        long num_bytes = num_gids * sizeof(int);
        long offset = 0;
        MPI_Exscan(&num_bytes, &offset, 1, MPI_OFFSET, MPI_SUM, MPI_COMM_WORLD);
        mpi_offset += offset;

        MPI_Status status;
        MPI_File_write_at_all(fh, mpi_offset, gids, num_gids, MPI_INT,  &status);
        MPI_File_close(&fh);
    ENDVERBATIM
}

: Write basic sim settings from Run block of BlueConfig
PROCEDURE write_sim_config() {
VERBATIM
    // should be done by rank 0 only
    if(nrnmpi_myid == 0) {
        FILE *fp = open_file(SIM_CONFIG_FILE, "w");
        fprintf(fp, "--outpath %s\n", hoc_gargstr(1));
        fprintf(fp, "--datpath %s\n", hoc_gargstr(2));
        fprintf(fp, "--tstop %lf\n", *getarg(3));
        fprintf(fp, "--dt %lf\n", *getarg(4));
        fprintf(fp, "--forwardskip %lf\n", *getarg(5));
        fprintf(fp, "--prcellgid %d\n", (int)*getarg(6));
        fprintf(fp, "--report %s\n", REPORT_CONFIG_FILE);
        fprintf(fp, "-mpi\n");
        fclose(fp);
    }
ENDVERBATIM
}

: Write report count as first line
PROCEDURE write_report_count() {
VERBATIM
    // should be done by rank 0 only
    if(nrnmpi_myid == 0) {
        FILE *fp = open_file(REPORT_CONFIG_FILE, "w");
        fprintf(fp, "%d", (int)*getarg(1));
        fclose(fp);
    }
ENDVERBATIM
}
