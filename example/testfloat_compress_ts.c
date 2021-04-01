/**
 *  @file test_compress_ts.c
 *  @author Sheng Di
 *  @date May, 2018
 *  @brief This is an example of using compression interface
 *  (C) 2015 by Mathematics and Computer Science (MCS), Argonne National Laboratory.
 *      See COPYRIGHT in top-level directory.
 */


#include <stdio.h>
#include <stdlib.h>
#include "sz.h"
#include "rw.h"
#include <libgen.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>

struct timeval startTime;
struct timeval endTime;  /* Start and end times */
struct timeval costStart; /*only used for recording the cost*/
double totalCost = 0;


void cost_start() {
    totalCost = 0;
    gettimeofday(&costStart, NULL);
}

void cost_end() {
    double elapsed;
    struct timeval costEnd;
    gettimeofday(&costEnd, NULL);
    elapsed = ((costEnd.tv_sec * 1000000 + costEnd.tv_usec) - (costStart.tv_sec * 1000000 + costStart.tv_usec)) /
              1000000.0;
    totalCost += elapsed;
}

void verify(float *ori_data, float *data, size_t num_elements, double *psnr, double *nrmse, double *diffMax) {
    size_t i = 0;
    double Max = ori_data[0];
    double Min = ori_data[0];
    *diffMax = fabs(data[0] - ori_data[0]);
    double diff_sum = 0;
    double maxpw_relerr = fabs((data[0] - ori_data[0]) / ori_data[0]);
    double sum1 = 0, sum2 = 0;
    for (i = 0; i < num_elements; i++) {
        sum1 += ori_data[i];
        sum2 += data[i];
    }
    double mean1 = sum1 / num_elements;
    double mean2 = sum2 / num_elements;

    double sum3 = 0, sum4 = 0;
    double sum = 0, prodSum = 0, relerr = 0;

    double *diff = (double *) malloc(num_elements * sizeof(double));

    for (i = 0; i < num_elements; i++) {
        diff[i] = data[i] - ori_data[i];
        diff_sum += data[i] - ori_data[i];
        if (Max < ori_data[i]) Max = ori_data[i];
        if (Min > ori_data[i]) Min = ori_data[i];
        double err = fabs(data[i] - ori_data[i]);
        if (ori_data[i] != 0) {
            relerr = err / fabs(ori_data[i]);
            if (maxpw_relerr < relerr)
                maxpw_relerr = relerr;
        }

        if (*diffMax < err)
            *diffMax = err;
        prodSum += (ori_data[i] - mean1) * (data[i] - mean2);
        sum3 += (ori_data[i] - mean1) * (ori_data[i] - mean1);
        sum4 += (data[i] - mean2) * (data[i] - mean2);
        sum += err * err;
    }

    double mse = sum / num_elements;
    double range = Max - Min;
    *psnr = 20 * log10(range) - 10 * log10(mse);
    *nrmse = sqrt(mse) / range;

    free(diff);
}

int main(int argc, char *argv[]) {
    int i = 0;
    size_t r1 = 0;
    char *varName;
    if (argc < 6) {
        printf("Test case: testfloat_compress_ts [file] [timesteps] [dimension] [reb] [timestep_interval]\n");
//        printf("Example: testfloat_compress_ts sz.config QCLOUDf /home/sdi/Data/Hurricane-ISA/consecutive-steps 48 500 500 100\n");
        exit(0);
    }

    char cfgFile[600];
    sprintf(cfgFile, "%s/code/sz2/example/sz.config", getenv("HOME"));
    printf("read config file from %s\n", cfgFile);
    int status = SZ_Init(cfgFile);
    if (status == SZ_NSCS)
        exit(0);

    confparams_cpr->errorBoundMode = REL;

    varName = argv[1];
    size_t timesteps = 0;
    int argp = 2;
    timesteps = atoi(argv[argp++]);
    r1 = atoi(argv[argp++]);
    confparams_cpr->relBoundRatio = atof(argv[argp++]);
    confparams_cpr->snapshotCmprStep = atoi(argv[argp++]);


    float *data = (float *) malloc(sizeof(float) * r1);
    SZ_registerVar(1, varName, SZ_FLOAT, data, confparams_cpr->errorBoundMode, confparams_cpr->absErrBound,
                   confparams_cpr->relBoundRatio, confparams_cpr->pw_relBoundRatio, 0, 0, 0, 0, r1);

    size_t outSize, totalOutSize = 0;
    unsigned char *bytes = NULL;
    size_t nbEle;
    float *data_all = readFloatData(varName, &nbEle, &status);
    float *dec_all = (float *) malloc(sizeof(float) * r1 * timesteps);
    double total_compress_time = 0;
    double total_decompress_time = 0;
    for (i = 0; i < timesteps; i++) {
        printf("simulation time step %d\n", i);

        memcpy(data, &data_all[i * r1], r1 * sizeof(float));

        cost_start();
        SZ_compress_ts(SZ_PERIO_TEMPORAL_COMPRESSION, &bytes, &outSize);
        cost_end();
        int currentStep = sz_tsc->currentStep;
        printf("Compress_time = %f\n", totalCost);
        total_compress_time += totalCost;
        totalOutSize += outSize;

        cost_start();
        SZ_decompress_ts(bytes, outSize);
        cost_end();
        sz_tsc->currentStep = currentStep;
        printf("Decompress_time = %f\n", totalCost);
        total_decompress_time += totalCost;
        memcpy(&dec_all[i * r1], data, r1 * sizeof(float));

        printf("Compression Ratio = %.3f\n", r1 * sizeof(float) * 1.0 / outSize);

        free(bytes);
//        free(data_);
    }
    char varNameOutput[600];
    char *folder = dirname(strdup(varName));
    char *filename = basename(strdup(varName));
    sprintf(varNameOutput, "%s/sztime", folder);
    struct stat st = {0};
    if (stat(varNameOutput, &st) == -1) {
        mkdir(varNameOutput, 0700);
    }
    sprintf(varNameOutput, "%s/%s.b%d.%.1e.out", varNameOutput, filename, confparams_cpr->snapshotCmprStep,
            confparams_cpr->relBoundRatio);
    printf("write decompressed file to %s\n", varNameOutput);
    writeByteData((unsigned char *) dec_all, r1 * timesteps * sizeof(float), varNameOutput, &status);

    double psnr, nrmse, max_diff;
    verify(data_all, dec_all, r1, &psnr, &nrmse, &max_diff);
    printf("file=%s, block=%d, compression_ratio=%.3f, reb=%.1e, eb=%.6f, psnr=%.3f, nsmse=%e, compress_time=%.3f, decompress_time=%.3f\n",
           varName, confparams_cpr->snapshotCmprStep,
           timesteps * r1 * sizeof(float) * 1.0 / totalOutSize,
           confparams_cpr->relBoundRatio,
           max_diff, psnr, nrmse,
           total_compress_time, total_decompress_time);

    free(data);
    free(data_all);
    SZ_Finalize();

    return 0;
}
