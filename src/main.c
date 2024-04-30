#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <sys/time.h> 
#include <complex.h>
#include <fftw3.h>

#include "rtl-sdr.h"

void dump_dev_list(uint32_t cnt) {
    for (uint32_t i = 0; i < cnt; i++) {
        const char *name = rtlsdr_get_device_name(i);
        printf("dev #%u, name: %s\n", i, name);
        char manufact[256] = {};
        char product[256] = {};
        char serial[256] = {};
        rtlsdr_get_device_usb_strings(i, manufact, product, serial);
        printf("manufact: '%s', product: '%s', serial: '%s'\n", manufact, product, serial);
    }
}

struct timeval START_TS, STOP_TS;

#define START_MEASURE gettimeofday(&START_TS, NULL);
#define REPORT_MEASURE(_MSG) { \
    gettimeofday(&STOP_TS, NULL); \
    double elapsedTime = (STOP_TS.tv_sec - START_TS.tv_sec) * 1000.0; \
    elapsedTime += (STOP_TS.tv_usec - START_TS.tv_usec) / 1000.0; \
    printf("Duration: %f ms (" # _MSG ").\n", elapsedTime); \
}

// https://superuser.com/questions/80763/is-there-any-usb2-0-data-transfer-chunk-size-limit
// for bulk transfers — 512 bytes for high-speed endpoints, 8, 16, 32, or 64 bytes for full-speed endpoints (and low-speed bulk endpoints are not allowed at all);
#define MAX_USB_BULK_SIZE (512)
#define MAX_READ_SIZE (MAX_USB_BULK_SIZE*2)

/* returns # of IQ pair read or -1 in case of error. */
int read_data_to_void(rtlsdr_dev_t *dev, int npairs) {
    unsigned remains = npairs*2;
    uint8_t buf[MAX_READ_SIZE];
    unsigned total = 0;
    while (remains > 0) {
        unsigned len = MAX_READ_SIZE;
        if (remains < MAX_READ_SIZE) {
            len = remains;
        }
        remains -= len;
        int n_read;
        int ret = rtlsdr_read_sync(dev, buf, len, &n_read);
        if (ret != 0) {
            printf("Cannot read data, error: %d\n", ret);
            return -1;
        } else if (n_read != len) {
            printf("Cannot read data, asked: %d, read: %d\n", len, n_read);
            return -1;
        }
        total += len;
    }
    return total/2;
}

/* returns # of IQ pair read or -1 in case of error. */
int read_data_to_file(rtlsdr_dev_t *dev, FILE *fout, int npairs) {
    unsigned remains = npairs*2;
    uint8_t buf[MAX_READ_SIZE];
    unsigned total = 0;
    while (remains > 0) {
        unsigned len = MAX_READ_SIZE;
        if (remains < MAX_READ_SIZE) {
            len = remains;
        }
        remains -= len;
        int n_read;
        int ret = rtlsdr_read_sync(dev, buf, len, &n_read);
        if (ret != 0) {
            printf("Cannot read data, error: %d\n", ret);
            return -1;
        } else if (n_read != len) {
            printf("Cannot read data, asked: %d, read: %d\n", len, n_read);
            return -1;
        }
        fwrite(buf, 1, len, fout);
        total += len;
    }
    return total/2;
}

/* mem should be at least npairs*2 large.
 * returns # of IQ pair read or -1 in case of error. */
int read_data_to_mem(rtlsdr_dev_t *dev, uint8_t *mem, int npairs) {
    unsigned remains = npairs*2;
    unsigned total = 0;
    while (remains > 0) {
        unsigned len = (remains < MAX_READ_SIZE) ? remains : MAX_READ_SIZE;
        remains -= len;
        int n_read;
        int ret = rtlsdr_read_sync(dev, mem, len, &n_read);
        if (ret != 0) {
            printf("Cannot read data, error: %d\n", ret);
            return -1;
        } else if (n_read != len) {
            printf("Cannot read data, asked: %d, read: %d\n", len, n_read);
            return -1;
        }
        mem += len;
        total += len;
    }
    return total/2;
}

/* converts internal representation to -1 : +1 range */
void normalizef(uint8_t *in, float *out, int niq) {
    for (int i = 0; i < niq; i++) {
        int in_i = (int)in[0] - 127;
        int in_q = (int)in[1] - 127;
        out[0] = in_i / 128.0;
        out[1] = in_q / 128.0;
        in += 2;
        out += 2;
    }
}

void normalize(uint8_t *in, int niq) {
    int8_t *out = (int8_t *)in;
    for (int i = 0; i < niq; i++) {
        int in_i = (int)in[0] - 127;
        int in_q = (int)in[1] - 127;
        out[0] = in_i;
        out[1] = in_q;
        in += 2;
        out += 2;
    }
}

void convert_rtl_to_complex(uint8_t *in, fftw_complex *out, int niq) {
    for (int i = 0; i < niq; i++) {
        __real__ out[0] = in[0] - 127;
        __imag__ out[0] = in[1] - 127;
        in += 2;
        out++;
    }
}

void stat(float *in, int niq, float *imean, float *qmean, float *powerDb) {
    float sumi = 0.0;
    float sumq = 0.0;
    float sump = 0.0;
    for (int i = 0; i < niq; i++) {
        float i = in[0], q = in[1];
        sumi += i;
        sumq += q;
        sump += i*i + q*q;
        in += 2;
    }
    *imean = sumi / (float) niq;
    *qmean = sumq / (float) niq;
    *powerDb = 10.0 * log10(sump / (float) niq);
}

#define SAMP_RATE (2400 * 1000)
#define IQ_100MS (SAMP_RATE / 10)
#define IQ_1S (SAMP_RATE)


#define FREQ (102*1000*1000)

int main(int argv, char **argc) {
    START_MEASURE
    uint32_t cnt = rtlsdr_get_device_count();
    REPORT_MEASURE("rtlsdr_get_device_count")
    printf("Detected devices: %u\n", cnt);
    if (cnt < 0) {
        puts("No device detected");
        return 1;
    }
    dump_dev_list(cnt);

    rtlsdr_dev_t *dev;
    START_MEASURE
    int ret = rtlsdr_open(&dev, 0);
    REPORT_MEASURE("rtlsdr_open")
    if (dev == NULL) {
        printf("Cannot open device #0, code: %d", ret);
        return 1;
    }

    uint32_t rtl_freq;
    uint32_t tuner_freq;
    START_MEASURE
    ret = rtlsdr_get_xtal_freq(dev, &rtl_freq, &tuner_freq);
    REPORT_MEASURE("rtlsdr_get_xtal_freq")
    if (ret != 0) {
        puts("Cannot read XTAL freq");
        goto exit;
    }
    printf("rtl_freq: %u, tuner_freq: %u\n", rtl_freq, tuner_freq);

    uint32_t freq = FREQ;
    START_MEASURE
    ret = rtlsdr_set_center_freq(dev, freq);
    REPORT_MEASURE("rtlsdr_set_center_freq")
    if (ret != 0) {
        printf("Cannot set center freq: %u", freq);
        goto exit;
    }
    START_MEASURE
    uint32_t dev_freq = rtlsdr_get_center_freq(dev);
    REPORT_MEASURE("rtlsdr_get_center_freq")
    if (dev_freq == 0) {
        printf("Cannot read center freq\n");
        goto exit;
    }
    printf("center freq: %u\n", dev_freq);

    START_MEASURE
    int freq_corr_ppm = rtlsdr_get_freq_correction(dev);
    REPORT_MEASURE("rtlsdr_get_freq_correction")
    printf("freq correction (PPM): %d\n", freq_corr_ppm);

    uint32_t samp_rate = SAMP_RATE;
    START_MEASURE
    ret = rtlsdr_set_sample_rate(dev, samp_rate);
    REPORT_MEASURE("rtlsdr_set_sample_rate")
    if (ret != 0) {
        printf("Cannot set sample rate: %u", samp_rate);
        goto exit;
    }
    // samp_rate = rtlsdr_get_sample_rate(dev);
    // if (ret == 0) {
    //     puts("Cannot read sample rate");
    //     goto exit;
    // }
    printf("samp_rate: %u\n", samp_rate);

    int gains[256];
    int gains_cnt;
    START_MEASURE
    gains_cnt = rtlsdr_get_tuner_gains(dev, gains);  // correctly is to call with (dev, NULL) to obtain the number first...
    REPORT_MEASURE("rtlsdr_get_tuner_gains")
    if (gains_cnt <= 0) {
        printf("Cannot get list of gains\n");
        goto exit;
    }
    printf("gains:");
    for (int i = 0; i < gains_cnt; i++) {
        printf(" %i", gains[i]);
    }
    printf("\n");

    START_MEASURE
    ret = rtlsdr_set_tuner_gain_mode(dev, 1);
    REPORT_MEASURE("rtlsdr_set_tuner_gain_mode")
    if (ret != 0) {
        printf("Cannot set manual gain mode\n");
        goto exit;
    }
    START_MEASURE
    ret = rtlsdr_set_agc_mode(dev, 0);
    REPORT_MEASURE("rtlsdr_set_agc_mode")
    if (ret != 0) {
        printf("Cannot disable internal digital AGC\n");
        goto exit;
    }

    START_MEASURE
    rtlsdr_set_direct_sampling(dev, RTLSDR_DS_IQ);
    REPORT_MEASURE("rtlsdr_set_direct_sampling")
    if (ret != 0) {
        printf("Unable to disable direct sampling\n");
        goto exit;
    }

    // Should set after turning off audo AGC, see rtlsdr_set_tuner_gain_mode().
    START_MEASURE
    ret = rtlsdr_set_tuner_gain(dev, gains[2]);
    REPORT_MEASURE("rtlsdr_set_tuner_gain")
    if (ret != 0) {
        printf("Cannot set gain\n");
        goto exit;
    }
    START_MEASURE
    int dev_gain = rtlsdr_get_tuner_gain(dev);
    REPORT_MEASURE("rtlsdr_get_tuner_gain")
    if (dev_gain == 0) {
        printf("Cannot get gain\n");
        goto exit;
    }
    printf("tuner gain: %d\n", dev_gain);

    START_MEASURE
    ret = rtlsdr_reset_buffer(dev);
    REPORT_MEASURE("rtlsdr_reset_buffer")
    if (ret != 0) {
        printf("Cannot reset buffer\n");
        goto exit;
    }

    // START_MEASURE
    // ret = rtlsdr_read_sync(dev, BUF, sizeof(BUF), &n_read);
    // REPORT_MEASURE("rtlsdr_read_sync")
    // if (ret != 0) {
    //     printf("Cannot read data\n");
    //     goto exit;
    // }
    // for (int i = 0; i < n_read; i++) {
    //     printf("%02x", BUF[i]);
    // }
    // size_t FILE_SIZE = 10*1024*1024;
    // FILE *fout = fopen("out.cu8", "wb");
    // START_MEASURE
    // for (int i = 0; i < FILE_SIZE/sizeof(BUF); i++) {
    //     ret = rtlsdr_read_sync(dev, BUF, sizeof(BUF), &n_read);
    //     if (ret != 0) {
    //         printf("Cannot read data\n");
    //         goto exit;
    //     }
    //     if (n_read != sizeof(BUF)) {
    //         printf("Incomplete read\n");
    //         goto exit;
    //     }
    //     fwrite(BUF, sizeof(uint8_t), sizeof(BUF), fout);
    // }
    // fclose(fout);
    // REPORT_MEASURE("writing 10Mb of IQ")

    // FILE *fout_fm = fopen("out_fm.cu8", "wb");
    // int fm_gain = gains[10];
    // rtlsdr_set_tuner_gain(dev, fm_gain);
    // printf("FM gain: %d\n", fm_gain);
    // read_data_to_void(dev, IQ_100MS);
    // read_data_to_file(dev, fout_fm, IQ_1S*10);
    // fclose(fout_fm);

    // FILE *fout_silence = fopen("out_silence.cu8", "wb");
    // for (int i = 0; i < gains_cnt; i++) {
    //     ret = rtlsdr_set_tuner_gain(dev, gains[i]);
    //     if (ret != 0) {
    //         printf("Cannot set gain: %d\n", gains[i]);
    //         goto exit;
    //     }
    //     read_data_to_void(dev, IQ_100MS);
    //     read_data_to_file(dev, fout_silence, IQ_100MS);
    // }
    // fclose(fout_silence);

    // int niq = 100000*2;
    // uint8_t *ubuf = malloc(niq*2);
    // int8_t *ibuf = (int8_t *)ubuf;
    // for (int i = 0; i < gains_cnt; i++) {
    //     int gain = gains[i];
    //     ret = rtlsdr_set_tuner_gain(dev, gain);
    //     if (ret != 0) {
    //         printf("Cannot set gain: %d\n", gain);
    //         goto exit;
    //     }
    //     read_data_to_void(dev, IQ_100MS);
    //     read_data_to_mem(dev, ubuf, niq);
    //     int mini = 0, maxi = 0, minq = 0, maxq = 0;
    //     normalize(ubuf, niq);
    //     for (int i = 0; i < niq; i++) {
    //         if (mini > ibuf[i*2]) { mini = ibuf[i*2]; }
    //         if (maxi < ibuf[i*2]) { maxi = ibuf[i*2]; }
    //         if (minq > ibuf[i*2+1]) { minq = ibuf[i*2+1]; }
    //         if (maxq < ibuf[i*2+1]) { maxq = ibuf[i*2+1]; }
    //     }
    //     printf("gain: %d, diff(i/q): %d/%d, max(i/q): %d/%d, min(i/q): %d/%d\n", gain, maxi-mini, maxq-minq, maxi, maxq, mini, minq);
    // }

    // uint8_t *out_silence_buf = malloc(2*IQ_100MS);
    // float *data_silence = malloc(2*IQ_100MS*sizeof(float));
    // for (int i = 0; i < gains_cnt - 10; i++) {
    //     int gain = gains[i];
    //     ret = rtlsdr_set_tuner_gain(dev, gain);
    //     if (ret != 0) {
    //         printf("Cannot set gain: %d\n", gain);
    //         goto exit;
    //     }
    //     read_data_to_void(dev, IQ_100MS);
    //     read_data_to_mem(dev, out_silence_buf, IQ_100MS);

    //     normalizef(out_silence_buf, data_silence, IQ_100MS);
    //     float imean, qmean, powerDb;
    //     stat(data_silence, IQ_100MS, &imean, &qmean, &powerDb);
    //     printf("Gain: %d, instant power(dB): %f, imean: %f, qmean: %f\n", gain, powerDb, imean, qmean);
    // }


    ret = rtlsdr_set_tuner_gain(dev, gains[18]);
    START_MEASURE
    ret = rtlsdr_set_center_freq(dev, 100*1000*1000);
    REPORT_MEASURE("rtlsdr_set_center_freq")
    if (ret != 0) {
        printf("Cannot set center freq: %u", freq);
        goto exit;
    }
    read_data_to_void(dev, 12*1000); // time to stabilize oscilator: .005s * 2.4e6sps = 12'000samples
    int N = 1024;
    uint8_t *raw_data = malloc(N*2);
    read_data_to_mem(dev, raw_data, N);
    fftw_complex *in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
    fftw_complex *out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
    fftw_plan plan = fftw_plan_dft_1d(N, in, out, FFTW_FORWARD, FFTW_ESTIMATE);
    double *mag = malloc(sizeof(double)*N);
    
    for (int i = 0; i < 10000; i++) {
        START_MEASURE
        read_data_to_mem(dev, raw_data, N);
        convert_rtl_to_complex(raw_data, in, N);
        fftw_execute(plan);
        REPORT_MEASURE("analysis")
        // for (int j = 0; j < N; j++) {
        //     printf("%d:<%.0f,%.0f>, ", j, __real__ in[j], __imag__ in[j]);
        // }
        // printf("\n");
        for (int j = 0; j < N; j++) {
            double power = __real__ out[j] * __real__ out[j] + __imag__ out[j] * __imag__ out[j];
            double magnitude = sqrt(power)/N;
            mag[j] = magnitude;
            // printf("%d:%.3f, ", j, magnitude);
        }
        // printf("\n");

        // collect level around "good" signal
        // 180kHz in both sides is 90kHz around center. How many bind in 90kHz?
        int good_bins = SAMP_RATE / 90000;
        // skip zero bin with DC + SDR spike.
        double p = 0.0;
        for (int j = 1; j < good_bins; j++) {
            p += mag[j];
            p += mag[N - j];
        }
        double mean_good_mag = p / 2 / good_bins;

        // interference on the left is from N/2 to N-good_bins, assume a gap 200kHz wide free of interfence before the signal spectrum
        int gap_bins = SAMP_RATE / 20000;
        int bad_bins = SAMP_RATE / 200000;
        double max_bad_mag = 0.0;
        int trig_idx = 0;
        for (int j = 0; j < bad_bins; j++) {
            int idx = N - good_bins - gap_bins - j;
            double m = mag[idx];
            if (max_bad_mag < m) {
                max_bad_mag = m;
                trig_idx = idx;
            }
        }
        printf("interference: %d, good mag: %f, bad mag: %f\n", mean_good_mag > max_bad_mag, mean_good_mag, max_bad_mag);
        if (mean_good_mag < max_bad_mag) {
            printf("bad idx: %d", trig_idx);
            for (int j = 0; j < N; j++) {
                printf("%.3f, ", mag[i]);
            }
            printf("\n");
            // for (int j = 0; j < N; j++) {
            //     printf("%d:<%.0f,%.0f>, ", j, __real__ in[j], __imag__ in[j]);
            // }
            // printf("\n");
        }
        usleep(500000);
    }

    fftw_destroy_plan(plan);
    fftw_free(in);
    fftw_free(out);
    free(raw_data);

exit:
    rtlsdr_close(dev);
    return 0;
}
