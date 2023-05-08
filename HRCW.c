/*
   Copyright 2018 Olivier SCHMITT

   This program is free software: you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#include "HRCW.h"

#include <stdio.h>
#include <stdlib.h>
#include <GL/glut.h>
#include <alsa/asoundlib.h>
#include <complex.h>
#include <fftw3.h>
#include <math.h>
#include <string.h>

#include <iostream>     // std::cout
#include <algorithm>    // std::random_shuffle
#include <vector>       // std::vector
#include <ctime>        // std::time
#include <cstdlib>      // std::rand, std::srand

#include <sys/time.h>

#include <pthread.h>


int SOUND_RATE = 48000;
char *SOUND_DEVICE = "default";
int SOUND_SAMPLES_PER_TURN = 2048;
float *HANNINGWINDOWS;
double filterth[500][2]={0};
int *peak = NULL;

int *aleath = {0};
int threshold_fact = 4;

// pthread_t decoders[max_index_thread];
pthread_t *decoders;
pthread_mutex_t mutex_global_thplan;
struct thplanst {
        int index;
        bool buffer_ready;
        bool is_running;
        double currentcolG[2048];
        double currentTH[2048];
        bool currentdTH[2048];
};
thplanst *thplan;
// int thplan[max_index_thread][3] = {false,0,false};

struct thargs {
        int index;
        int fftindexhz;
};


/* Informations about the window, display options. */
struct interactionInfo
{
        int width;
        int height;

        int update;

        bool showOvertones;
        bool doPanning;
        bool forceOverview;
        bool showMainGrid;
        bool showWaveform;
        bool showFrequency;
        bool frequencyLabelLeft;

        int showthread;

        int lastMouseDownBS[2];
        int lastMouseDownES[2];
        double lastMouseDownBW[2];
        double lastMouseDownEW[2];

        double offsetX, lastOffsetX;
        double scaleX;
} interaction;

/* Global sound info. */
struct soundInfo
{
        snd_pcm_t *handle;

        char *buffer, *bufferLast;
        snd_pcm_uframes_t bufferSizeFrames;
        snd_pcm_uframes_t bufferFill;
        int bufferReady;

        int reprepare;
} sound;

/* Global fftw info. */
struct fftwInfo
{
        fftw_complex *in;
        fftw_complex *out;
        fftw_plan plan;
        int outlen;
        double binWidth;

        double *currentLine;

        unsigned char *textureData;
        GLuint textureHandle;
        int textureWidth, textureHeight;
} fftw;





/* Check for OpenGL errors. */
void checkError(int line)
{
        GLenum err = glGetError();
        switch (err)
        {
        case GL_NO_ERROR:
                break;

        case GL_INVALID_ENUM:
                fprintf(stderr, "GL_INVALID_ENUM: %d\n", line);
                break;
        case GL_INVALID_VALUE:
                fprintf(stderr, "GL_INVALID_VALUE: %d\n", line);
                break;
        case GL_INVALID_OPERATION:
                fprintf(stderr, "GL_INVALID_OPERATION: %d\n", line);
                break;
        case GL_STACK_OVERFLOW:
                fprintf(stderr, "GL_STACK_OVERFLOW: %d\n", line);
                break;
        case GL_STACK_UNDERFLOW:
                fprintf(stderr, "GL_STACK_UNDERFLOW: %d\n", line);
                break;
        case GL_OUT_OF_MEMORY:
                fprintf(stderr, "GL_OUT_OF_MEMORY: %d\n", line);
                break;
        case GL_TABLE_TOO_LARGE:
                fprintf(stderr, "GL_TABLE_TOO_LARGE: %d\n", line);
                break;
        }
}

/* Get i'th sample from buffer and convert to short int. */
short int getFrame(char *buffer, int i, short int channel)
{
        return (buffer[(4 * i) + channel] & 0xFF) + ((buffer[(4 * i) + 1 + channel] & 0xFF) << 8);
}

/* Return the environment variable "name" or "def" if it's unset. */
char *getenvDefault(char *name, char *def)
{
        char *val = getenv(name);
        if (val == NULL)
                return def;
        else
                return val;
}

void threadinit(){
        for (int z=0; z<max_index_thread; z++) {

                thplan[z].index=0;
                thplan[z].buffer_ready=false;
                thplan[z].is_running=false;
                memset(thplan[z].currentcolG, 0, sizeof(thplan[z].currentcolG));
                memset(thplan[z].currentTH, 0, sizeof(thplan[z].currentcolG));
                memset(thplan[z].currentdTH, false, sizeof(thplan[z].currentdTH));

                // thplan[z][0]=false; //is_running
                // thplan[z][1]=0; //index
                // thplan[z][2]=false; //buffer_ready

        }
}

/* Open and init the default recording device. */
void audioInit(void)
{
        int rc;
        int size;
        snd_pcm_hw_params_t *params;
        unsigned int val;
        int dir = 0;

        /* Open PCM device for recording (capture). */
        rc = snd_pcm_open(&sound.handle, getenvDefault(SOUND_DEVICE_ENV,
                                                       SOUND_DEVICE),
                          SND_PCM_STREAM_CAPTURE, 0);
        if (rc < 0)
        {
                fprintf(stderr, "unable to open pcm device: %s\n", snd_strerror(rc));
                exit(EXIT_FAILURE);
        }

        /* Allocate a hardware parameters object. */
        snd_pcm_hw_params_alloca(&params);

        /* Fill it in with default values. */
        snd_pcm_hw_params_any(sound.handle, params);

        /* Set the desired hardware parameters. */

        /* Interleaved mode. */
        snd_pcm_hw_params_set_access(sound.handle, params,
                                     SND_PCM_ACCESS_RW_INTERLEAVED);

        /* Signed 16-bit little-endian format. */
        snd_pcm_hw_params_set_format(sound.handle, params,
                                     SND_PCM_FORMAT_S16_LE);

        /* One channel (mono). */
        snd_pcm_hw_params_set_channels(sound.handle, params, 2);

        /* SOUND_RATE bits/second sampling rate (CD quality). */
        val = SOUND_RATE;
        snd_pcm_hw_params_set_rate_near(sound.handle, params, &val, &dir);

        /* Set period size. It's best to set this to the same value as
         * SOUND_SAMPLES_PER_TURN. A lower value would generate more
         * hardware interrupts and thus a lower latency but that's of no use
         * since we have to wait for an amount of SOUND_SAMPLES_PER_TURN
         * samples anyway. */
        snd_pcm_uframes_t frames = SOUND_SAMPLES_PER_TURN;
        snd_pcm_hw_params_set_period_size_near(sound.handle, params,
                                               &frames, &dir);

        /* Write the parameters to the driver. */
        rc = snd_pcm_hw_params(sound.handle, params);
        if (rc < 0)
        {
                fprintf(stderr, "unable to set hw parameters: %s\n",
                        snd_strerror(rc));
                exit(EXIT_FAILURE);
        }

        /* Acquire n frames per turn. */
        sound.bufferSizeFrames = SOUND_SAMPLES_PER_TURN;
        size = sound.bufferSizeFrames * 2 * 2; /* 2 bytes/sample, 2 channel */

        /* Initialize the buffer. */
        sound.buffer = (char *)malloc(size);
        sound.bufferLast = (char *)malloc(size);
        sound.bufferFill = 0;
        sound.bufferReady = 0;


        /* Try to switch to non-blocking mode for reading. If that fails,
         * print a warning and continue in blocking mode. */
        rc = snd_pcm_nonblock(sound.handle, 1);
        if (rc < 0)
        {
                fprintf(stderr, "Could not switch to non-blocking mode: %s\n",
                        snd_strerror(rc));
        }

        /* Prepare in audioRead() for the first time. */
        sound.reprepare = 1;
}

/* Read as far as you can (when in non-blocking mode) or until our
 * buffer is filled (when in blocking mode). */
int audioRead(void)
{
        if (sound.reprepare)
        {
                int ret;
                ret = snd_pcm_drop(sound.handle);
                if (ret < 0)
                {
                        fprintf(stderr, "Error while dropping samples: %s\n",
                                snd_strerror(ret));
                }

                ret = snd_pcm_prepare(sound.handle);
                if (ret < 0)
                {
                        fprintf(stderr, "Error while preparing to record: %s\n",
                                snd_strerror(ret));
                }

                sound.reprepare = 0;
        }

        /* Request
         *   "size - fill" frames
         * starting at
         *   "base + numFramesFilled * 2" bytes.
         * Do "* 2" because we get two bytes per sample.
         *
         * When in blocking mode, this always fills the buffer to its
         * maximum capacity.
         */
        snd_pcm_sframes_t rc;
        rc = snd_pcm_readi(sound.handle, sound.buffer + (sound.bufferFill * 2 *2),
                           sound.bufferSizeFrames - sound.bufferFill);
        if (rc == -EPIPE)
        {
                /* EPIPE means overrun */
                snd_pcm_recover(sound.handle, rc, 0);
        }
        else if (rc == -EAGAIN)
        {
                /* Not ready yet. Come back again later. */
        }
        else if (rc < 0)
        {
                fprintf(stderr, "error from read: %s\n", snd_strerror(rc));
        }
        else
        {
                sound.bufferFill += rc;
                if (sound.bufferFill == sound.bufferSizeFrames)
                {
                        /* Buffer full. display() can add this to the history. */
                        sound.bufferFill = 0;
                        sound.bufferReady = 1;
                }
        }

        return rc;
}

/*hanning function*/
/*hanning function*/
float *hanningInit(int N, short itype = 0)
{
        int half, i, idx, n;
        float *w;

        w = (float*) calloc(N, sizeof(float));
        memset(w, 0, N*sizeof(float));

        if(itype==1) //periodic function
                n = N-1;
        else
                n = N;

        if(n%2==0)
        {
                half = n/2;
                for(i=0; i<half; i++) //CALC_HANNING   Calculates Hanning window samples.
                        w[i] = 0.5 * (1 - cos(2*M_PI*(i+1) / (n+1)));

                idx = half-1;
                for(i=half; i<n; i++) {
                        w[i] = w[idx];
                        idx--;
                }
        }
        else
        {
                half = (n+1)/2;
                for(i=0; i<half; i++) //CALC_HANNING   Calculates Hanning window samples.
                        w[i] = 0.5 * (1 - cos(2*M_PI*(i+1) / (n+1)));

                idx = half-2;
                for(i=half; i<n; i++) {
                        w[i] = w[idx];
                        idx--;
                }
        }

        if(itype==1) //periodic function
        {
                for(i=N-1; i>=1; i--)
                        w[i] = w[i-1];
                w[0] = 0.0;
        }
        return(w);
}


int *randomth(int N)
{
        srand(time(NULL));
        bool intused[N]={0};
        int *number;
        number = (int*) calloc(N, sizeof(int));
        memset(number, 0, N*sizeof(int));

        int p=0;
        for (int i = 0; i < N; i++) {
                while(intused[p]) {p = rand()%N;}
                intused[p]=1;
                number[i]=p;
        }
        return number;
}


/* Shutdown audio device. */
void audioDeinit(void)
{
        snd_pcm_drop(sound.handle);
        snd_pcm_close(sound.handle);
        free(sound.buffer);
        free(sound.bufferLast);
}

/* Create FFTW-plan, allocate buffers. */
void fftwInit(void)
{
        fftw.outlen = sound.bufferSizeFrames;

        fftw.in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * sound.bufferSizeFrames);
        fftw.out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (sound.bufferSizeFrames+1));

//	fftw.in = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * sound.bufferSizeFrames);
//	fftw.out = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * (fftw.outlen + 1));

        fftw.plan = fftw_plan_dft_1d(sound.bufferSizeFrames, fftw.in, fftw.out,FFTW_FORWARD,
                                     FFTW_MEASURE);

        fftw_init_threads();
        fftw_plan_with_nthreads(FFTW_NO_NONTHREADED);

        fftw.currentLine = (double *)malloc(sizeof(double) * fftw.outlen);
        memset(fftw.currentLine, 0, sizeof(double) * fftw.outlen);

        fftw.textureWidth = fftw.outlen;
        fftw.textureHeight = FFTW_HISTORY_SIZE;
        fftw.textureData = (unsigned char *)malloc(sizeof(unsigned char)
                                                   * fftw.textureWidth
                                                   * fftw.textureHeight * 3);
        memset(fftw.textureData, 0, sizeof(unsigned char)
               * fftw.textureWidth * fftw.textureHeight * 3);

        /* How many hertz does one "bin" comprise? */
        fftw.binWidth = (double)SOUND_RATE / (double)fftw.outlen;
        printf ("FFT Resolution = %f\n",fftw.binWidth);
}

/* Free any FFTW resources. */
void fftwDeinit(void)
{
        fftw_destroy_plan(fftw.plan);
        fftw_free(fftw.in);
        fftw_free(fftw.out);
        free(fftw.currentLine);
        free(fftw.textureData);
        fftw_cleanup();
}

void reverse(char* str) {
        int len = strlen(str);
        for (int i = 0; i < len / 2; i++) {
                char temp = str[i];
                str[i] = str[len - i - 1];
                str[len - i - 1] = temp;
        }
}

void morse_to_text(char* morse, char* buffer) {
        const char* morse_table[] = {".-", "-...", "-.-.", "-..", ".", "..-.", "--.", "....", "..", ".---", "-.-", ".-..", "--", "-.", "---", ".--.", "--.-", ".-.", "...", "-", "..-", "...-", ".--", "-..-", "-.--", "--..", ".----", "..---", "...--", "....-", ".....", "-....", "--...", "---..", "----.", "-----", "--..--", ".-.-.-", "..--..", "-.-.--", "-....-", ".-.-.", ".----.", "-..-.", "-.--.", "-.--.-", "/"}; // table de correspondance morse/ASCII
        const char* ascii_table[] = {"A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K", "L", "M", "N", "O", "P", "Q", "R", "S", "T", "U", "V", "W", "X", "Y", "Z", "1", "2", "3", "4", "5", "6", "7", "8", "9", "0", ",", ".", "?", "!", "-", "/", "(", ")", "&", ":", ";", "=", "+", "_", "$", "@", " "}; // table de correspondance ASCII/morse

        // initialisation du buffer
        buffer[0] = '\0';

        // traitement de la chaîne de caractères morse
        char* token = strtok(morse, " "); // séparation des caractères en morse
        while (token != NULL) {
                if (strcmp(token, "~") == 0) { // ajout d'un espace si le caractère est ~
                        if (strlen(buffer) + 1 < 2048) { // vérification de la taille du buffer
                                strcat(buffer, " ");
                        }
                } else {
                        int i;
                        for (i = 0; i < 47; i++) {
                                if (strcmp(token, morse_table[i]) == 0) { // comparaison avec la table de correspondance morse/ASCII
                                        if (strlen(buffer) + 1 < 2048) { // vérification de la taille du buffer
                                                strcat(buffer, ascii_table[i]); // ajout du caractère ASCII correspondant au buffer
                                        }
                                        break;
                                }
                        }
                }
                token = strtok(NULL, " "); // passage au caractère suivant
        }
}


void detect_morse_code(bool signals[], char* outbuffer) {

        int distrizero[fftw.outlen];
        int distrione[fftw.outlen];
        // Détection des points et traits
        int j=0;
        int k=0;
        int l=0;
        bool lastsig = 0;
        for (int i = 0; i < fftw.outlen; i++) {
                if (signals[i] == lastsig) {
                        j++;
                }else{
                        lastsig=signals[i];
                        if(lastsig==1) {
                                distrione[j]++;
                        }
                        if(lastsig==0) {
                                distrizero[j]++;
                        }
                        j=0;
                }
        }

        int val = 0;
        int maxval = 0;
        for (int i = 1; i < fftw.outlen; i++) {
                if(distrione[i] > maxval) {maxval=distrione[i]; val=i;}
        }


        int val2 = val*2;
        int maxval2 = 0;
        for (int i = val2; i < fftw.outlen; i++) {
                if(distrione[i] > maxval2) {maxval2=distrione[i]; val2=i;}
        }

        fflush(stdout);

        int val3 = val/2;
        int maxval3 = 0;
        for (int i = val3; i > 0; i--) {
                if(distrione[i] > maxval3) {maxval3=distrione[i]; val3=i;}
        }

        int thresoldcw = 0;

        if(maxval2 > maxval3) {
                thresoldcw =  val+val2/2;
                // printf("dit %i dash %i th %i\n",val,val2,thresoldcw);
        }
        else{
                thresoldcw =  val+val3/2;
                // printf("dit %i dash %i th %i\n",val3,val,thresoldcw);
        }


        char cw[2048] = "";
        cw[0] = '\0';
        int cwindex = 0;
        for (int i = 0; i < fftw.outlen; i++) {
                if ((signals[i] == lastsig) | (j<3)) {
                        // if(lastsig==1) {printf("°");}
                        // if(lastsig==0) {printf("_");}
                        j++;
                }else{
                        if(lastsig==1) {
                                if(j>=(thresoldcw)) {
                                        // printf("-");
                                        cw[cwindex]=(char)45;
                                        cwindex++;
                                }
                                else{
                                        // printf(".");
                                        cw[cwindex]=(char)46;
                                        cwindex++;
                                }

                        }
                        else if(lastsig==0) {
                                if(j>thresoldcw*2) {
                                        // printf(" ");
                                        cw[cwindex]=(char)32;
                                        cwindex++;
                                        cw[cwindex]=(char)126;
                                        cwindex++;
                                        cw[cwindex]=(char)32;
                                        cwindex++;
                                }else if(j>thresoldcw) {
                                        // printf(" ");
                                        cw[cwindex]=(char)32;
                                        cwindex++;
                                }
                        }

                        lastsig=signals[i];
                        j=0;
                }
        }
        cw[cwindex]='\0';
        reverse(cw);
        morse_to_text(cw,outbuffer);
}

double compute_noise_level(double* measurements, int num_measurements) {
        double mean = 0.0, variance = 0.0, stddev = 0.0, noise_level = 0.0;
        int i;

        // Calculer la moyenne
        for(i = 0; i < num_measurements; i++) {
                mean += measurements[i];
        }
        mean /= num_measurements;

        // Calculer l'écart type
        for(i = 0; i < num_measurements; i++) {
                variance += pow(measurements[i] - mean, 2.0);
        }
        variance /= num_measurements;
        stddev = sqrt(variance);

        // Retirer les valeurs qui dépassent 3 fois l'écart type
        int num_valid_measurements = 0;
        double valid_measurements[num_measurements];
        for(i = 0; i < num_measurements; i++) {
                if(fabs(measurements[i] - mean) <= 3.0 * stddev) {
                        valid_measurements[num_valid_measurements++] = measurements[i];
                }
        }

        // Calculer la moyenne des valeurs restantes pour obtenir la valeur de bruit de fond
        for(i = 0; i < num_valid_measurements; i++) {
                noise_level += valid_measurements[i];
        }
        noise_level /= num_valid_measurements;

        return noise_level;
}

void* decoderthfc(void *input)
{

        int maxwaiketime = 60000; //temps de vie du thread

        //on parse les infos
        int index = ((struct thargs*)input)->index;
        int fftindex = ((struct thargs*)input)->fftindexhz;

        printf("new th %i on %f\n", index, (float)(((struct thargs*)input)->fftindexhz * fftw.binWidth));

        double currentcol[fftw.outlen];
        int sizeofcurrentcol = (fftw.outlen-1) * sizeof(double);
        int sizeofcurrentcolbool = (fftw.outlen-1) * sizeof(bool);

        double threasold = 0;

        int zz=0;

        while (maxwaiketime > 0)
        {
                while(!thplan[index].buffer_ready) {
                        usleep(1000);
                        maxwaiketime--;
                }

                memmove(&currentcol[1], &currentcol[0], sizeofcurrentcol);
                memmove(&thplan[index].currentcolG[1], &thplan[index].currentcolG[0], sizeofcurrentcol);
                thplan[index].currentcolG[0]=fftw.currentLine[fftindex];

                memmove(&thplan[index].currentTH[1], &thplan[index].currentTH[0], sizeofcurrentcol);

                currentcol[0]=fftw.currentLine[fftindex];

                threasold =  compute_noise_level(currentcol,fftw.outlen/8);
                thplan[index].currentTH[0]=threasold;

                memmove(&thplan[index].currentdTH[1], &thplan[index].currentdTH[0], sizeofcurrentcolbool);

                if(currentcol[0]>threasold) {thplan[index].currentdTH[0]=1;}
                else{thplan[index].currentdTH[0]=0;}

                char buffer[2048];
                if(zz % 1024) {detect_morse_code(thplan[index].currentdTH,buffer);}
                printf("%i) %s\n",index,buffer);
                fflush(stdout);
                zz++;

                usleep(1000);
                maxwaiketime--;
                ////
                thplan[index].buffer_ready=false;
        }
        thplan[index].is_running = false;
        thplan[index].index = 0;
        return 0;
}

/* Read from audio device and display current buffer. */
void drawindic(float pd, char c)
{
        float posi = 2*pd / (fftw.outlen) - 1;
        float lineS = 0;

        switch ( c )
        {
        case 'b':
                glColor3ub(0, 0, 255);//blue
                lineS = 0.04;
                break;
        case 'g':
                glColor3ub(0, 255, 0);//green
                lineS = 0.08;
                break;
        case 'r':
                glColor3ub(255, 0, 0);//red
                break;
        default:
                break;
        }

        glBegin(GL_QUADS); // debut du dessin
        //14
        //23
        glVertex2f((posi-0.0002), -0.54 - lineS);
        glVertex2f((posi-0.0002), -0.5 - lineS);
        glVertex2f((posi + 0.0002), -0.5 - lineS);
        glVertex2f((posi + 0.0002), -0.54 - lineS);
        glEnd(); // fin du dessin
}

void updateDisplay(void)
{
        float bgcolor[3] = DISPLAY_BACKGROUND_COLOR;
        glClearColor(bgcolor[0], bgcolor[1], bgcolor[2], 1);
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

        if (interaction.update)
        {
                /* Try again *now* if it failed. */
                while (audioRead() < 0);
        }

        if (sound.bufferReady)
        {

                /* The buffer is marked as "full". We can now read it. After the
                 * texture has been updated, the buffer gets marked as "not
                 * ready". */

                /* First, copy the current buffer to the secondary buffer. We
                 * will show that second buffer if the first buffer is not yet
                 * ready. */
                memmove(sound.bufferLast, sound.buffer, sound.bufferSizeFrames * 4);
                /*on copy la ligne pour le traitement*/

                for (int i = 0; i < (int)sound.bufferSizeFrames; i++)
                {
                        short int valq = getFrame(sound.buffer, i, CLEFT);
                        short int vali = getFrame(sound.buffer, i, CRIGHT);
                        fftw.in[i][_Q_]=((2 * (double)vali / (256 * 256)) * (HANNINGWINDOWS[i]));
                        fftw.in[i][_I_]=((2 * (double)valq / (256 * 256)) * (HANNINGWINDOWS[i]));
                }

                fftw_execute(fftw.plan);

                for (int z=0; z<max_index_thread; z++) {
                        if(thplan[z].buffer_ready==false) {
                                thplan[z].buffer_ready=true;
                        }
                }

                /* Draw history into a texture. First, move old texture one line up. */
                memmove(fftw.textureData + (3 * fftw.textureWidth), fftw.textureData,(fftw.textureHeight - 1) * fftw.textureWidth * 3);

                int ha = 0, ta = 0;
                double histramp[][4] = DISPLAY_SPEC_HISTORY_RAMP;
                double moyenne = 0;
                for (int p = 0; p < fftw.outlen; p++)
                {
                        int i = 0;
                        if(p<(fftw.outlen/2)) {i=(fftw.outlen/2)-p;}
                        else{i=(fftw.outlen)-(p-(fftw.outlen/2));}

                        double val = sqrt(fftw.out[i][_Q_] * fftw.out[i][_Q_] + fftw.out[i][_I_] * fftw.out[i][_I_]) / FFTW_SCALE;


                        val = val > 1.0 ? 1.0 : val;

                        /* Save current line for current spectrum. */
                        fftw.currentLine[ha++] = val;

                        moyenne = moyenne + val;

                        /* Find first index where "val" is outside that color
                         * interval. */
                        int colat = 1;
                        while (colat < DISPLAY_SPEC_HISTORY_RAMP_NUM && val > histramp[colat][0]) colat++;

                        colat--;

                        /* Scale "val" into this interval. */
                        double span = histramp[colat + 1][0] - histramp[colat][0];
                        val -= histramp[colat][0];
                        val /= span;

                        /* Interpolate those two colors linearly. */
                        double colnow[3];
                        colnow[0] = histramp[colat][1] * (1 - val) + val * histramp[colat + 1][1];
                        colnow[1] = histramp[colat][2] * (1 - val) + val * histramp[colat + 1][2];
                        colnow[2] = histramp[colat][3] * (1 - val) + val * histramp[colat + 1][3];

                        /* Write this line into new first line of the texture. */
                        fftw.textureData[ta++] = (unsigned char)(colnow[0] * 255);
                        fftw.textureData[ta++] = (unsigned char)(colnow[1] * 255);
                        fftw.textureData[ta++] = (unsigned char)(colnow[2] * 255);
                }

                /*On calcule la puissance moyenne et le squelch*/
                moyenne /= fftw.outlen;

                double threshold = 0;
                threshold = (double)moyenne * threshold_fact;



                /*On cherche les peak > au squelch, max max_index_peak */
                int peakindex[max_index_peak] = {0};
                int lastpe = 0;
                for (int i = start_band; i < stop_band; ++i)
                {
                        if (fftw.currentLine[i] > threshold)
                        {
                                while (fftw.currentLine[i+1] > fftw.currentLine[i]) {i++;} //On cherche le point haut du pic
                                peakindex[lastpe] = i;
                                lastpe++;
                                while (fftw.currentLine[i+1] < fftw.currentLine[i]) {i++;} //On cherche le point bas du pic
                        }
                }


                for (int d = 0; d < lastpe; ++d)
                {
                        drawindic(peakindex[d], 'r');
                }

                for (unsigned int d = 0; d < max_index_peak; ++d)
                {
                        if(peakindex[d]!=0) {
                                int endpeak = d;
                                while ((peakindex[endpeak+1] == ((peakindex[endpeak])+1))) {endpeak++;}
                                d = d + ((endpeak - d)/2);
                                peak[(int)peakindex[d]]+= 40;
                                if(peak[(int)peakindex[d]]>400) peak[(int)peakindex[d]]= 400;
                                d=endpeak;
                        }
                }

                for ( int d = 0; d < SOUND_SAMPLES_PER_TURN; ++d)
                {
                        if(peak[d] > 0) peak[d]--;
                }


                for ( int r = 0; r < SOUND_SAMPLES_PER_TURN; ++r)
                {
                        int d = aleath[r];
                        if(peak[d]>200) {

                                bool is_not_in_array = true;
                                int witch_is_notset = -1;
                                for (int z=0; z<max_index_thread; z++) {
                                        if(thplan[z].index == d) {
                                                is_not_in_array = false;
                                                break;
                                        }
                                        if(thplan[z].is_running == false) {witch_is_notset=z;}
                                }


                                if((witch_is_notset != -1) && (is_not_in_array)) {
                                        thplan[witch_is_notset].is_running = true;
                                        //si un thread est libre
                                        //creation du thread
                                        thplan[witch_is_notset].index = d;

                                        struct thargs *argvvv = (struct thargs *)malloc(sizeof(struct thargs));

                                        argvvv->index=witch_is_notset;
                                        argvvv->fftindexhz = d;

                                        pthread_create(&decoders[witch_is_notset],NULL,&decoderthfc,(void*)argvvv);
                                }
                                if(is_not_in_array) {drawindic(d, 'b');} // bleu
                                else{drawindic(d, 'g');} //vert
                        }
                }
        }

        /* Enable texturing for the quad/history. */
        glEnable(GL_TEXTURE_2D);
        glBindTexture(GL_TEXTURE_2D, fftw.textureHandle);
        if (sound.bufferReady)
        {
                /* Update texture. */
                glTexSubImage2D(GL_TEXTURE_2D, 0, 0, 0,
                                fftw.textureWidth, fftw.textureHeight,
                                GL_RGB, GL_UNSIGNED_BYTE, fftw.textureData);
                checkError(__LINE__);

                /* Reset buffer state. The buffer is no longer ready and we
                 * can't update the texture from it until audioRead() re-marked
                 * it as ready. */
                sound.bufferReady = 0;
        }

        /* Apply zoom and panning. */
        glMatrixMode(GL_MODELVIEW);
        glLoadIdentity();
        if (!interaction.forceOverview)
        {
                glScaled(interaction.scaleX, 1, 1);
                glTranslated(interaction.offsetX, 0, 0);
        }

        /* Draw a textured quad. */
        glColor3f(1, 1, 1);
        glBegin(GL_QUADS);
        /* The texture must be moved half the width of a bin to the left to
         * match the line spectrogram. (Yes, these "0.5"s cancel out. Let
         * the compiler do this. It's easier to understand this way.) */
        double halfBin = (0.5 * fftw.binWidth) / (0.5 * SOUND_RATE);
        glTexCoord2d(0 + halfBin, 0);  glVertex2f(-1, -0.5);
        glTexCoord2d(1 + halfBin, 0);  glVertex2f( 1, -0.5);
        glTexCoord2d(1 + halfBin, 1);  glVertex2f( 1,  1);
        glTexCoord2d(0 + halfBin, 1);  glVertex2f(-1,  1);
        glEnd();
        glDisable(GL_TEXTURE_2D);

        /* Show current spectrum. */
        if (!interaction.showWaveform)
        {
                // float curcol[3] = DISPLAY_SPEC_CURRENT_COLOR;
                // glColor3fv(curcol);

                if(interaction.showthread>=0) {
                        char txt[3];
                        sprintf(txt, "%d", interaction.showthread);
                        glRasterPos2f(0, 0);
                        int len = (int)strlen(txt);
                        for (int i = 0; i < len; i++) {
                                glutBitmapCharacter(GLUT_BITMAP_HELVETICA_12, txt[i]);
                        }

                        glBegin(GL_LINE_STRIP);
                        glColor3ub(0, 255, 0);//green
                        for (int i = 0; i < fftw.outlen; i++)
                        {
                                /* relX will be in [-1, 1], relY will be in [0, 1]. */
                                double relX = 2 * ((double)i / fftw.outlen) - 1;
                                // double relY = fftw.currentLine[i];

                                double relY = thplan[interaction.showthread].currentcolG[i];
                                /* Move relY so it'll be shown at the bottom of the screen. */
                                relY *= 0.5;
                                relY -= 1;
                                glVertex2f(relX, relY);
                        }
                        glEnd();

                        glBegin(GL_LINE_STRIP);
                        glColor3ub(255, 0, 0);//red
                        for (int i = 0; i < fftw.outlen; i++)
                        {
                                /* relX will be in [-1, 1], relY will be in [0, 1]. */
                                double relX = 2 * ((double)i / fftw.outlen) - 1;
                                // double relY = fftw.currentLine[i];
                                double relY = thplan[interaction.showthread].currentTH[i];
                                /* Move relY so it'll be shown at the bottom of the screen. */
                                relY *= 0.5;
                                relY -= 1;
                                glVertex2f(relX, relY);

                        }
                        glEnd();

                        glBegin(GL_LINE_STRIP);
                        glColor3ub(0, 0, 255);//red
                        for (int i = 0; i < fftw.outlen; i++)
                        {
                                /* relX will be in [-1, 1], relY will be in [0, 1]. */
                                double relX = 2 * ((double)i / fftw.outlen) - 1;
                                // double relY = fftw.currentLine[i];
                                double relY = (double)thplan[interaction.showthread].currentdTH[i];
                                /* Move relY so it'll be shown at the bottom of the screen. */
                                relY *= 0.5;
                                relY -= 1;
                                glVertex2f(relX, relY);

                        }
                        glEnd();
                }
                else{
                        glBegin(GL_LINE_STRIP);
                        glColor3ub(0, 255, 0);//green
                        for (int i = 0; i < fftw.outlen; i++)
                        {
                                /* relX will be in [-1, 1], relY will be in [0, 1]. */
                                double relX = 2 * ((double)i / fftw.outlen) - 1;
                                double relY = fftw.currentLine[i];
                                /* Move relY so it'll be shown at the bottom of the screen. */
                                relY *= 0.5;
                                relY -= 1;
                                glVertex2f(relX, relY);
                        }
                        glEnd();

                }


        }
        else
        {
                glPushMatrix();
                glLoadIdentity();
                float curcol[3] = DISPLAY_WAVEFORM_COLOR;
                glColor3fv(curcol);
                glBegin(GL_LINE_STRIP);
                for (int i = 0; i < (int)sound.bufferSizeFrames; i++)
                {
                        /* relX will be in [-1, 1], relY will be in [-s, s] where s
                         * is WAVEFORM_SCALE. */
                        short int val = getFrame(sound.bufferLast, i, CLEFT);
                        double relX = 2 * ((double)i / sound.bufferSizeFrames) - 1;
                        double relY = 2 * WAVEFORM_SCALE * (double)val / (256 * 256);

                        /* Clamp relY ... WAVEFORM_SCALE may be too high. */
                        relY = relY > 1 ? 1 : relY;
                        relY = relY < -1 ? -1 : relY;

                        /* Move relY so it'll be shown at the bottom of the screen. */
                        relY *= 0.25;
                        relY -= 0.75;
                        glVertex2f(relX, relY);
                }
                glEnd();
                glPopMatrix();
        }

        float lineYStart = -1;
        if (interaction.showWaveform)
                lineYStart = -0.5;

        /* Current line and overtones? */
        if (interaction.showOvertones)
        {
                glBegin(GL_LINES);

                /* Crosshair. */
                float colcross[3] = DISPLAY_LINECOLOR_CROSS;
                glColor3fv(colcross);
                glVertex2f(interaction.lastMouseDownEW[0], lineYStart);
                glVertex2f(interaction.lastMouseDownEW[0], 1);

                glColor3fv(colcross);
                glVertex2f(-1, interaction.lastMouseDownEW[1]);
                glVertex2f( 1, interaction.lastMouseDownEW[1]);

                /* Indicate overtones at all multiples of the current frequency
                 * (... this draws unneccssary lines when zoomed in). Don't draw
                 * these lines if they're less than 5 pixels apart. */
                float colover[3] = DISPLAY_LINECOLOR_OVERTONES;
                glColor3fv(colover);
                double nowscale = interaction.forceOverview ? 1 : interaction.scaleX;
                double xInitial = interaction.lastMouseDownEW[0] + 1;
                if (xInitial * interaction.width * nowscale > 5)
                {
                        double x = xInitial * 2;
                        while (x - 1 < 1)
                        {
                                glVertex2f(x - 1, lineYStart);
                                glVertex2f(x - 1, 1);
                                x += xInitial;
                        }
                }

                /* Undertones until two lines are less than 2 pixels apart. */
                double x = xInitial;
                while ((0.5 * x * interaction.width * nowscale)
                       - (0.25 * x * interaction.width * nowscale) > 2)
                {
                        x /= 2;
                        glVertex2f(x - 1, lineYStart);
                        glVertex2f(x - 1, 1);
                }

                glEnd();
        }
        else if (interaction.showMainGrid)
        {
                /* Show "main grid" otherwise. */

                /*glBegin(GL_LINES);

                   float colgrid1[3] = DISPLAY_LINECOLOR_GRID_1;
                   glColor3fv(colgrid1);
                   glVertex2f(0, lineYStart);
                   glVertex2f(0, 1);

                   float colgrid2[3] = DISPLAY_LINECOLOR_GRID_2;
                   glColor3fv(colgrid2);
                   glVertex2f(0.5, lineYStart);
                   glVertex2f(0.5, 1);

                   glVertex2f(-0.5, lineYStart);
                   glVertex2f(-0.5, 1);

                   glEnd();*/
        }

        if (interaction.showFrequency)
        {
                /* Scale from [-1, 1] to [0, fftw.outlen). */
                double t = (interaction.lastMouseDownEW[0] + 1) / 2.0;
                int bin = (int)round(t * fftw.outlen);
                bin = (bin < 0 ? 0 : bin);
                bin = (bin >= fftw.outlen ? fftw.outlen - 1 : bin);

                /* Where exactly is this bin displayed? We want to snap our
                 * guide line to that position. */
                double snapX = ((double)bin / fftw.outlen) * 2 - 1;

                /* SOUND_RATE and SOUND_SAMPLES_PER_TURN determine the "size" of
                 * each "bin" (see calculation of binWidth). Each bin has a size
                 * of some hertz. The i'th bin corresponds to a frequency of i *
                 * <that size> Hz. Note that the resolution is pretty low on
                 * most setups, so it doesn't make any sense to display decimal
                 * places. */
                int freq = (int)(fftw.binWidth * bin);

                /* Draw frequency -- left or right of the guide line. */
                float coltext[3] = DISPLAY_TEXTCOLOR;
                glColor3fv(coltext);

                double nowscale = interaction.forceOverview ? 1 : interaction.scaleX;
                double nowoffX = interaction.forceOverview ? 0 : interaction.offsetX;
                double screenX = (interaction.lastMouseDownEW[0] + nowoffX) * nowscale;

                /* Flipping the label could be done at exactly 50% of the
                 * screen. But we only flip it if it's some pixels away from the
                 * center. */
                if (screenX < -0.25)
                {
                        interaction.frequencyLabelLeft = 1;
                }
                else if (screenX > 0.25)
                {
                        interaction.frequencyLabelLeft = 0;
                }

                char freqstr[256] = "";
                if (interaction.frequencyLabelLeft)
                {
                        glRasterPos2d(snapX, interaction.lastMouseDownEW[1]);
                        snprintf(freqstr, 256, " <- approx. %d Hz", freq);
                }
                else
                {
                        snprintf(freqstr, 256, "approx. %d Hz -> ", freq);
                        glRasterPos2d(snapX - 10 * (double)strlen(freqstr)
                                      / interaction.width / nowscale,
                                      interaction.lastMouseDownEW[1]);
                }

                size_t i;
                for (i = 0; i < strlen(freqstr); i++)
                        glutBitmapCharacter(GLUT_BITMAP_HELVETICA_10, freqstr[i]);

                /* Show guideline for this frequency. */
                float colcross[3] = DISPLAY_LINECOLOR_CROSS;
                glColor3fv(colcross);
                glBegin(GL_LINES);
                glVertex2f(snapX, lineYStart);
                glVertex2f(snapX, 1);
                glEnd();
        }

        /* Separator between current spectrum and history; border. */
        glBegin(GL_LINES);

        float colborder[3] = DISPLAY_LINECOLOR_BORDER;
        glColor3fv(colborder);
        glVertex2f(-1, -0.5);
        glVertex2f( 1, -0.5);

        glVertex2f(-1, lineYStart);
        glVertex2f(-1, 1);

        glVertex2f( 1, lineYStart);
        glVertex2f( 1, 1);

        glVertex2f((((double)start_band/SOUND_SAMPLES_PER_TURN)*2)-1, -1);
        glVertex2f((((double)start_band/SOUND_SAMPLES_PER_TURN)*2)-1, 1);

        glVertex2f((((double)stop_band/SOUND_SAMPLES_PER_TURN)*2)-1, -1);
        glVertex2f((((double)stop_band/SOUND_SAMPLES_PER_TURN)*2)-1, 1);

        glEnd();


        glutSwapBuffers();
}

/* Simple orthographic projection. */
void reshape(int w, int h)
{
        interaction.width = w;
        interaction.height = h;

        glViewport(0, 0, w, h);
        glMatrixMode(GL_PROJECTION);
        glLoadIdentity();
        glOrtho(-1, 1, -1, 1, -4, 4);

        glMatrixMode(GL_MODELVIEW);
        glLoadIdentity();

        glutPostRedisplay();
}

/* Keyboard interaction. */
void keyboard(unsigned char key,
              int x __attribute__((unused)),
              int y __attribute__((unused)))
{
        switch (key)
        {

        case 27:
        case 'q':
                exit(EXIT_SUCCESS);

        case ' ':
                interaction.update = !interaction.update;
                sound.reprepare = 1;
                break;

        case 'r':
                interaction.offsetX = 0;
                interaction.lastOffsetX = 0;
                interaction.scaleX = 1;
                break;

        case 'o':
                interaction.forceOverview = !interaction.forceOverview;
                break;

        case 'j':
                interaction.scaleX *= 2;
                break;

        case 'k':
                interaction.scaleX /= 2;
                break;

        case 'h':
                interaction.offsetX += 0.5 / interaction.scaleX;
                interaction.lastOffsetX = interaction.offsetX;
                break;

        case 'l':
                interaction.offsetX -= 0.5 / interaction.scaleX;
                interaction.lastOffsetX = interaction.offsetX;
                break;

        case 'H':
                interaction.scaleX = 4;
                interaction.offsetX = 0.75;
                interaction.lastOffsetX = interaction.offsetX;
                break;

        case 'g':
                interaction.showMainGrid = !interaction.showMainGrid;
                break;

        case 't':

                if(interaction.showthread+1>=max_index_thread) {
                        interaction.showthread=-1;
                }else{interaction.showthread++;}
                break;


        case 'w':
                interaction.showWaveform = !interaction.showWaveform;
                break;
        }
}


void vSpecial(int key, int x, int y)
{
        switch (key)
        {
        case GLUT_KEY_UP:
                interaction.scaleX *= 2;
                break;

        case GLUT_KEY_DOWN:
                interaction.scaleX /= 2;
                break;

        case GLUT_KEY_LEFT:
                interaction.offsetX += 0.5 / interaction.scaleX;
                interaction.lastOffsetX = interaction.offsetX;
                break;

        case GLUT_KEY_RIGHT:
                interaction.offsetX -= 0.5 / interaction.scaleX;
                interaction.lastOffsetX = interaction.offsetX;
                break;

        default:
                printf(" Autre Touche Speciale\n ");
                break;
        }
}

/* Convert 2D screen coordinates into world coordinates. */
void worldCoord(int *screen, double *world)
{
        world[0] = (double)screen[0] / interaction.width;
        world[1] = (double)screen[1] / interaction.height;

        world[0] *= 2;
        world[1] *= 2;

        world[0] -= 1;
        world[1] -= 1;

        world[1] *= -1;

        /* Panning and scaling only on X axis. */
        if (!interaction.forceOverview)
        {
                world[0] /= interaction.scaleX;
                world[0] -= interaction.lastOffsetX;
        }
}

/* Mouse clicks. */
void mouse(int button, int state, int x, int y)
{
        if (state == GLUT_DOWN)
        {
                /* Save mouse positions for everything but zooming. */
                if (button == GLUT_LEFT_BUTTON || button == GLUT_RIGHT_BUTTON
                    || button == GLUT_MIDDLE_BUTTON)
                {
                        interaction.lastMouseDownBS[0] = x;
                        interaction.lastMouseDownBS[1] = y;
                        worldCoord(interaction.lastMouseDownBS,
                                   interaction.lastMouseDownBW);
                        interaction.lastMouseDownEW[0] = interaction.lastMouseDownBW[0];
                        interaction.lastMouseDownEW[1] = interaction.lastMouseDownBW[1];
                }

                if (button == GLUT_LEFT_BUTTON)
                {
                        interaction.showOvertones = 1;
                }
                else if (button == GLUT_RIGHT_BUTTON && !interaction.forceOverview)
                {
                        interaction.doPanning = 1;
                        interaction.lastOffsetX = interaction.offsetX;
                }
                else if (button == INTERACTION_ZOOM_IN)
                {
                        interaction.scaleX *= INTERACTION_ZOOM_SPEED;
                }
                else if (button == INTERACTION_ZOOM_OUT)
                {
                        interaction.scaleX /= INTERACTION_ZOOM_SPEED;
                        if (interaction.scaleX < 1)
                                interaction.scaleX = 1;
                }
                else if (button == GLUT_MIDDLE_BUTTON)
                {
                        interaction.showFrequency = 1;
                }
        }
        else
        {
                /* Copy new offset if we were panning. */
                if (interaction.doPanning)
                {
                        double dx = interaction.lastMouseDownEW[0]
                                    - interaction.lastMouseDownBW[0];
                        interaction.offsetX = interaction.lastOffsetX + dx;
                        interaction.lastOffsetX = interaction.offsetX;
                }

                interaction.showOvertones = 0;
                interaction.doPanning = 0;
                interaction.showFrequency = 0;
        }
}

/* Mouse movements/drags. */
void motion(int x, int y)
{
        if (!interaction.showOvertones && !interaction.doPanning
            && !interaction.showFrequency)
                return;

        interaction.lastMouseDownES[0] = x;
        interaction.lastMouseDownES[1] = y;
        worldCoord(interaction.lastMouseDownES, interaction.lastMouseDownEW);

        if (interaction.doPanning)
        {
                double dx = interaction.lastMouseDownEW[0]
                            - interaction.lastMouseDownBW[0];
                interaction.offsetX = interaction.lastOffsetX + dx;
        }
}

/* Create the window, set up callbacks and interaction parameters. */
void displayInit(int argc, char *argv[])
{
        interaction.width = DISPLAY_INITIAL_WIDTH;
        interaction.height = DISPLAY_INITIAL_HEIGHT;
        interaction.update = 1;
        interaction.showOvertones = 0;
        interaction.doPanning = 0;
        interaction.forceOverview = 0;
        interaction.showMainGrid = 1;
        interaction.showWaveform = 0;
        interaction.scaleX = 1;
        interaction.offsetX = 0;
        interaction.lastOffsetX = 0;
        interaction.showFrequency = 0;
        interaction.frequencyLabelLeft = 1;

/*#ifdef INTERACTION_ZOOM_STARTUP_FIRST_QUARTER
   interaction.scaleX = 4;
   interaction.offsetX = 0.75;
   interaction.lastOffsetX = interaction.offsetX;
 #endif*/

        glutInit(&argc, argv);
        glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGBA | GLUT_DEPTH);
        glutInitWindowSize(interaction.width, interaction.height);
        glutCreateWindow("HRCW");

        glutDisplayFunc(updateDisplay);
        glutReshapeFunc(reshape);
        glutKeyboardFunc(keyboard);
        glutSpecialFunc(vSpecial);
        glutMouseFunc(mouse);
        glutMotionFunc(motion);
        glutIdleFunc(updateDisplay);
}

/* Create an initial texture (name + data). */
void textureInit(void)
{
        glEnable(GL_TEXTURE_2D);
        glGenTextures(1, &fftw.textureHandle);
        glBindTexture(GL_TEXTURE_2D, fftw.textureHandle);
        glTexImage2D(GL_TEXTURE_2D, 0, 3,
                     fftw.textureWidth, fftw.textureHeight, 0,
                     GL_RGB, GL_UNSIGNED_BYTE, fftw.textureData);
        checkError(__LINE__);

        /* "Smooth" texture filtering. */
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);

        /* No texture wrapping! See display(), we have to move the texture a
         * little to the left. Texture wrapping would result in a wrong
         * spectrogram. */
        glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP);
        glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP);

        glDisable(GL_TEXTURE_2D);
}

/* Delete the texture. */
void textureDeinit(void)
{
        glEnable(GL_TEXTURE_2D);
        glDeleteTextures(1, &fftw.textureHandle);
        glDisable(GL_TEXTURE_2D);
        checkError(__LINE__);
}





int main(int argc, char *argv[])
{
        int c;
        while ((c = getopt (argc, argv, "hd:r:s:e:n:t:")) != -1)
                switch (c)
                {
                case 'h':
                        printf ("HRCW -d plughw:CARD=PCH,DEV=0 -r 192000 -n 10 -s 170000 -e 180000\n");
                        exit(0);
                        break;
                case 'd':
                        SOUND_DEVICE = optarg;
                        break;
                        return 1;
                case 'r':
                        SOUND_RATE = atoi(optarg);
                        break;
                case 's':
                        start_band = atoi(optarg);
                        break;
                case 'e':
                        stop_band = atoi(optarg);
                        break;
                case 'n':
                        max_index_thread = atoi(optarg);
                        break;
                case 't':
                        threshold_fact = atoi(optarg);
                        break;
                        return 1;
                default:
                        abort ();
                }

        SOUND_SAMPLES_PER_TURN = (SOUND_RATE * 512 / 48000);

        if((start_band<0) | (start_band > stop_band) | (stop_band>SOUND_SAMPLES_PER_TURN))
        {
                printf ("Please correct start or stop band in hz\n");
        }

        if(start_band==0) start_band=0;
        if(stop_band==0) stop_band = SOUND_RATE;

        start_band = (start_band*SOUND_SAMPLES_PER_TURN)/SOUND_RATE;
        stop_band = (stop_band*SOUND_SAMPLES_PER_TURN)/SOUND_RATE;

        printf ("SOUND_DEVICE = %s\nSOUND_RATE = %d\nSOUND_SAMPLES_PER_TURN = %d\n",SOUND_DEVICE, SOUND_RATE, SOUND_SAMPLES_PER_TURN);

        HANNINGWINDOWS = hanningInit(SOUND_SAMPLES_PER_TURN);

        peak = (int *) malloc(sizeof(int) * SOUND_SAMPLES_PER_TURN);
        aleath = randomth(SOUND_SAMPLES_PER_TURN);

        displayInit(argc, argv);
        thplan = (thplanst*)malloc(sizeof(thplanst)*max_index_thread);
        decoders = (pthread_t*)malloc(sizeof(pthread_t)*max_index_thread);

        threadinit();
        audioInit();
        fftwInit();
        printf ("FFT resolution: %f",fftw.binWidth);
        interaction.showthread=-1;
        textureInit();

        atexit(audioDeinit);
        atexit(fftwDeinit);
        atexit(textureDeinit);

        glutMainLoop();
        return 0; /* Not reached. */
}
