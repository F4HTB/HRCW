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
int threshold_fact = 2;
pthread_t decoders[max_index_thread];
int thplan[max_index_thread][3] = {false,0,false};
pthread_mutex_t mutex_global_thplan;

bool flagk = false;

struct thargs {
    int index;
		double frequency;
};


/* Informations about the window, display options. */
struct interactionInfo
{
	int width;
	int height;

	int update;

	int showOvertones;
	int doPanning;
	int forceOverview;
	int showMainGrid;
	int showWaveform;
	int showFrequency;
	int frequencyLabelLeft;

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
    thplan[z][0]=false;
    thplan[z][1]=0;
    thplan[z][2]=false;

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
	size = sound.bufferSizeFrames * 2 * 2;  /* 2 bytes/sample, 2 channel */

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

    if(itype==1)    //periodic function
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

    if(itype==1)    //periodic function
    {
        for(i=N-1; i>=1; i--)
            w[i] = w[i-1];
        w[0] = 0.0;
    }
    return(w);
}


float **oscInit(int N, float f, float Fs)
{
    float **y;
    float w = 2 * M_PI * f / Fs;

    y = (float**)calloc(N, sizeof(float)*2);
   for(int i = 0; i < N; i++){
      y[i] = (float*)calloc(2, sizeof(float));
      //memset(y[i][_Q_], 0, sizeof(float));
      //memset(y[i][_I_], 0, sizeof(float));
    }



    float dr = cos(w); // dr,di are used to rotate the vector
    float di = sin(w);
    y[0][_Q_] = cos(0); // initial vector at phase phi
    y[0][_I_] = sin(0);

    for (int z=1; z<N; z++) {
         y[z][_Q_] = dr * y[z-1][_Q_] - di * y[z-1][_I_];
         y[z][_I_] = dr * y[z-1][_I_] + di * y[z-1][_Q_];
    }
    return(y);
}


int *randomth(int N)
{
    srand(time(NULL));
    bool intused[N]={0};
    int *number;
    number = (int*) calloc(N, sizeof(int));
    memset(number, 0, N*sizeof(int));

    int p=0;
    for (int i = 0; i < N; i++){
      while(intused[p]){p = rand()%N;}
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

char docodem(char * codem){
    if (codem[0] == '.'){             // .
        if (codem[1] == '.'){         // ..
            if (codem[2] == '.'){     // ...
                if (codem[3] == '.'){ // ....
                    if (strcmp(codem,"...."   ) == 0){ return('H');}
                    if (strcmp(codem,"....."  ) == 0){ return('5'); }
                    if (strcmp(codem,"....-"  ) == 0){ return('4'); }
                } else if (codem[3] == '-'){     // ...-
                    if (codem[4] == '.'){        // ...-.
                        if (strcmp(codem,"...-."  ) == 0)
                                                    { return(126); }
                        if (strcmp(codem,"...-.-" ) == 0)
                                                    { return(62);  }
                        if (strcmp(codem,"...-..-") == 0)
                                                    { return(36);  }
                    } else if (codem[4] == '-'){ // ...--
                        if (strcmp(codem,"...--"  ) == 0)
                                                    { return('3'); }
                    } else {
                        if (strcmp(codem,"...-"   ) == 0)
                                                    { return('V'); }
                    }
                } else {                        // ...
                    if (strcmp(codem,"..."    ) == 0){ return('S'); }
                }
            } else if (codem[2] == '-'){ // ..-
                if (strcmp(codem,"..-"    ) == 0){ return('U');  }
                if (strcmp(codem,"..-."   ) == 0){ return('F');  }
                if (strcmp(codem,"..---"  ) == 0){ return('2');  }
                if (strcmp(codem,"..--.." ) == 0){ return(63);   }
            } else {                    // ..
                if (strcmp(codem,".."      ) == 0){ return('I');  }
            }
        } else if (codem[1] == '-'){         // .-
            if (codem[2] == '.'){            // .-.
                if (codem[3] == '.'){        // .-..
                    if (strcmp(codem,".-.."   ) == 0){ return('L'); }
                    if (strcmp(codem,".-..."  ) == 0){ return(95);  }
                } else if (codem[3] == '-'){ // .-.-
                    if (strcmp(codem,".-.-"   ) == 0){ return(3);   }
                    if (strcmp(codem,".-.-."  ) == 0){ return(60);  }
                    if (strcmp(codem,".-.-.-" ) == 0){ return(46);  }
                } else {                    // .-.
                    if (strcmp(codem,".-."    ) == 0){ return('R'); }
                }
            } else if (codem[2] == '-'){     // .--
                if (codem[3] == '.'){        // .--.
                    if (strcmp(codem,".--."   ) == 0){ return('P'); }
                    if (strcmp(codem,".--.-"  ) == 0){ return(6);   }
                    if (strcmp(codem,".--.-." ) == 0){ return(64);  }
                } else if (codem[3] == '-'){ // .---
                    if (strcmp(codem,".---"   ) == 0){ return('J'); }
                    if (strcmp(codem,".----"  ) == 0){ return('1'); }
                } else {                    // .--
                    if (strcmp(codem,".--"    ) == 0){ return('W'); }
                }
            } else {                        // .-
                if (strcmp(codem,".-") == 0){ return('A'); }
            }
        } else {    // .
            if (strcmp(codem,".") == 0){ return('E'); }
        }
    } else if (codem[0] == '-'){             // -
        if (codem[1] == '.'){                // -.
            if (codem[2] == '.'){            // -..
                if (codem[3] == '.'){        // -...
                    if (strcmp(codem,"-..."   ) == 0){ return('B'); }
                    if (strcmp(codem,"-...."  ) == 0){ return('6'); }
                    if (strcmp(codem,"-....-" ) == 0){ return(45);  }
                } else if (codem[3] == '-'){ // -..-
                    if (strcmp(codem,"-..-"   ) == 0){ return('X'); }
                    if (strcmp(codem,"-..-."  ) == 0){ return(47);  }
                } else {
                    if (strcmp(codem,"-.."    ) == 0){ return('D'); }
                }
            } else if (codem[2] == '-'){     // -.-
                if (codem[3] == '.'){        // -.-.
                    if (strcmp(codem,"-.-."   ) == 0){ return('C'); }
                    if (strcmp(codem,"-.-.--" ) == 0){ return(33);  }
                } else if (codem[3] == '-'){ // -.--
                    if (strcmp(codem,"-.--"   ) == 0){ return('Y'); }
                    if (strcmp(codem,"-.--."  ) == 0){ return(40);  }
                    if (strcmp(codem,"-.--.-" ) == 0){ return(41);  }
                } else {                    // -.-
                    if (strcmp(codem,"-.-"    ) == 0){ return('K'); }
                }
            } else {                        // -.
                if (strcmp(codem,"-.") == 0){ return('N'); }
            }
        } else if (codem[1] == '-'){         // -
            if (codem[2] == '.'){            // --.
                if (strcmp(codem,"--."    ) == 0){ return('G'); }
                if (strcmp(codem,"--.."   ) == 0){ return('Z'); }
                if (strcmp(codem,"--.-"   ) == 0){ return('Q'); }
                if (strcmp(codem,"--..."  ) == 0){ return('7'); }
                if (strcmp(codem,"--..--" ) == 0){ return(44);  }
            } else if (codem[2] == '-'){     // ---
                if (codem[3] == '.'){        // ---.
                    if (strcmp(codem,"---.."  ) == 0){ return('8'); }
                    if (strcmp(codem,"---."   ) == 0){ return(4);   }
                    if (strcmp(codem,"---..." ) == 0){ return(58);  }
                } else if (codem[3] == '-'){ // ----
                    if (strcmp(codem,"----."  ) == 0){ return('9'); }
                    if (strcmp(codem,"-----"  ) == 0){ return('0'); }
                } else {                    // ---
                    if (strcmp(codem,"---"    ) == 0){ return('O'); }
                }
            } else {        // --
                if (strcmp(codem,"--") == 0){ return('M'); }
            }
        } else {    // -
            if (strcmp(codem,"-") == 0){ return('T'); }
        }
    }
}

void* decoderthfc(void *input)
{


    //on parse les infos
  	int index = ((struct thargs*)input)->index;
  	float frequency = ((struct thargs*)input)->frequency;
    printf("new th %i on %f\n", index, frequency);

    int maxwaiketime = 60000; //temps de vie du thread
    long long timebuffer = 0;//sert a avoir une notion de temps en fonction des samples traités

    int frequencyr = 0;
    if(frequency>96000){frequencyr = frequency-96000+1000;}
    else{frequencyr = frequency+96000+1000;}
    printf("realfreq %f\n", frequencyr);

    //on prepare le buffer du thread (I et Q donc *2)
    double bufferth[sound.bufferSizeFrames][2];

  //on calcul la translation en frequence to DC
   float **thxlating;
   thxlating = oscInit(sound.bufferSizeFrames,(frequencyr),SOUND_RATE);

   //on prepare le LPF a 1000 hz
   double RC = 1.0/(1000*2*M_PI);
   double dt = 1.0/192000;
   double alpha = dt/(RC+dt);

   //on calcule la cw
   double      magnitude = 0;
   double      magnitudelimit = 0.0001;
   double      magnitudelimit_low = 0.0001;

   bool        realstate = _LOW_;
   bool        realstatebefore = _LOW_;
   bool        filteredstate = _LOW_;
   bool        filteredstatebefore = _LOW_;

   long     highduration = 0;
   long     laststarttime = 0;
   long     lasthighduration = 0;

   double      q0=0;
   double     nbtime = 10;


   int last_index = 0;
   double acc[2] = {0};

    while (maxwaiketime > 0)
      {
        while(thplan[index][2]==false){
          usleep(1000);
          maxwaiketime--;

        }



        //on transpose la frequence
        for (int i = 0; i < (int)sound.bufferSizeFrames ; i++)
    		{
          short int valq = getFrame(sound.bufferLast, i, CLEFT);
          short int vali = getFrame(sound.bufferLast, i, CRIGHT);
          bufferth[i][_Q_]=((400 * (double)vali / (1 * 256)) * (HANNINGWINDOWS[i]))*thxlating[i][_Q_];
          bufferth[i][_I_]=((400 * (double)valq / (1 * 256)) * (HANNINGWINDOWS[i]))*thxlating[i][_I_];
        }



        //on applique le lpf
        for (int i = 1; i < (int)sound.bufferSizeFrames ; i++)
        {
          bufferth[i][_Q_]=bufferth[i-1][_Q_] + (alpha*(bufferth[i][_Q_] - bufferth[i-1][_Q_]));
          bufferth[i][_I_]=bufferth[i-1][_I_] + (alpha*(bufferth[i][_I_] - bufferth[i-1][_I_]));

        }

        //on downsample a 8khz
        int downfact = SOUND_RATE/8000;
        double bufferthdown[8000][2];
        for (int i = 1; i < (int)(sound.bufferSizeFrames/downfact) ; i++)
        {
          bufferthdown[i][_Q_]=bufferth[i*downfact][_Q_];
          bufferthdown[i][_I_]=bufferth[i*downfact][_I_];
        }




      /*  for (int r = 1; r < (int)(sound.bufferSizeFrames/downfact)/SAMPLEFILTER_TAP_NUM ; r++)
        {
          int index = last_index, i;
          for(i = 0; i < SAMPLEFILTER_TAP_NUM; ++i) {
            index = index != 0 ? index-1 : SAMPLEFILTER_TAP_NUM-1;
            bufferthdown[i*r][_Q_] *= filter_taps[i];
            bufferthdown[i*r][_I_] *= filter_taps[i];
            last_index++;
            if(last_index == SAMPLEFILTER_TAP_NUM)last_index = 0;
          }
        }*/

      for (int c = 1; c < 8 ; c++){
        timebuffer += ((sound.bufferSizeFrames/(downfact*8))*1000)/(SOUND_RATE/downfact);
        for (int z = 1; z < ((int)sound.bufferSizeFrames/downfact)/8 ; z++)
        {
          int i=c*8+z;

          double magnitudeSquared = (bufferthdown[i][_I_]*bufferthdown[i][_I_]) + (bufferthdown[i][_Q_]*bufferthdown[i][_Q_]);
          magnitude = sqrt(magnitudeSquared);
        }
        //printf("magnitude %f\n", magnitude);
        if (magnitude > magnitudelimit){
          magnitudelimit = (magnitudelimit +((magnitude - magnitudelimit)/6));  /// moving average filter
        }

        if (magnitudelimit < magnitudelimit_low)magnitudelimit = magnitudelimit_low;

        if(magnitude > magnitudelimit*0.6)
        {realstate = _HIGH_;}
        else
        {realstate = _LOW_;}
        //printf("magnitude %f magnitudelimit %f\nrealstate %i", magnitude,magnitudelimit,realstate);

        if (realstate != realstatebefore){

            laststarttime = timebuffer;
          }

          if ((timebuffer-laststarttime)> nbtime){
            if (realstate != filteredstate){
              filteredstate = realstate;
            }
          }

          realstatebefore = realstate;
          lasthighduration = highduration;
          filteredstatebefore = filteredstate;
          printf("%i", realstate);
      }

        usleep(1000);
        maxwaiketime--;
        ////
        thplan[index][2]=false;

      }
      free(thxlating);
      thplan[index][0] = false;
      thplan[index][1] = 0;

}

void firdes_lowpass_f(double output[][2], int num_taps, double fc, double fs)
{

  double m_lambda = M_PI * fc / (fs);

  int n;
	double mm;

	for(n = 0; n < num_taps; n++){
		mm = n - (num_taps - 1.0) / 2.0;
		if( mm == 0.0 ){
      output[n][_Q_] = m_lambda / M_PI;
      output[n][_I_] = m_lambda / M_PI;
    }
		else {
      output[n][_Q_] = cos( mm * m_lambda ) / (mm * M_PI);
      output[n][_I_] = sin( mm * m_lambda ) / (mm * M_PI);
    }
  }

	return;

}


void FilterWithFIR2(double FirCoeff[][2], int NumTaps, double Signal[][2], double FilteredSignal[][2], int NumSigPts)
{
 int j, k;
 double y[2], Reg[NumTaps][2];

 for(j=0; j<NumTaps; j++){Reg[j][_I_] = 0.0;Reg[j][_Q_] = 0.0;} // Init the delay registers.

 for(j=0; j<NumSigPts; j++)
 {
  // Shift the register values down and set Reg[0].
  for(k=NumTaps; k>1; k--){
    Reg[k-1][_I_] = Reg[k-2][_I_];
    Reg[k-1][_Q_] = Reg[k-2][_Q_];
  }
  Reg[0][_I_] = Signal[j][_I_];
  Reg[0][_Q_] = Signal[j][_Q_];

  y[_I_] = 0.0;
  y[_Q_] = 0.0;
  for(k=0; k<NumTaps; k++)
  {
    y[_I_] += FirCoeff[k][_I_] * Reg[k][_I_];
    y[_Q_] += FirCoeff[k][_Q_] * Reg[k][_Q_];
  }
  FilteredSignal[j][_I_] = y[_I_];
  FilteredSignal[j][_Q_] = y[_Q_];
 }

}



/* Read from audio device and display current buffer. */
void updateDisplay(void)
{
	//int i;
	float posi = 0;
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

    for (int z=0; z<max_index_thread; z++) {
      if(thplan[z][2]==false){
        thplan[z][2]=true;
      }
    }

    float **xlating;
    xlating = oscInit(sound.bufferSizeFrames,5000,192000);



		for (int i = 0; i < (int)sound.bufferSizeFrames ; i++)
		{
			short int valq = getFrame(sound.buffer, i, CLEFT);
      short int vali = getFrame(sound.buffer, i, CRIGHT);

      if(!flagk){
        fftw.in[i][_Q_]=((2 * (double)vali / (256 * 256)) * (HANNINGWINDOWS[i]));
        fftw.in[i][_I_]=((2 * (double)valq / (256 * 256)) * (HANNINGWINDOWS[i]));
      }
      else{
        vali=0;
        fftw.in[i][_Q_]=((2 * (double)vali / (256 * 256)) * (HANNINGWINDOWS[i]))*xlating[i][_Q_];
        fftw.in[i][_I_]=((2 * (double)valq / (256 * 256)) * (HANNINGWINDOWS[i]))*xlating[i][_I_];
      }
  }

if(flagk){
fftw_complex *fftchg;
fftchg = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * sound.bufferSizeFrames);
FilterWithFIR2(filterth, 500, fftw.in, fftchg, 8192);
memmove(fftw.in, fftchg, sizeof(fftw_complex) * sound.bufferSizeFrames);
}




free(xlating);

		fftw_execute(fftw.plan);



    for (int z=0; z<max_index_thread; z++) {
      if(thplan[z][2]==false){
        thplan[z][2]=true;
      }
    }

		/* Draw history into a texture. First, move old texture one line up. */
		memmove(fftw.textureData + (3 * fftw.textureWidth), fftw.textureData,
		        (fftw.textureHeight - 1) * fftw.textureWidth * 3);

		int ha = 0, ta = 0;
		double histramp[][4] = DISPLAY_SPEC_HISTORY_RAMP;
		double moyenne = 0;
		for (int p = 0; p < fftw.outlen; p++)
		{
      int i = 0;
      if(p<(fftw.outlen/2)){i=(fftw.outlen/2)-p;}
      else{i=(fftw.outlen)-(p-(fftw.outlen/2));}

			double val = sqrt(fftw.out[i][_Q_] * fftw.out[i][_Q_]
			                  + fftw.out[i][_I_] * fftw.out[i][_I_]) / FFTW_SCALE;


			val = val > 1.0 ? 1.0 : val;



			/* Save current line for current spectrum. */
			fftw.currentLine[ha++] = val;

			moyenne = moyenne + val;

			/* Find first index where "val" is outside that color
			 * interval. */
			int colat = 1;
			while (colat < DISPLAY_SPEC_HISTORY_RAMP_NUM
			       && val > histramp[colat][0])
				colat++;

			colat--;

			/* Scale "val" into this interval. */
			double span = histramp[colat + 1][0] - histramp[colat][0];
			val -= histramp[colat][0];
			val /= span;

			/* Interpolate those two colors linearly. */
			double colnow[3];
			colnow[0] = histramp[colat][1] * (1 - val)
			            + val * histramp[colat + 1][1];
			colnow[1] = histramp[colat][2] * (1 - val)
			            + val * histramp[colat + 1][2];
			colnow[2] = histramp[colat][3] * (1 - val)
			            + val * histramp[colat + 1][3];

			/* Write this line into new first line of the texture. */
			fftw.textureData[ta++] = (unsigned char)(colnow[0] * 255);
			fftw.textureData[ta++] = (unsigned char)(colnow[1] * 255);
			fftw.textureData[ta++] = (unsigned char)(colnow[2] * 255);
		}

	/*On calcule la puissance moyenne et le squelch*/
  	moyenne /= fftw.outlen;
		double threshold = (double)moyenne * threshold_fact;

	/*On cherche les peak > au squelch, max max_index_peak */
		int peakindex[max_index_peak] = {0};
		double max = threshold;
		int lastpe = 0;
		for (int i = 0; i < fftw.outlen; ++i)
		{
		    if (fftw.currentLine[i] > max)
		    {
				while (fftw.currentLine[i+1] > fftw.currentLine[i]){i++;} //On cherche le point haut du pic
						peakindex[lastpe] = i;
						lastpe++;
				while (fftw.currentLine[i+1] < fftw.currentLine[i]){i++;} //On cherche le point bas du pic
		    }
		}


	for (int d = 0; d < lastpe; ++d)
	{
		posi = 2*(float)peakindex[d] / (fftw.outlen) -1 ; // divise par deux car complex?
		glColor3ub(255, 0, 0); // rouge
		glBegin(GL_QUADS); // debut du dessin
			//14
			//23
			glVertex2f((posi), -0.54); // vertex 1
			glVertex2f((posi), -0.5);
		    glVertex2f((posi + 0.0005), -0.5); // vertex 3
		    glVertex2f((posi + 0.0005), -0.54); // vertex 4
		glEnd(); // fin du dessin
	}

	for (unsigned int d = 0; d < max_index_peak; ++d)
	{
		if(peakindex[d]!=0){
			int endpeak = d;
			while ((peakindex[endpeak+1] == ((peakindex[endpeak])+1))){endpeak++;}
			d = d + ((endpeak - d)/2);
			peak[(int)peakindex[d]]+= 40;
			if(peak[(int)peakindex[d]]>400) peak[(int)peakindex[d]]= 400;
			d=endpeak;
		}
	}

	for ( int d = 0; d < SOUND_SAMPLES_PER_TURN; ++d)
	{
		if(peak[d] > 0)peak[d]--;
	}


	for ( int r = 0; r < SOUND_SAMPLES_PER_TURN; ++r)
	{
    int d = aleath[r];
		if(peak[d]>200){

      bool is_not_in_array = true;
      int witch_is_notset = -1;
			for (int z=0; z<max_index_thread; z++) {
					  if(thplan[z][1] == d){
                is_not_in_array = false;
						    break;
			       }
             if(thplan[z][0] == false){witch_is_notset=z;}
      }


      if((witch_is_notset != -1) && (is_not_in_array)){
        thplan[witch_is_notset][0] = true;
					//si un thread est libre
					//creation du thread
				  thplan[witch_is_notset][1] = d;

					struct thargs *argvvv = (struct thargs *)malloc(sizeof(struct thargs));

					argvvv->index=witch_is_notset;
					argvvv->frequency = ((float)d * fftw.binWidth );

	        thplan[witch_is_notset][0] = true;

          pthread_create(&decoders[witch_is_notset],NULL,&decoderthfc,(void*)argvvv);

	    }


			posi = 2*(float)d / (fftw.outlen) - 1;
			if(is_not_in_array){glColor3ub(0, 0, 255);} // bleu
      else{glColor3ub(0, 255, 0);} //vert
			glBegin(GL_QUADS); // debut du dessin
				//14
				//23
				glVertex2f((posi), -0.58); // vertex 1
				glVertex2f((posi), -0.55);
			  glVertex2f((posi + 0.0005), -0.55); // vertex 3
			  glVertex2f((posi + 0.0005), -0.58); // vertex 4
			glEnd(); // fin du dessin
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
		float curcol[3] = DISPLAY_SPEC_CURRENT_COLOR;
		glColor3fv(curcol);
		glBegin(GL_LINE_STRIP);
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

		case 'f':
      flagk = !flagk;
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

		case 'w':
			interaction.showWaveform = !interaction.showWaveform;
			break;
	}
}


void vSpecial(int key, int x, int y)
{
	switch (key)
	{
		case GLUT_KEY_UP :
			interaction.scaleX *= 2;
			break;

		case GLUT_KEY_DOWN :
			interaction.scaleX /= 2;
			break;

		case GLUT_KEY_LEFT :
			interaction.offsetX += 0.5 / interaction.scaleX;
			interaction.lastOffsetX = interaction.offsetX;
			break;

		case GLUT_KEY_RIGHT :
			interaction.offsetX -= 0.5 / interaction.scaleX;
			interaction.lastOffsetX = interaction.offsetX;
			break;

		default :
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
	while ((c = getopt (argc, argv, "hd:r:t:")) != -1)
	    switch (c)
	      {
	      case 'h':
		printf ("HRCW -r 192000 -d plughw:CARD=PCH,DEV=0\n");
		exit(0);
		break;
	      case 'd':
		SOUND_DEVICE = optarg;
		break;
		return 1;
	      case 'r':
		SOUND_RATE = atoi(optarg);
		break;
	      case 't':
		threshold_fact = atoi(optarg);
		break;
		return 1;
	      default:
		abort ();
	      }

	SOUND_SAMPLES_PER_TURN = (SOUND_RATE * 2048 / 48000);

	printf ("SOUND_DEVICE = %s\nSOUND_RATE = %d\nSOUND_SAMPLES_PER_TURN = %d\n",SOUND_DEVICE, SOUND_RATE, SOUND_SAMPLES_PER_TURN);

	HANNINGWINDOWS = hanningInit(SOUND_SAMPLES_PER_TURN);

  firdes_lowpass_f(filterth, 100, 5000, SOUND_RATE);


	peak = (int *) malloc(sizeof(int) * SOUND_SAMPLES_PER_TURN);
  aleath = randomth(SOUND_SAMPLES_PER_TURN);

	displayInit(argc, argv);
  threadinit();
	audioInit();
	fftwInit();
	printf ("FFT resolution: %f",fftw.binWidth);
	textureInit();

	atexit(audioDeinit);
	atexit(fftwDeinit);
	atexit(textureDeinit);

	glutMainLoop();
	return 0;  /* Not reached. */
}
