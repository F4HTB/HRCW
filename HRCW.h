#ifndef CONFIG_H
#define CONFIG_H

/**variables**/

//Limit to start and stop to scan cw
int start_band = 0;
int stop_band  = 0;

//Numer of thread traitement of cw
int max_index_thread = 1;

// Définir la valeur maximale représentable pour un échantillon signé de 16 bits
#define MAX_SAMPLE_VALUE 31103

// Définir le pourcentage maximum de saturation autorisé
#define MAX_SATURATION_PERCENTAGE 0.01

// Définir la valeur RMS minimale acceptable
#define MIN_RMS_THRESHOLD 0.01

/* Default recording device, sample rate. Sound is always recorded in
 * mono. The device can be overriden at runtime using the environment
 * variable SOUND_DEVICE_ENV. */

#define SOUND_DEVICE_ENV "RTSPECCY_CAPTURE_DEVICE"

#define CLEFT 0
#define CRIGHT 2

#define _Q_ 0
#define _I_ 1

#define _LOW_ 0
#define _HIGH_ 1

/* A higher value means larger latency and more memory consumption. On
 * the other hand, a higher value also results in a higher resolution of
 * the spectrogram.
 *
 * Note: This value must be 2^n for some integer n > 1. */


/* Number of lines in the spectrogram history (upper half of the screen).
 *
 * Note: This value must be 2^m for some integer m > 0. */

#define FFTW_HISTORY_SIZE 128

#define max_index_peak 500

/* The fourier transformation produced by libfftw3 is not normalized.
 * Hence, you need to provide a scaling factor. This is a sane default.
 * Increase this value if you see too much noise. */
#define FFTW_SCALE (0.0125 * fftw.outlen)

/* "Home audio" microphones tend to be very quiet. Increase this value
 * if the waveform's amplitude is too low (or adjust your mixer). */
#define WAVEFORM_SCALE 4

/* Initial size of the window. */
#define DISPLAY_INITIAL_WIDTH 1024
#define DISPLAY_INITIAL_HEIGHT 512

/* Background color. Applies to the background of the "current"
 * spectrogram as well as any "exterior" areas. */
#define DISPLAY_BACKGROUND_COLOR { 0, 0, 0 }

/* Color of the "current" spectrogram (bottom of the screen). */
#define DISPLAY_SPEC_CURRENT_COLOR { 0, 1, 0 }

/* Color of the waveform (if shown). */
#define DISPLAY_WAVEFORM_COLOR { 1, 0, 0 }

#define DISPLAY_SPEC_HISTORY_RAMP_NUM 3
/* 2D array: { pos, red, green, blue }. That is, the value 0.1 will get
 * a color between { 0, 0, 0 } and { 0, 0, 1 } because it's less than
 * 0.125 (the "pos" value of the second entry. Whereas a value of 0.8
 * will get a color between { 0, 0, 1 } and { 1, 0, 1 }. Colors are
 * interpolated linearly. */
#define DISPLAY_SPEC_HISTORY_RAMP { \
	{ 0,     0, 0, 0 }, \
	{ 0.125, 0, 0, 1 }, \
	{ 1,     1, 0, 1 }, \
	}

/* Other colors. */
#define DISPLAY_LINECOLOR_CROSS { 0.7, 0, 0 }
#define DISPLAY_LINECOLOR_OVERTONES { 0.7, 0.7, 0.7 }
#define DISPLAY_LINECOLOR_BORDER { 1, 1, 1 }
#define DISPLAY_LINECOLOR_GRID_1 { 1, 1, 1 }
#define DISPLAY_LINECOLOR_GRID_2 { 0.3, 0.3, 0.3 }
#define DISPLAY_TEXTCOLOR { 1, 0, 0 }

/* Zooming. */
#define INTERACTION_ZOOM_SPEED 1.05
#define INTERACTION_ZOOM_IN 4        /* Usually mouse wheel down */
#define INTERACTION_ZOOM_OUT 3       /* Usually mouse wheel up */

/* Most of the interesting stuff (when you're singing or speaking)
 * happens in the first quarter of the spectrogram. We will zoom into
 * that area by default if the following constant is defined. */
#define INTERACTION_ZOOM_STARTUP_FIRST_QUARTER

#endif /* CONFIG_H */

/* vim: set ft=cpp : */
