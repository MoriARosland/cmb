#include <math.h>
#include <string.h>
#include <stdlib.h>

#include <omp.h>

#include "ppm.h"

// Image from:
// http://7-themes.com/6971875-funny-flowers-pictures.html

typedef struct
{
	double red, green, blue;
} AccuratePixel;

typedef struct
{
	int x, y;
	AccuratePixel *data;
} AccurateImage;

// Convert ppm to high precision format.
AccurateImage *convertToAccurateImage(PPMImage *image)
{
	// Make a copy
	AccurateImage *imageAccurate;
	imageAccurate = (AccurateImage *)malloc(sizeof(AccurateImage));
	imageAccurate->data = (AccuratePixel *)malloc(image->x * image->y * sizeof(AccuratePixel));
	for (int i = 0; i < image->x * image->y; i++)
	{
		imageAccurate->data[i].red = (double)image->data[i].red;
		imageAccurate->data[i].green = (double)image->data[i].green;
		imageAccurate->data[i].blue = (double)image->data[i].blue;
	}
	imageAccurate->x = image->x;
	imageAccurate->y = image->y;

	return imageAccurate;
}

PPMImage *convertToPPPMImage(AccurateImage *imageIn)
{
	PPMImage *imageOut;
	imageOut = (PPMImage *)malloc(sizeof(PPMImage));
	imageOut->data = (PPMPixel *)malloc(imageIn->x * imageIn->y * sizeof(PPMPixel));

	imageOut->x = imageIn->x;
	imageOut->y = imageIn->y;

	for (int i = 0; i < imageIn->x * imageIn->y; i++)
	{
		imageOut->data[i].red = imageIn->data[i].red;
		imageOut->data[i].green = imageIn->data[i].green;
		imageOut->data[i].blue = imageIn->data[i].blue;
	}
	return imageOut;
}

// blur one color channel
void blurIteration(AccurateImage *imageOut, AccurateImage *imageIn, int size)
{

	AccuratePixel *kernelIntermediate = malloc(sizeof(AccuratePixel) * imageIn->x * imageIn->y);
	if (kernelIntermediate == NULL)
	{
		fprintf(stderr, "Error allocating memory for kernelIntermediate\n");
		exit(1);
	}

	// Horizontal blur pass
	for (int y = 0; y < imageIn->y; y++)
	{
		for (int x = 0; x < imageIn->x; x++)
		{
			// Discard edges if kernel does not fit
			int start_x = (x - size < 0) ? 0 : x - size;
			int end_x = (x + size >= imageIn->x) ? imageIn->x - 1 : x + size;
			double sum_red = 0.0;
			double sum_green = 0.0;
			double sum_blue = 0.0;

			// Sum pixel values in the horizontal window
			for (int currentX = start_x; currentX <= end_x; currentX++)
			{
				int offset = y * imageIn->x + currentX;

				sum_red += imageIn->data[offset].red;
				sum_green += imageIn->data[offset].green;
				sum_blue += imageIn->data[offset].blue;
			}

			// Compute average and store in intermediate array
			int count = end_x - start_x + 1;
			int offset = y * imageIn->x + x;

			kernelIntermediate[offset].red = sum_red / count;
			kernelIntermediate[offset].green = sum_green / count;
			kernelIntermediate[offset].blue = sum_blue / count;
		}
	}

	// Vertical blur pass
	for (int x = 0; x < imageIn->x; x++)
	{
		for (int y = 0; y < imageIn->y; y++)
		{
			// Discard edges if kernel does not fit
			int start_y = (y - size < 0) ? 0 : y - size;
			int end_y = (y + size >= imageIn->y) ? imageIn->y - 1 : y + size;
			double sum_red = 0.0;
			double sum_green = 0.0;
			double sum_blue = 0.0;

			// Sum pixel values in the vertical window from kernelIntermediate array
			for (int currentY = start_y; currentY <= end_y; currentY++)
			{
				int offset = currentY * imageIn->x + x;
				sum_red += kernelIntermediate[offset].red;
				sum_green += kernelIntermediate[offset].green;
				sum_blue += kernelIntermediate[offset].blue;
			}

			// Compute average and write to output image
			int count = end_y - start_y + 1;
			int offset = y * imageIn->x + x;

			imageOut->data[offset].red = sum_red / count;
			imageOut->data[offset].green = sum_green / count;
			imageOut->data[offset].blue = sum_blue / count;
		}
	}

	free(kernelIntermediate);
}

// Perform the final step, and return it as ppm.
PPMImage *imageDifference(AccurateImage *imageInSmall, AccurateImage *imageInLarge)
{
	PPMImage *imageOut;
	imageOut = (PPMImage *)malloc(sizeof(PPMImage));
	imageOut->data = (PPMPixel *)malloc(imageInSmall->x * imageInSmall->y * sizeof(PPMPixel));

	imageOut->x = imageInSmall->x;
	imageOut->y = imageInSmall->y;

	for (int i = 0; i < imageInSmall->x * imageInSmall->y; i++)
	{
		double value = (imageInLarge->data[i].red - imageInSmall->data[i].red);
		if (value > 255)
			imageOut->data[i].red = 255;
		else if (value < -1.0)
		{
			value = 257.0 + value;
			if (value > 255)
				imageOut->data[i].red = 255;
			else
				imageOut->data[i].red = floor(value);
		}
		else if (value > -1.0 && value < 0.0)
		{
			imageOut->data[i].red = 0;
		}
		else
		{
			imageOut->data[i].red = floor(value);
		}

		value = (imageInLarge->data[i].green - imageInSmall->data[i].green);
		if (value > 255)
			imageOut->data[i].green = 255;
		else if (value < -1.0)
		{
			value = 257.0 + value;
			if (value > 255)
				imageOut->data[i].green = 255;
			else
				imageOut->data[i].green = floor(value);
		}
		else if (value > -1.0 && value < 0.0)
		{
			imageOut->data[i].green = 0;
		}
		else
		{
			imageOut->data[i].green = floor(value);
		}

		value = (imageInLarge->data[i].blue - imageInSmall->data[i].blue);
		if (value > 255)
			imageOut->data[i].blue = 255;
		else if (value < -1.0)
		{
			value = 257.0 + value;
			if (value > 255)
				imageOut->data[i].blue = 255;
			else
				imageOut->data[i].blue = floor(value);
		}
		else if (value > -1.0 && value < 0.0)
		{
			imageOut->data[i].blue = 0;
		}
		else
		{
			imageOut->data[i].blue = floor(value);
		}
	}
	return imageOut;
}

int main(int argc, char **argv)
{
	// read image
	PPMImage *image;
	// select where to read the image from
	if (argc > 1)
	{
		// from file for debugging (with argument)
		image = readPPM("flower.ppm");
	}
	else
	{
		// from stdin for cmb
		image = readStreamPPM(stdin);
	}

	AccurateImage *imageAccurate1_tiny = convertToAccurateImage(image);
	AccurateImage *imageAccurate2_tiny = convertToAccurateImage(image);

	// Process the tiny case:
	int size = 2;
	blurIteration(imageAccurate2_tiny, imageAccurate1_tiny, size);
	blurIteration(imageAccurate1_tiny, imageAccurate2_tiny, size);
	blurIteration(imageAccurate2_tiny, imageAccurate1_tiny, size);
	blurIteration(imageAccurate1_tiny, imageAccurate2_tiny, size);
	blurIteration(imageAccurate2_tiny, imageAccurate1_tiny, size);

	AccurateImage *imageAccurate1_small = convertToAccurateImage(image);
	AccurateImage *imageAccurate2_small = convertToAccurateImage(image);

	// Process the small case:
	size = 3;
	blurIteration(imageAccurate2_small, imageAccurate1_small, size);
	blurIteration(imageAccurate1_small, imageAccurate2_small, size);
	blurIteration(imageAccurate2_small, imageAccurate1_small, size);
	blurIteration(imageAccurate1_small, imageAccurate2_small, size);
	blurIteration(imageAccurate2_small, imageAccurate1_small, size);

	// an intermediate step can be saved for debugging like this
	//    writePPM("imageAccurate2_tiny.ppm", convertToPPPMImage(imageAccurate2_tiny));

	AccurateImage *imageAccurate1_medium = convertToAccurateImage(image);
	AccurateImage *imageAccurate2_medium = convertToAccurateImage(image);

	// Process the medium case:
	size = 5;
	blurIteration(imageAccurate2_medium, imageAccurate1_medium, size);
	blurIteration(imageAccurate1_medium, imageAccurate2_medium, size);
	blurIteration(imageAccurate2_medium, imageAccurate1_medium, size);
	blurIteration(imageAccurate1_medium, imageAccurate2_medium, size);
	blurIteration(imageAccurate2_medium, imageAccurate1_medium, size);

	AccurateImage *imageAccurate1_large = convertToAccurateImage(image);
	AccurateImage *imageAccurate2_large = convertToAccurateImage(image);

	size = 8;
	blurIteration(imageAccurate2_large, imageAccurate1_large, size);
	blurIteration(imageAccurate1_large, imageAccurate2_large, size);
	blurIteration(imageAccurate2_large, imageAccurate1_large, size);
	blurIteration(imageAccurate1_large, imageAccurate2_large, size);
	blurIteration(imageAccurate2_large, imageAccurate1_large, size);

	// calculate difference
	PPMImage *final_tiny = imageDifference(imageAccurate2_tiny, imageAccurate2_small);
	PPMImage *final_small = imageDifference(imageAccurate2_small, imageAccurate2_medium);
	PPMImage *final_medium = imageDifference(imageAccurate2_medium, imageAccurate2_large);
	// Save the images.
	if (argc > 1)
	{
		writePPM("flower_tiny.ppm", final_tiny);
		writePPM("flower_small.ppm", final_small);
		writePPM("flower_medium.ppm", final_medium);
	}
	else
	{
		writeStreamPPM(stdout, final_tiny);
		writeStreamPPM(stdout, final_small);
		writeStreamPPM(stdout, final_medium);
	}
}
