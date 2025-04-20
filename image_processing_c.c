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

#pragma omp parallel for
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

/// Blurs all color channels of an image
void blurIteration(AccurateImage *imageOut, AccurateImage *imageIn, AccuratePixel *summedAreaTable, int size)
{
	// Compute summed area table
	for (int y = 0; y < imageIn->y; y++)
	{
		for (int x = 0; x < imageIn->x; x++)
		{
			int offset = y * imageIn->x + x;
			summedAreaTable[offset].red = imageIn->data[offset].red;
			summedAreaTable[offset].green = imageIn->data[offset].green;
			summedAreaTable[offset].blue = imageIn->data[offset].blue;

			if (x > 0)
			{
				summedAreaTable[offset].red += summedAreaTable[offset - 1].red;
				summedAreaTable[offset].green += summedAreaTable[offset - 1].green;
				summedAreaTable[offset].blue += summedAreaTable[offset - 1].blue;
			}
			if (y > 0)
			{
				summedAreaTable[offset].red += summedAreaTable[offset - imageIn->x].red;
				summedAreaTable[offset].green += summedAreaTable[offset - imageIn->x].green;
				summedAreaTable[offset].blue += summedAreaTable[offset - imageIn->x].blue;
			}
			if (x > 0 && y > 0)
			{
				summedAreaTable[offset].red -= summedAreaTable[offset - imageIn->x - 1].red;
				summedAreaTable[offset].green -= summedAreaTable[offset - imageIn->x - 1].green;
				summedAreaTable[offset].blue -= summedAreaTable[offset - imageIn->x - 1].blue;
			}
		}
	}

	// Compute output using summedAreaTable
	for (int y = 0; y < imageIn->y; y++)
	{
		for (int x = 0; x < imageIn->x; x++)
		{
			int x_start = (x - size < 0) ? 0 : x - size;
			int x_end = (x + size >= imageIn->x) ? imageIn->x - 1 : x + size;
			int y_start = (y - size < 0) ? 0 : y - size;
			int y_end = (y + size >= imageIn->y) ? imageIn->y - 1 : y + size;

			double sum_red = summedAreaTable[y_end * imageIn->x + x_end].red;
			double sum_green = summedAreaTable[y_end * imageIn->x + x_end].green;
			double sum_blue = summedAreaTable[y_end * imageIn->x + x_end].blue;

			if (x_start > 0)
			{
				sum_red -= summedAreaTable[y_end * imageIn->x + (x_start - 1)].red;
				sum_green -= summedAreaTable[y_end * imageIn->x + (x_start - 1)].green;
				sum_blue -= summedAreaTable[y_end * imageIn->x + (x_start - 1)].blue;
			}
			if (y_start > 0)
			{
				sum_red -= summedAreaTable[(y_start - 1) * imageIn->x + x_end].red;
				sum_green -= summedAreaTable[(y_start - 1) * imageIn->x + x_end].green;
				sum_blue -= summedAreaTable[(y_start - 1) * imageIn->x + x_end].blue;
			}
			if (x_start > 0 && y_start > 0)
			{
				sum_red += summedAreaTable[(y_start - 1) * imageIn->x + (x_start - 1)].red;
				sum_green += summedAreaTable[(y_start - 1) * imageIn->x + (x_start - 1)].green;
				sum_blue += summedAreaTable[(y_start - 1) * imageIn->x + (x_start - 1)].blue;
			}

			int count = (x_end - x_start + 1) * (y_end - y_start + 1);
			int offset = y * imageIn->x + x;
			imageOut->data[offset].red = sum_red / count;
			imageOut->data[offset].green = sum_green / count;
			imageOut->data[offset].blue = sum_blue / count;
		}
	}
}

// Helper function to process a single color channel
static inline unsigned char process_channel(double value)
{
	if (value >= 255.0)
		return 255;
	if (value < -1.0)
	{
		value = 257.0 + value;
		return (value >= 255.0) ? 255 : (unsigned char)floor(value);
	}
	if (value < 0.0)
		return 0;
	return (unsigned char)floor(value);
}

PPMImage *imageDifference(const AccurateImage *imageInSmall, const AccurateImage *imageInLarge)
{
	// Allocate output image
	PPMImage *imageOut = malloc(sizeof(PPMImage));

	// Initialize dimensions
	imageOut->x = imageInSmall->x;
	imageOut->y = imageInSmall->y;

	// Allocate pixel data
	size_t pixelCount = imageOut->x * imageOut->y;
	imageOut->data = malloc(pixelCount * sizeof(PPMPixel));

// Check each pixel with OpenMP enabled
#pragma omp parallel for
	for (size_t i = 0; i < pixelCount; i++)
	{
		imageOut->data[i].red = process_channel(
			imageInLarge->data[i].red - imageInSmall->data[i].red);
		imageOut->data[i].green = process_channel(
			imageInLarge->data[i].green - imageInSmall->data[i].green);
		imageOut->data[i].blue = process_channel(
			imageInLarge->data[i].blue - imageInSmall->data[i].blue);
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

	AccurateImage *accurateImages[4][2];

#pragma omp parallel for
	for (int i = 0; i < 4; i++)
	{
		accurateImages[i][0] = convertToAccurateImage(image);
		accurateImages[i][1] = convertToAccurateImage(image);
	}

	const int blurSizes[] = {2, 3, 5, 8};
	const int numBlurs = 4;

#pragma omp parallel for
	for (int i = 0; i < numBlurs; i++)
	{
		AccuratePixel *summedAreaTable = malloc(sizeof(AccuratePixel) * image->x * image->y);

		const int size = blurSizes[i];

		blurIteration(accurateImages[i][1], accurateImages[i][0], summedAreaTable, size);
		blurIteration(accurateImages[i][0], accurateImages[i][1], summedAreaTable, size);
		blurIteration(accurateImages[i][1], accurateImages[i][0], summedAreaTable, size);
		blurIteration(accurateImages[i][0], accurateImages[i][1], summedAreaTable, size);
		blurIteration(accurateImages[i][1], accurateImages[i][0], summedAreaTable, size);

		free(summedAreaTable);
	}

	// calculate difference
	PPMImage *final_tiny = imageDifference(accurateImages[0][1], accurateImages[1][1]);
	PPMImage *final_small = imageDifference(accurateImages[1][1], accurateImages[2][1]);
	PPMImage *final_medium = imageDifference(accurateImages[2][1], accurateImages[3][1]);

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
