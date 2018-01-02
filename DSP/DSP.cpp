/*
// DSP.cpp : 콘솔 응용 프로그램에 대한 진입점을 정의합니다.
//

#include "stdafx.h"
#include<stdio.h>
#include<stdlib.h>
#include<math.h>

#define PI 3.141592

#define width  512
#define height 512

void FileOpen(char* fp, char** buf);
void FileWrite(char** buf);
void DFT(char** fbuf);

void FileOpen(char* fp, char **buf) { //파일 열고 버퍼에 정보저장
	FILE *file;
	fopen_s(&file, fp, "rb");
	if (file == NULL) { printf("File open error!!"); return; }
	fread(buf, sizeof(char), width * height, file);
	fclose(file);

	DFT(buf);	
}

void FileWrite(char **buf) { //파일 만들기
	FILE* newfile;
	fopen_s(&newfile, "dft.raw", "wb");
	fwrite(buf, sizeof(char), width*height, newfile);

	fclose(newfile);
}

void DFT(char **fbuf) { 

	printf("DFT Start\n");

	double** temp_re = (double**)malloc(8*width*height); //저장할 공간 할당
	double** temp_im = (double**)malloc(8 * width*height);
//	double** temp_re_2 = (double**)malloc(8 * width*height);
//	double** temp_im_2 = (double**)malloc(8 * width*height);

//	double temp_re[width][height] = { 0, };
//	double temp_im[width][height] = { 0, };
	
//double temp_re_2[width][height] = { 0, };
	//double temp_im_2[width][height] = { 0, };

/*	if (temp_re == NULL || temp_im ==NULL)
	{
		printf("malloc error\n");	
		exit(0);
	}///*
	
	double sum_1, sum_2;
	compl
	for (int j = 0; j < height; j++) {
		for(int i=0;i<width; i++)

			temp_re[j][i] = fbuf[j][i]
	}
	
	for (int j = 0; j < height; j++) {
		for (int i = 0; i < width; i++) {
			sum_1 = sum_2 = 0;
			for (int y = 0; y < height; y++) {
				for (int x = 0; x < width; x++)
				{
					double w = (-2)* PI*((double)i*x / width + (double)j*y / height);
					sum_1 += temp_re[y][x] * cos(w) - temp_im[y][x] * sin(w);
					sum_2 += temp_re[y][x] * sin(w) + temp_im[y][x] * sin(w);
				}
			}
			temp_re[j][i] = sum_1;
			temp_im[j][i] = sum_2;
		}
	}
	/*
	for (int i = 0, sum_1 = 0, sum_2 = 0; i < width; i++) {
		
		for (int j = 0; j < height; j++) {
			//newData[i][j] = 0;
			for (int k = 0; k < height; k++) {
				double w = (-2)* PI*j*k / height;
				sum_1 += (double)fbuf[i][k] * cos(w);
				sum_2 += (double)fbuf[i][k] * sin(w);
			}
			temp_re[i][j] = sum_1;
			temp_im[i][j] = sum_2;
		}
		printf("ROW DFT : %d%%\n", (int)((double)i / width * 100));
	}

	for (int i = 0, sum_1 = 0, sum_2 = 0; i < height; i++) {
		for (int j = 0; j < width; j++) {
			for (int k = 0; k < width; k++) {
				double w = (-2)*PI*j*k / width;

				 sum_1 += temp_re[k][i] * cos(w) - temp_im[k][i] * sin(w);
				sum_2 += temp_re[k][i] * sin(w) + temp_im[k][i] * cos(w);
			}
			temp_re_2[i][j] = sum_1 / sqrt(width);
			temp_im_2[i][j] = sum_2 / sqrt(width);
		}
		printf("Column DFT : %d%%\n", (int)((double)i / width * 100));
	}

	for (int i = 0; i<width; ++i) {
		for (int j = 0; j<height; ++j) {
			double judge = sqrt(pow(temp_re_2[i][j], 2.0) - pow(temp_im_2[i][j], 2.0)) / 100;
			if (judge>0) judge > 255 ? fbuf[i][j] = 255 : fbuf[i][j] = (unsigned char)judge;
		}
		printf("DFT Judge : %d%%\n", (int)((double)i / width * 100));
	}
	//*

	for (int j = 0; j < height; j++)
		for (int i = 0; i < width; i++)
			fbuf[j][i] = (unsigned char)(temp_re[j][i] + temp_im[j][i]);

	FileWrite(fbuf);
	free(temp_re);
	free(temp_im);
//	free(newData);
}

int main()
{
	
	char** fbuf = (char**)malloc(width*height);

	FileOpen("C:\\bear.raw", fbuf);


	free(fbuf);
    return 0;
}
*/
#include"stdafx.h"
#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include <time.h>

#define row 512//행
#define col 512//열
#define PI 3.141592
typedef unsigned char UCHAR;

double arr_real[row][col] = { 0 }, arr_im[row][col] = { 0 };
double arr_temp_1[row][col] = { 0 }, arr_temp_2[row][col] = { 0 };
void DFT(UCHAR(*origin_img_data)[col], UCHAR(*new_img_data)[col]); //DFT function
void IDFT(UCHAR(*new_img_data)[col]); //IDFT function
void LPF(UCHAR(*new_img_data)[col]); //LPF function
void LPF(double(*complex)[col]);
void GaussianNoise(int amount, UCHAR(*origin_data)[col], UCHAR(*new_img_data)[col]);
void fileWrite(char* filename, UCHAR(*write_data)[col]);
void shuffling(UCHAR(*dft_img_data)[col]); //shuffle function
void shuffling(double(*complex)[col]); //shuffle function

int main()
{
	FILE *fp;
	UCHAR origin_data[row][col] = { 0, }, new_img_data[row][col] = { 0, };
	UCHAR noise_data[row][col] = { 0, };
	fopen_s(&fp, "C:\\bear.raw", "rb");
	
	//read
	for (int i = 0; i<col; ++i)	fread(origin_data[i], row, 1, fp); 
	//gaussian
	GaussianNoise(10, origin_data, noise_data);
	fileWrite("result_gaussian.raw", noise_data);
	//dft
	DFT(noise_data, new_img_data);
	fileWrite("result_dft.raw", new_img_data);
	//lpf
	LPF(new_img_data);
	fileWrite("result_lpf.raw", new_img_data);
	//idft
	IDFT(noise_data);
	fileWrite("result_idft.raw", noise_data);

	fclose(fp);

	return 0;
}
void fileWrite(char* filename, UCHAR(*write_data)[col]) {
	FILE* file;
	fopen_s(&file, filename, "wb");

	for (int i = 0; i < col; i++)
	fwrite(write_data[i], row, 1, file);

	fclose(file);
	return;
}
void GaussianNoise(int amount, UCHAR(*origin_data)[col], UCHAR(*new_img_data)[col]) {
	srand((unsigned int)time(NULL));
	double noise;
	double rand1, rand2;
	double gn;
	for (int i = 0; i < row; ++i) {
		for (int j = 0; j < col; ++j) {
			rand1 = (double)rand() / RAND_MAX;
			rand2 = (double)rand() / RAND_MAX;
			gn = amount*sqrt((-2)*log(rand1))*cos(2 * PI*rand2);
			double tmp = origin_data[i][j] + gn;
			if (tmp < 0)
				new_img_data[i][j] = 0;
			else if (tmp > 255)
				new_img_data[i][j] = 255;
			else
				new_img_data[i][j] = (UCHAR)tmp;
		}
		printf("Gaussian : %d%%\n", (int)((double)i / row * 100));
	}
}
void DFT(UCHAR(*origin_img_data)[col], UCHAR(*new_img_data)[col]) //DFT function
{
	double temp_sum_1 = 0, temp_sum_2 = 0, temp = 0;

	//x축
	for (int i = 0; i<row; ++i) {
		for (int j = 0; j<col; ++j) {
			temp_sum_1 = temp_sum_2 = 0;
			for (int k = 0; k<col; ++k) {
				temp_sum_1 += origin_img_data[i][k] * cos((-2)*k*PI*j / col);
				temp_sum_2 += origin_img_data[i][k] * sin((-2)*k*PI*j / col);
			}
			arr_temp_1[i][j] = temp_sum_1;
			arr_temp_2[i][j] = temp_sum_2;
		}
		printf("ROW DFT : %d%%\n", (int)((double)i / row * 100));
	}
	//y축
	for (int i = 0; i<col; ++i) {
		for (int j = 0; j<row; ++j) {
			temp_sum_1 = temp_sum_2 = 0;
			for (int k = 0; k<row; ++k) {
				temp_sum_1 += arr_temp_1[k][i] * cos((-2)*k*PI*j / row) - arr_temp_2[k][i] * sin((-2)*k*PI*j / row);
				temp_sum_2 += arr_temp_1[k][i] * sin((-2)*k*PI*j / row) + arr_temp_2[k][i] * cos((-2)*k*PI*j / row);
			}
			arr_real[i][j] = temp_sum_1;
			arr_im[i][j] = temp_sum_2;
		}
		printf("Column DFT : %d%%\n", (int)((double)i / row * 100));
	}

	for (int i = 0; i<row; ++i) {
		for (int j = 0; j<col; ++j) {
			temp = sqrt(pow(arr_real[i][j], 2.0) - pow(arr_im[i][j], 2.0)) / 100;
			if (temp < 0)
				new_img_data[i][j] = 0;
			else if(temp>255)
				new_img_data[i][j] = 255;
			else
				new_img_data[i][j] = (UCHAR)temp;
		}
	}
	return;
}
void IDFT(UCHAR(*new_img_data)[col]) //IDFT function
{
	shuffling(arr_real); shuffling(arr_im);
	LPF(arr_real); LPF(arr_im);
	shuffling(arr_real); shuffling(arr_im);
	double temp_sum_1 = 0, temp_sum_2 = 0, temp = 0;

	for (int j = 0; j < row; ++j) {
		for (int i = 0; i < col; ++i) {
			temp_sum_1 = 0;
			temp_sum_2 = 0;
			for (int k = 0; k < col; ++k) {
				temp_sum_1 += (arr_real[j][k] * cos((double)(2 * i *PI*k) / col) - arr_im[j][k] * sin((double)(2 * i *PI*k) / col)) / col;
				temp_sum_2 += (arr_real[j][k] * sin((double)(2 * i *PI*k) / col) + arr_im[j][k] * cos((double)(2 * i *PI*k) / col)) / col;
			}
			arr_temp_1[j][i] = temp_sum_1;
			arr_temp_2[j][i] = temp_sum_2;
		}
		printf("Column IDFT : %d%%\n", (int)((double)j / row * 100));
	}
	for (int i = 0; i < col; ++i) {
		for (int j = 0; j < row; ++j) {
			temp_sum_1 = 0;
			temp_sum_2 = 0;
			for (int k = 0; k < row; ++k) {
				temp_sum_1 += (arr_temp_1[k][i] * cos((double)(2 * j  *PI*k) / row) - arr_temp_2[k][i] * sin((double)(2 * j  *PI*k) / row)) / row;
				temp_sum_2 += (arr_temp_1[k][i] * sin((double)(2 * j  *PI*k) / row) + arr_temp_2[k][i] * cos((double)(2 * j  *PI*k) / row)) / row;
			}
			arr_real[i][j] = temp_sum_1;
			arr_im[i][j] = temp_sum_2;
		}
		printf("ROW IDFT : %d%%\n", (int)((double)i / row * 100));
	}
	for (int i = 0; i < row; ++i) {
		for (int j = 0; j < col; ++j) {
			temp = sqrt(pow(arr_real[i][j], 2.0) - pow(arr_im[i][j], 2.0));
			if (temp < 0)
				new_img_data[i][j] = 0;
			else if (temp>255)
				new_img_data[i][j] = 255;
			else
				new_img_data[i][j] = (UCHAR)temp;
		}
	}
	
	return;
}
void LPF(UCHAR(*new_img_data)[col]) //LPF function
{
	for (int i = 0; i < row; ++i)
		for (int j = 0; j < col; ++j)
			if (i<(int)(row * 3 / 8) || i >= (int)(row * 5 / 8) || j<(int)(col * 3 / 8) || j >= (int)(col * 5 / 8))	new_img_data[i][j] = 0;
	return;
}
void LPF(double(*complex)[col]) //LPF function
{
	for (int i = 0; i < row; ++i)
		for (int j = 0; j < col; ++j)
			if (i<(int)(row * 3 / 8) || i >= (int)(row * 5 / 8) || j<(int)(col * 3 / 8) || j >= (int)(col * 5 / 8))	complex[i][j] = 0;
	return;
}
void shuffling(UCHAR(*new_img_data)[col]) //shuffle function
{
	UCHAR arr_temp[row][col] = { 0 };
	for (int i = 0; i<row / 2; ++i) {
		memcpy(arr_temp[i], new_img_data[row / 2 + i], sizeof(UCHAR)*col);
		memcpy(arr_temp[row / 2 + i], new_img_data[i], sizeof(UCHAR)*col);
	}
	for (int i = 0; i<row; ++i)	memcpy(new_img_data[i], arr_temp[i], sizeof(UCHAR)*col);
	for (int i = 0; i<row; ++i) {
		for (int j = 0; j<col / 2; ++j) {
			memcpy(arr_temp[i], new_img_data[i] + (col / 2 - 1), sizeof(UCHAR)*col / 2);
			memcpy(arr_temp[i] + (col / 2 - 1), new_img_data[i], sizeof(UCHAR)*col / 2);
		}
	}
	for (int i = 0; i<row; ++i)	memcpy(new_img_data[i], arr_temp[i], sizeof(UCHAR)*col);
	return;
}

void shuffling(double(*complex)[col]) //shuffle function
{
	for (int i = 0; i<row / 2; ++i) {
		for (int j = 0; j < col; ++j) {
			arr_temp_1[i][j] = complex[row / 2 + i][j];
			arr_temp_1[row / 2 + i][j] = complex[i][j];
		}
	}
	for (int i = 0; i<row; ++i)
		for (int j = 0; j < col; ++j)	complex[i][j] = arr_temp_1[i][j];
	for (int i = 0; i<row; ++i) {
		for (int j = 0; j<col / 2; ++j) {
			arr_temp_1[i][j] = complex[i][j + col / 2];
			arr_temp_1[i][j + col / 2] = complex[i][j];
		}
	}
	for (int i = 0; i<row; ++i)
		for (int j = 0; j < col; ++j)	complex[i][j] = arr_temp_1[i][j];

	return;
}