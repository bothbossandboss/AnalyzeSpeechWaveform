#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <iostream>

#define OUTPUT "result_analysis.dat"
#define CEF "cef_result.dat"
#define GOSEI "gosei.dat"
#define SAMPLING_F 16000.0

using namespace std;

/**
 * xr[] : 時間信号の実数部分
 * xi[] : 時間信号の虚数部分
 * Xr[] : 周波数信号の実数部分
 * Xi[] : 周波数信号の虚数部分
 * 	N 	: サンプル数
 */
void dft(double xr[], double xi[], double Xr[], double Xi[], int N){
	for(int k=0;k<N;k++){
		Xr[k] = Xi[k] = 0.0;
		for(int n=0;n<N;n++){
			double theta = 2 * M_PI * n * k / N;
			double Wr = cos(theta);
			double Wi = sin(theta);
			Xr[k] += xr[n] * Wr + xi[n] * Wi;
			Xi[k] += xi[n] * Wr - xr[n] * Wi;
		}
	}
}

void idft(double Xr[], double Xi[], double xr[], double xi[], int N){
	for(int n=0;n<N;n++){
		xr[n] = xi[n] = 0.0;
		for(int k=0;k<N;k++){
			double theta = 2 * M_PI * k * n / N;
			double Wr = cos(theta);
			double Wi = sin(theta);
			xr[n] += Xr[k] * Wr - Xi[k] * Wi;
			xi[n] += Xi[k] * Wr + Xr[k] * Wi;
		}
		xr[n] /= N;
		xi[n] /= N;
	}
}

int main(int argc, char *argv[]){
	double skipSample; //ファイル先頭から解析開始位置をどれだけずらすか。
	int frameLength; //解析するサンプル数(窓長)
	char inputFile[128];
	short *rawData1, *rawData2;
	double *xr, *xi, *Xr, *Xi;
	FILE *input, *output, *cef, *gosei;

    cout << "input file : ";
	cin >> inputFile;
    cout << "skip sample : ";
	cin >> skipSample;
	cout << "frame length : ";
	cin >> frameLength;

	if( (input = fopen(inputFile, "r")) == NULL ){
		perror("open input file");
		return -1;
	}
	if( (output = fopen(OUTPUT, "w")) == NULL ){
		perror("open output file");
		return -1;
	}
	if( (cef = fopen(CEF, "w")) == NULL ){
		perror("open cef file");
		return -1;
	}
	if( (gosei = fopen(GOSEI, "w")) == NULL ){
		perror("open gosei file");
		return -1;
	}
	rawData1 = (short*)calloc(sizeof(short), frameLength);
	rawData2 = (short*)calloc(sizeof(short), frameLength);
	xr = (double*)calloc(sizeof(double), frameLength);
	xi = (double*)calloc(sizeof(double), frameLength);
	Xr = (double*)calloc(sizeof(double), frameLength);
	Xi = (double*)calloc(sizeof(double), frameLength);

	fread(rawData1, sizeof(short), frameLength, input);
	//基本周波数の半分だけずらしてみる。
	printf("%f\n", SAMPLING_F / frameLength);
	skipSample = 4.0 * SAMPLING_F / frameLength;
	printf("%f\n", skipSample);
	fseek(input, skipSample * sizeof(short), SEEK_SET);
	fread(rawData2, sizeof(short), frameLength, input);

	/**
	 * 切り取ったサンプル群が周期的に繰り返される信号だと捉えるため、
	 * フレーム両端の重みを減少させるために窓関数を掛け合わせる。
	 * 今回はハミング窓をかける。
	 */
	for(int n=0;n<frameLength;n++){
		double window = 0.54 - 0.46 * cos(2 * M_PI * n / (frameLength - 1));
		xr[n] = window * rawData1[n];
		xi[n] = 0.0;
	}

	dft(xr, xi, Xr, Xi, frameLength);

	/**
     * 対数パワースペクトルを出力
	 */
	double *Spec1 = (double*)calloc(sizeof(double), frameLength);
	for(int k=0;k<frameLength;k++){
		double power = Xr[k] * Xr[k] + Xi[k] * Xi[k];
		double spectrum = log10(power / (double)frameLength);
		fprintf(output, "%f %f\n", k*SAMPLING_F / frameLength, spectrum);
		Spec1[k] = spectrum;
	}

	/*2回目*/
	for(int n=0;n<frameLength;n++){
		double window = 0.54 - 0.46 * cos(2 * M_PI * n / (frameLength - 1));
		xr[n] = window * rawData2[n];
		xi[n] = 0.0;
	}
	dft(xr, xi, Xr, Xi, frameLength);
	double *Spec2 = (double*)calloc(sizeof(double), frameLength);
	for(int k=0;k<frameLength;k++){
		double power = Xr[k] * Xr[k] + Xi[k] * Xi[k];
		double spectrum = log10(power / (double)frameLength);
		fprintf(cef, "%f %f\n", k*SAMPLING_F / frameLength, spectrum);
		Spec2[k] = spectrum;
	}
	for(int k=0;k<frameLength;k++){
		double spectrum = Spec1[k] + Spec2[k];
		fprintf(gosei, "%f %f\n", k*SAMPLING_F / frameLength, spectrum);
	}

	fclose(input);
	fclose(output);
	fclose(cef);
	fclose(gosei);
	return 0;
}