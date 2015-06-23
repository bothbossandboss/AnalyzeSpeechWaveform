/**
 * LPC(Linear Predictive Coding, 線形予測分析)を連続で適用し、LPC係数を保存する。
 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <iostream>

#define SAMPLING_F 16000.0

using namespace std;

double *lpc(double r[], double alpha[], double q, int depth, int N){
	//終了条件
	if(depth == N +1){
		return alpha;
	}

	//kの更新
	double sum = 0.0;
	for(int i = 1; i <= depth - 1; i++){
		sum += alpha[i] * r[depth - i];
	}
	double k = -(r[depth] + sum) / q;

	//alphaの更新(1〜depthの配列を有効にするため、確保はdepth+1)
	double *nextAlpha = (double*)calloc(sizeof(double), depth + 1);
	for(int i = 1; i <= depth -1; i++){
		nextAlpha[i] = alpha[i] + k * alpha[depth - i];
	}
	nextAlpha[depth] = k;

	//qの更新
	q = (1.0 - k * k) * q;

	//再帰的に適用。
	lpc(r, nextAlpha, q, depth+1, N);
}

int main(int argc, char *argv[]){
	int N; //LPC予測の点の数。
	int frameLength; //解析するサンプル数(窓長)
	int skipSample; //窓のシフト長。frameLengthの半分
	char inputFile[128];
	char outputFile[128];
	short *rawData;
	double *xr, *rr, *lpcc;
	FILE *input, *output;

    cout << "input file : ";
	cin >> inputFile;
    cout << "output file : ";
	cin >> outputFile;
    cout << "N : ";
	cin >> N;
	cout << "frame length : ";
	cin >> frameLength;
	skipSample = frameLength / 2;
	printf("skip sample : %d\n", skipSample);

	if( (input = fopen(inputFile, "r")) == NULL ){
		perror("open input file");
		return -1;
	}
	if( (output = fopen(outputFile, "w")) == NULL ){
		perror("open output file");
		return -1;
	}
	rawData = (short*)calloc(sizeof(short), frameLength);
	xr = (double*)calloc(sizeof(double), frameLength);
	rr = (double*)calloc( sizeof(double), frameLength);
	lpcc = (double*)calloc( sizeof(double), N+1);  //配列インデックスを1〜Nにするため。
	int r, sum, now;
	sum = now = 0;
	while( (r = fread(rawData, sizeof(short), frameLength, input)) > 0){
		sum += r * (int)sizeof(short);
		printf("%d : r = %d, sum = %d\n", now, r, sum);
		if(r != frameLength) frameLength = r;
		//窓掛け
		for(int n=0;n<frameLength;n++){
			double window = 0.54 - 0.46 * cos(2 * M_PI * n / (frameLength - 1));
			xr[n] = window * rawData[n];
		}
		//自己相関関数
		for(int tau=0;tau<frameLength;tau++){
			for(int i=0;i<frameLength - tau;i++){
				rr[tau] += xr[i] * xr[i+tau];
			}
		}
		//自己相関の正規化
		double r0 = rr[0];
		for(int tau=0;tau<frameLength;tau++){
			rr[tau] = rr[tau] / r0;
		}
		//LPC係数を決定
		double alpha1[2];
		alpha1[1] = - rr[1] / rr[0];
		double k1 = alpha1[1];
		double q1 = rr[0] + rr[1] * alpha1[1];
		lpcc = lpc(rr, alpha1, q1, 2, N);
		//ファイルに係数を出力
		fprintf(output, "%d : (", now);
		for(int i=1;i<N;i++){
			fprintf(output, "%f, ", lpcc[i]);
		}
		fprintf(output, "%f)\n", lpcc[N]);
		//解析開始位置をずらす
		++now;
		int tmp = skipSample * now;
		fseek(input, tmp * sizeof(short), SEEK_SET);
	}

	fclose(input);
	fclose(output);
	return 0;
}