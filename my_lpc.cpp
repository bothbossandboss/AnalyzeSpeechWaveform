/**
 * LPC(Linear Predictive Coding, 線形予測分析)を実装したもの。
 * スペクトル包絡を求めることを念頭においている。
 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <iostream>

#define OUTPUT "result_lpc.dat"
#define SELF "self.dat"
#define HORAKU "horaku.dat"
#define SAMPLING_F 16000.0

using namespace std;

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

//LPC係数を求める。係数を求めることで、元の波形を予測できる。
/**
 * r[] : 信号の自己相関関数 (r[1]〜r[N]が有効な値。)
 * k -> alpha[] -> qの順に更新
 * depth : 現在解くべきユール・ウォーカー方程式の次元
 * N : 予測に用いる点の数(窓長とは異なることに注意。)
 */
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
	char inputFile[128];
	short *rawData;
	double *xr, *xi, *Xr, *Xi, *rr, *lpcc;
	FILE *input, *output, *self, *horaku;

    cout << "input file : ";
	cin >> inputFile;
    cout << "N : ";
	cin >> N;
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
	if( (self = fopen(SELF, "w")) == NULL ){
		perror("open self file");
		return -1;
	}
	if( (horaku = fopen(HORAKU, "w")) == NULL ){
		perror("open horaku file");
		return -1;
	}
	rawData = (short*)calloc(sizeof(short), frameLength);
	xr = (double*)calloc(sizeof(double), frameLength);
	xi = (double*)calloc(sizeof(double), frameLength);
	Xr = (double*)calloc(sizeof(double), frameLength);
	Xi = (double*)calloc(sizeof(double), frameLength);
	rr = (double*)calloc( sizeof(double), frameLength);
	lpcc = (double*)calloc( sizeof(double), N+1);  //配列インデックスを1〜Nにするため。
	fread(rawData, sizeof(short), frameLength, input);

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
		fprintf(self, "%d %f\n", tau, rr[tau]);
	}

	//LPC係数を決定
	double alpha1[2];
	alpha1[1] = - rr[1] / rr[0];
	double k1 = alpha1[1];
	double q1 = rr[0] + rr[1] * alpha1[1];
	lpcc = lpc(rr, alpha1, q1, 2, N);

	//LPCで推定した音声波形をプロット
	for(int t=0;t<frameLength;t++){
		double hatX = 0.0;
		for(int n=1;n<=N;n++){
			if(t - n < 0) continue;
			hatX += - lpcc[n] * xr[t-n];
		}
		fprintf(output, "%d %f %f\n", t, xr[t], hatX);
	}

	//dftによるスペクトルを求める
	dft(xr, xi, Xr, Xi, frameLength);

	double *fineSpec = (double*)calloc(sizeof(double), frameLength);
	for(int k=0;k<frameLength;k++){
		double power = Xr[k] * Xr[k] + Xi[k] * Xi[k];
		fineSpec[k] = power / (double)frameLength;
		fineSpec[k] = log10(power / (double)frameLength);
	}

	//LPCによる、スペクトル包絡を求め、詳細なスペクトルと共にプロット
	double gg = rr[0];
	for(int n=1;n<=N;n++){
		gg += lpcc[n] * rr[n];
	}

	for(int i=0;i<frameLength;i++){
		double re = 0.0;
		double im = 0.0;
		for(int n=1;n<=N;n++){
			double theata = 2 * M_PI * i * n / frameLength;
			re += lpcc[n] * cos(theata);
			im += lpcc[n] * sin(theata);
		}
		double lower = (re + 1) * (re + 1) + im * im;
		double spec = log10(gg / lower);
		fprintf(horaku, "%d %f %f\n", i, fineSpec[i], spec);
	}

	fclose(input);
	fclose(output);
	fclose(self);
	fclose(horaku);
	return 0;
}