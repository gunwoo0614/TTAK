#ifndef TEST__H
#define TEST__H

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <stdint.h>
void bit_string_conversion(int L, int l, int s[], int b[]);
void monobit(int b[], int N, double* P_value);
void poker(int b[], int N, double* P_value);
void run(int b[], int N, double* P_value1, double* P_value0);
void long_run(int b[], int N, double* P_value);
void NIST_long_run1(int b[], double* P_value);
void NIST_long_run2(int b[], double* P_value);
void autocorrelation(int b[], int N, double* P_value);

void adaptive_proportion(int s[], int L, int C, double* TF);
void repetition_count(int s[], int L, int C, double* TF);

void The_Most_Common_Value(int s[], int L, int k, double* H);
void Collision(int s[], int L, int k, double* H);
void Markov(int s[], int L, int k, double* H);
void Mutual_Information(int b[], int N, int T, int Q, int K, double* H);
void entropy_test(int b[], int N, int T, int Q, int K, double* H);
void Byte_Correlation(uint8_t X[256][8], int min_distub, int use_shannon, double* e_f);
void Compression(int s[], int L, int k, double* H);







#endif