#pragma once
#include <stdio.h>   
#include <stdint.h>
#include <stdbool.h>
void Compression(int s[], int L, int k, double* H);
void Collision(int s[], int L, int k, double* H);
void autocorrelation(int b[], int N, double* P_value);