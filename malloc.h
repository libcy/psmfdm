/*
 * const.h
 *
 *  Created on: Mar 7, 2014
 *      Author: widget
 */

#ifndef MALLOC_H_
#define MALLOC_H_

int* intVec(int n){
	int *vec=(int*)malloc(n*sizeof(int));
	return vec;
}

float* floatVec(int n){
	float *vec=(float*)malloc(n*sizeof(float));
	return vec;
}

int** intMat(int m,int n){
	int **mat=(int **)malloc(m*sizeof(int *));
	for(int i=0;i<m;i++){
		*(mat+i)=(int *)malloc(n*sizeof(int));
	}
	return mat;
}

float** floatMat(int m,int n){
	float **mat=(float **)malloc(m*sizeof(float *));
	float *data=(float *)malloc(m*n*sizeof(float));
	for(int i=0;i<m;i++){
		mat[i]=data+n*i;
	}
	return mat;
}

double** doubleMat(int m,int n){
	double **mat=(double **)malloc(m*sizeof(double *));
	double *data=(double *)malloc(m*n*sizeof(double));
	for(int i=0;i<m;i++){
		mat[i]=data+n*i;
	}
	return mat;
}

void printMat(float **mat,int nx,int ny){
	for(int i=0;i<nx;i++){
		for(int j=0;j<ny;j++){
			// if(abs(mat[i][j])>1e-6)
			printf("%f ", mat[i][j]);
		}
		printf("\n");
	}
}

void printVec(float *vec,int m){
	for(int i=0;i<m;i++){
		printf("%f ",vec[i]);
	}
	printf("\n");
}

void fprintMat(FILE *file,float **mat,int nx,int ny){
	for(int i=0;i<nx;i++){
		for(int j=0;j<ny;j++){
			fprintf(file,"%f\n",mat[i][j]);
		}
	}
}

void fscanfMat(FILE *file,float **mat,int nx,int ny){
	for(int i=0;i<nx;i++){
		for(int j=0;j<ny;j++){
			fscanf(file,"%f",*(mat+i)+j);
		}
	}
}

void fscanfDoubleMat(FILE *file,double **mat,int nx,int ny){
	for(int i=0;i<nx;i++){
		for(int j=0;j<ny;j++){
			fscanf(file,"%lf",*(mat+i)+j);
		}
	}
}

#endif /* MALLOC_H_ */
