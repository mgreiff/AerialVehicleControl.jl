/*******************************************************************************
* Copyright (C) Marcus Greiff 2020
*
* @file cont_types.c
* @author Marcus Greiff
* @date June 2020
*******************************************************************************/
#include "cont_types.h"

/** \cond INTERNAL */
void dbg_printf(const char *fmt, ...){
  va_list args;
  va_start(args, fmt);
  vfprintf(stderr, fmt, args);
  va_end(args);
}
/** \endcond */

void matrix_print(
  matrix_double_t *mat,
  char * variableName
){
  int i, j, julia = 0;
  for(i = 0; i<strlen(variableName); i++) printf("%c",variableName[i]);
  printf(" = [");
  for(i = 0; i < mat->numRows; i++){
    if(i > 0) for(j = 0; j<strlen(variableName)+4; j++) printf(" ");
    for(j = 0; j < mat->numCols; j++) printf("%f%s", matrix_get(mat, i, j), (j<(mat->numCols - 1))?((julia>0)?" ":", "):((i<(mat->numRows - 1))?";\n":"]\n"));
  }
}

double matrix_get(
  matrix_double_t *mat,
  int row,
  int column
){
  return mat->pData[row + column*mat->numRows];
}

void matrix_set(
  matrix_double_t *mat,
  int row,
  int column,
  double value
){
  mat->pData[row + column*mat->numRows] = value;
}

void matrix_zero(matrix_double_t *matrix){
  int i;
  for (i = 0; i < matrix->numRows*matrix->numCols; i++) matrix->pData[i] = 0.0;
}

void matrix_copy(matrix_double_t *Amat, matrix_double_t *Bmat){
  int i;
  for (i = 0; i < Amat->numRows*Amat->numCols; i++) Bmat->pData[i] = Amat->pData[i];
}

void matrix_identity(matrix_double_t *matrix){
  int i;
  matrix_zero(matrix);
  for (i = 0; i < matrix->numRows; i++) matrix->pData[i*matrix->numRows + i] = 1.0;
}

void matrix_allocate(matrix_double_t *matrix, int numRows, int numCols){
  matrix->numRows = numRows;
  matrix->numCols = numCols;
  matrix->pData   = (double*)malloc(numRows * numCols * sizeof(double));
  if(NULL == matrix->pData){
    TRACE(5,("Could not allocate memory in matrix_allocate()\n"));
    assert(NULL != matrix->pData);
  }
}

void matrix_define(matrix_double_t *matrix, int numRows, int numCols, double *data){
  matrix->numRows = numRows;
  matrix->numCols = numCols;
  matrix->pData   = data;
}
