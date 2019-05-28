#ifndef __QUANTUM_H__
#define __QUANTUM_H__

#include <iostream>
#include <math.h>
#include <stdio.h>
#include <string.h>
#include <gsl/gsl_sf.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_permutation.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_poly.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_sort.h>
#include <gsl/gsl_sort_vector.h>

#endif

void WignerTransform_au(gsl_matrix_complex * dest, gsl_vector_complex * source, int grades, double dx, double dp);
void WignerTransform_SI(gsl_matrix_complex * dest,gsl_vector_complex * source, int grades, double dx, double dp);