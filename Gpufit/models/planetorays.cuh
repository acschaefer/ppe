#ifndef GPUFIT_PLANETORAYS_CUH_INCLUDED
#define GPUFIT_PLANETORAYS_CUH_INCLUDED

/* Description of the calculate_planetorays function
* ===================================================
*
* This function calculates the radius of a ray intersecting a plane.
* See rayxpln.h for details.
*
* The parameters are the 3 ranges to define the plane, while the
* corresponding 3 support vectors are provided in the user_info.
* Further, the user_info contains a direction vector for each
* provided data point (radius).
*
* Hence, the user info is composed as:
*     [V1_1, V2_1, V3_1, v1_1, ..., vN_1,
*      ...
*      V1_M, V2_M, V3_M, v1_M, ..., vN_M]
* With:
*     VX : Plane support vectors
*     vX : Data point directional vectors
*     N  : Number of data points per fit = n_points
*     M  : Number of fits = n_fits
*
* Parameters:
*
* parameters     : The 3 ranges defining the plane.
* n_fits         : The number of fits.
* n_points       : The number of data points per fit.
* value          : Output radius of the ray-plane intersection.
* derivative     : An output vector of model function partial derivatives.
* point_index    : The data point index.
* fit_index      : The fit index.
* chunk_index    : The chunk index. Used for indexing of user_info.
* user_info      : An input vector containing user information, see above.
* user_info_size : The size of user_info in bytes.
*
*/

__device__
#include "rayxpln.h"

__device__ void calculate_planetorays(
    REAL const * parameters,
    int const n_fits,
    int const n_points,
    REAL * value,
    REAL * derivative,
    int const point_index,
    int const fit_index,
    int const chunk_index,
    char * user_info,
    std::size_t const user_info_size)
{
    // parameters:  3 plane ranges
    // user info:   3 plane vectors + n_points ray vectors (each vector has 3 REALs)
    // value:       1 range (corresponding to one of the n_points ray vectors)
    
    // indices

    int const v_size = 3;
    int const v_data_size = (3 + n_points) * v_size;
    int const fit_begin = (chunk_index * n_fits + fit_index) * v_data_size;
    REAL * vectors = (REAL*) user_info + fit_begin;

    // value and derivatives

    rayxpln(
      vectors + 0 * v_size,
      vectors + 1 * v_size,
      vectors + 2 * v_size,
      parameters[0],
      parameters[1],
      parameters[2],
      vectors + (3 + point_index) * v_size,
      value[point_index],
      derivative[point_index + 0 * n_points],
      derivative[point_index + 1 * n_points],
      derivative[point_index + 2 * n_points]);
}

#endif
