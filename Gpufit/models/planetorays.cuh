#ifndef GPUFIT_PLANETORAYS_CUH_INCLUDED
#define GPUFIT_PLANETORAYS_CUH_INCLUDED

/* Description of the calculate_linear1d function
* ===================================================
*
* This function calculates the values of one-dimensional linear model functions
* and their partial derivatives with respect to the model parameters. 
*
* This function makes use of the user information data to pass in the 
* independent variables (X values) corresponding to the data.  The X values
* must be of type REAL.
*
* Note that if no user information is provided, the (X) coordinate of the 
* first data value is assumed to be (0.0).  In this case, for a fit size of 
* M data points, the (X) coordinates of the data are simply the corresponding 
* array index values of the data array, starting from zero.
*
* There are three possibilities regarding the X values:
*
*   No X values provided: 
*
*       If no user information is provided, the (X) coordinate of the 
*       first data value is assumed to be (0.0).  In this case, for a 
*       fit size of M data points, the (X) coordinates of the data are 
*       simply the corresponding array index values of the data array, 
*       starting from zero.
*
*   X values provided for one fit:
*
*       If the user_info array contains the X values for one fit, then 
*       the same X values will be used for all fits.  In this case, the 
*       size of the user_info array (in bytes) must equal 
*       sizeof(REAL) * n_points.
*
*   Unique X values provided for all fits:
*
*       In this case, the user_info array must contain X values for each
*       fit in the dataset.  In this case, the size of the user_info array 
*       (in bytes) must equal sizeof(REAL) * n_points * nfits.
*
* Parameters:
*
* parameters: An input vector of model parameters.
*             p[0]: offset
*             p[1]: slope
*
* n_fits: The number of fits.
*
* n_points: The number of data points per fit.
*
* value: An output vector of model function values.
*
* derivative: An output vector of model function partial derivatives.
*
* point_index: The data point index.
*
* fit_index: The fit index.
*
* chunk_index: The chunk index. Used for indexing of user_info.
*
* user_info: An input vector containing user information.
*
* user_info_size: The size of user_info in bytes.
*
* Calling the calculate_linear1d function
* =======================================
*
* This __device__ function can be only called from a __global__ function or an other
* __device__ function.
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
