//============================================================================
// Name        : Pcv2test.cpp
// Author      : Ronny Haensch
// Version     : 1.0
// Copyright   : -
// Description : 
//============================================================================

#include "Pcv2.h"

namespace pcv2 {


/**
 * @brief get the conditioning matrix of given points
 * @param the points as matrix
 * @returns the condition matrix (already allocated)
 */
cv::Matx33f getCondition2D(const std::vector<cv::Vec3f> &points)
{
    // TO DO !!! Done
    cv::Scalar t = cv::mean(points);
    float tx = (float) t[0];
    float ty = (float) t[1];
    float Sx = 0;
    float Sy = 0;
    size_t N = points.size();
    for (const auto& point : points){
        Sx += cv::abs(point[0]-tx);
        Sy += cv::abs(point[1]-ty);
    }
    Sx /= N;
    Sy /= N;
    return cv::Matx33f(1/Sx, 0, -tx/Sx, 
                        0, 1/Sy, -ty/Sy, 
                        0, 0, 1);
}


/**
 * @brief define the design matrix as needed to compute 2D-homography
 * @param conditioned_base first set of conditioned points x' --> x' = H * x
 * @param conditioned_attach second set of conditioned points x --> x' = H * x
 * @returns the design matrix to be computed
 */
cv::Mat_<float> getDesignMatrix_homography2D(const std::vector<cv::Vec3f> &conditioned_base, const std::vector<cv::Vec3f> &conditioned_attach)
{
    // TO DO !!! Done
    int N = (int) conditioned_base.size();
    cv::Mat_<float> A = cv::Mat_<float>::zeros(2*N, 9);
    for (int i = 0; i < N; i++) {
        cv::Vec3f x_b = conditioned_base[i];
        cv::Vec3f x = conditioned_attach[i];

        A.at<float>(2*i, 0) = -x_b[2] * x[0];
        A.at<float>(2*i, 1) = -x_b[2] * x[1];
        A.at<float>(2*i, 2) = -x_b[2] * x[2];

        A.at<float>(2*i, 6) = x_b[0] * x[0];
        A.at<float>(2*i, 7) = x_b[0] * x[1];
        A.at<float>(2*i, 8) = x_b[0] * x[2];

        A.at<float>(2*i+1, 3) = -x_b[2] * x[0]; 
        A.at<float>(2*i+1, 4) = -x_b[2] * x[1]; 
        A.at<float>(2*i+1, 5) = -x_b[2] * x[2]; 

        A.at<float>(2*i+1, 6) = x_b[1] * x[0]; 
        A.at<float>(2*i+1, 7) = x_b[1] * x[1]; 
        A.at<float>(2*i+1, 8) = x_b[1] * x[2]; 
    }
    return A;
}


/**
 * @brief solve homogeneous equation system by usage of SVD
 * @param A the design matrix
 * @returns solution of the homogeneous equation system
 */
cv::Matx33f solve_dlt_homography2D(const cv::Mat_<float> &A)
{
    // TO DO !!! Done
    cv::SVD svd(A, cv::SVD::FULL_UV);
    cv::Matx33f H = cv::Matx33f::zeros();

    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            H(i, j) = svd.vt.at<float>(8, 3*i+j);
        }
    }
    return H;
}


/**
 * @brief decondition a homography that was estimated from conditioned point clouds
 * @param T_base conditioning matrix T' of first set of points x'
 * @param T_attach conditioning matrix T of second set of points x
 * @param H conditioned homography that has to be un-conditioned (in-place)
 */
cv::Matx33f decondition_homography2D(const cv::Matx33f &T_base, const cv::Matx33f &T_attach, const cv::Matx33f &H) 
{
    // TO DO !!! Done
    return T_base.inv() * H * T_attach;
}


/**
 * @brief compute the homography
 * @param base first set of points x'
 * @param attach second set of points x
 * @returns homography H, so that x' = Hx
 */
cv::Matx33f homography2D(const std::vector<cv::Vec3f> &base, const std::vector<cv::Vec3f> &attach)
{
    // TO DO !!!
    cv::Matx33f T_base = getCondition2D(base);
    cv::Matx33f T_attach = getCondition2D(attach);

    std::vector<cv::Vec3f> conditioned_points_base = applyH_2D(base, T_base, GEOM_TYPE_POINT);
    std::vector<cv::Vec3f> conditioned_points_attach = applyH_2D(attach, T_attach, GEOM_TYPE_POINT);

    cv::Mat_<float> A = getDesignMatrix_homography2D(conditioned_points_base, conditioned_points_attach);

    cv::Matx33f H = solve_dlt_homography2D(A);

    return decondition_homography2D(T_base, T_attach, H);
}



// Functions from exercise 1
// Reuse your solutions from the last exercise here

/**
 * @brief Applies a 2D transformation to an array of points or lines
 * @param H Matrix representing the transformation
 * @param geomObjects Array of input objects, each in homogeneous coordinates
 * @param type The type of the geometric objects, point or line. All are the same type.
 * @returns Array of transformed objects.
 */
std::vector<cv::Vec3f> applyH_2D(const std::vector<cv::Vec3f>& geomObjects, const cv::Matx33f &H, GeometryType type)
{
    std::vector<cv::Vec3f> result;
    
    /******* Small std::vector cheat sheet ************************************/
    /*
     *   Number of elements in vector:                 a.size()
     *   Access i-th element (reading or writing):     a[i]
     *   Resize array:                                 a.resize(count);
     *   Append an element to an array:                a.push_back(element);
     *     \-> preallocate memory for e.g. push_back:  a.reserve(count);
     */
    /**************************************************************************/

    // TO DO !!! Done

    switch (type) {
        case GEOM_TYPE_POINT: {
            for (const auto& obj : geomObjects) {
                result.push_back(H*obj);
            }
        } break;
        case GEOM_TYPE_LINE: {
            cv::Matx33f H_inv_T = H.inv().t();
            for (const auto& obj : geomObjects) {
                result.push_back(H_inv_T*obj);
            }
        } break;
        default:
            throw std::runtime_error("Unhandled geometry type!");
    }
    
    return result;
}


/**
 * @brief Convert a 2D point from Euclidean to homogeneous coordinates
 * @param p The point to convert (in Euclidean coordinates)
 * @returns The same point in homogeneous coordinates
 */
cv::Vec3f eucl2hom_point_2D(const cv::Vec2f& p)
{
    // TO DO !!! Done
    return cv::Vec3f(p[0], p[1], 1.0f);
}

}
