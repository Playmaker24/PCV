//============================================================================
// Name        : Pcv3.cpp
// Author      : Ronny Haensch
// Version     : 1.0
// Copyright   : -
// Description : Camera calibration
//============================================================================

#include "Pcv3.h"

namespace pcv3 {

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
 * @brief get the conditioning matrix of given points
 * @param the points as matrix
 * @returns the condition matrix (already allocated)
 */
cv::Matx44f getCondition3D(const std::vector<cv::Vec4f> &points)
{
    // TO DO !!! Done
    cv::Scalar t = cv::mean(points);
    float tx = (float) t[0];
    float ty = (float) t[1];
    float tz = (float) t[2];
    float Sx = 0;
    float Sy = 0;
    float Sz = 0;
    size_t N = points.size();
    for (const auto& point : points){
        Sx += cv::abs(point[0]-tx);
        Sy += cv::abs(point[1]-ty);
        Sz += cv::abs(point[2]-tz);
    }
    Sx /= N;
    Sy /= N;
    Sz /= N;
    return cv::Matx44f(1/Sx, 0, 0, -tx/Sx, 
                        0, 1/Sy, 0, -ty/Sy, 
                        0, 0, 1/Sz, -tz/Sz, 
                        0, 0, 0, 1);
}

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
 * @brief Applies a 3D transformation to an array of points
 * @param H Matrix representing the transformation
 * @param points Array of input points, each in homogeneous coordinates
 * @returns Array of transformed objects.
 */
std::vector<cv::Vec4f> applyH_3D_points(const std::vector<cv::Vec4f>& points, const cv::Matx44f &H)
{
    std::vector<cv::Vec4f> result;
    
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
    for (const auto& point : points) {
        result.push_back(H*point);
    }

    return result;
}



/**
 * @brief Define the design matrix as needed to compute projection matrix
 * @param points2D Set of 2D points within the image
 * @param points3D Set of 3D points at the object
 * @returns The design matrix to be computed
 */
cv::Mat_<float> getDesignMatrix_camera(const std::vector<cv::Vec3f>& points2D, const std::vector<cv::Vec4f>& points3D)
{
    // TO DO !!!
    int N = (int) points2D.size();
    cv::Mat_<float> A = cv::Mat_<float>::zeros(2*N, 12);

    for (int i = 0; i < N; i++) {
        cv::Vec3f x2D = points2D[i];
        cv::Vec4f X3D = points3D[i];

        float u = x2D[0];
        float v = x2D[1];
        float w = x2D[2];

        A.at<float>(2*i, 0) = -w * X3D[0];
        A.at<float>(2*i, 1) = -w * X3D[1];
        A.at<float>(2*i, 2) = -w * X3D[2];
        A.at<float>(2*i, 3) = -w * X3D[3];

        A.at<float>(2*i, 8) = u * X3D[0];
        A.at<float>(2*i, 9) = u * X3D[1];
        A.at<float>(2*i, 10) = u * X3D[2];
        A.at<float>(2*i, 11) = u * X3D[3];

        A.at<float>(2*i+1, 4) = -w * X3D[0]; 
        A.at<float>(2*i+1, 5) = -w * X3D[1]; 
        A.at<float>(2*i+1, 6) = -w * X3D[2]; 
        A.at<float>(2*i+1, 7) = -w * X3D[3]; 

        A.at<float>(2*i+1, 8) = v * X3D[0]; 
        A.at<float>(2*i+1, 9) = v * X3D[1]; 
        A.at<float>(2*i+1, 10) = v * X3D[2]; 
        A.at<float>(2*i+1, 11) = v * X3D[3]; 
    }
    
    return A;
}


/**
 * @brief Solve homogeneous equation system by usage of SVD
 * @param A The design matrix
 * @returns The estimated projection matrix
 */
cv::Matx34f solve_dlt_camera(const cv::Mat_<float>& A)
{
    // TO DO !!! Done
    cv::SVD svd(A, cv::SVD::FULL_UV);
    cv::Matx34f P = cv::Matx34f::zeros();
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 4; ++j) {
            P(i, j) = svd.vt.at<float>(11, 4*i+j);
        }
    }
    return P;
}

/**
 * @brief Decondition a projection matrix that was estimated from conditioned point clouds
 * @param T_2D Conditioning matrix of set of 2D image points
 * @param T_3D Conditioning matrix of set of 3D object points
 * @param P Conditioned projection matrix that has to be un-conditioned (in-place)
 */
cv::Matx34f decondition_camera(const cv::Matx33f& T_2D, const cv::Matx44f& T_3D, const cv::Matx34f& P)
{
    // TO DO !!! Done
    return T_2D.inv() * P * T_3D;
}


/**
 * @brief Estimate projection matrix
 * @param points2D Set of 2D points within the image
 * @param points3D Set of 3D points at the object
 * @returns The projection matrix to be computed
 */
cv::Matx34f calibrate(const std::vector<cv::Vec3f>& points2D, const std::vector<cv::Vec4f>& points3D)
{
    // TO DO !!! Done
    cv::Matx33f T2D = getCondition2D(points2D);
    cv::Matx44f T3D = getCondition3D(points3D);

    std::vector<cv::Vec3f> cond2D = applyH_2D(points2D, T2D, GEOM_TYPE_POINT);
    std::vector<cv::Vec4f> cond3D = applyH_3D_points(points3D, T3D);

    cv::Mat_<float> A = getDesignMatrix_camera(cond2D, cond3D);
    cv::Matx34f P = solve_dlt_camera(A);

    return decondition_camera(T2D, T3D, P);
}



/**
 * @brief Extract and prints information about interior and exterior orientation from camera
 * @param P The 3x4 projection matrix, only "input" to this function
 * @param K Matrix for returning the computed internal calibration
 * @param R Matrix for returning the computed rotation
 * @param info Structure for returning the interpretation such as principal distance
 */
void interprete(const cv::Matx34f &P, cv::Matx33f &K, cv::Matx33f &R, ProjectionMatrixInterpretation &info)
{
    // TO DO !!!
    cv::Matx33f M = P.get_minor<3, 3>(0, 0);
    cv::Matx33f R_temp, K_temp;
    RQDecomp3x3(M, K_temp, R_temp);

    for (int i = 0; i < 3; ++i) {
        if (K_temp(i, i) < 0) {
            K_temp.row(i) = K_temp.row(i) * -1;
            R_temp.row(i) = R_temp.row(i) * -1;
        }
    }

    K = (1.0f / K_temp(2,2))*K_temp;
    R = R_temp;

    // Principal distance or focal length
    info.principalDistance = K(0,0);
    
    // Skew as an angle and in degrees
    info.skew = atan(-K(0, 0) / K(0, 1)) * 180.0f / (float) CV_PI;
    
    // Aspect ratio of the pixels
    info.aspectRatio = K(1, 1) / K(0, 0);
    
    // Location of principal point in image (pixel) coordinates
    info.principalPoint(0) = K(0, 2);
    info.principalPoint(1) = K(1, 2);
    
    // Camera rotation angle 1/3
    info.omega = (float) atan2(-R(2, 1), R(2, 2)) * 180.0f / (float) CV_PI;
    
    // Camera rotation angle 2/3
    info.phi = (float) asin(R(2,0)) * 180.0f / (float) CV_PI;
    
    // Camera rotation angle 3/3
    info.kappa = (float) atan2(-R(1, 0), R(0, 0)) * 180.0f / (float) CV_PI;
    
    cv::Matx<float, 3, 1> C = -M.inv() * P.col(3);  
    // Directly assign values from the Matx to the info struct
    info.cameraLocation(0) = C(0, 0);
    info.cameraLocation(1) = C(1, 0);
    info.cameraLocation(2) = C(2, 0);

}




}
