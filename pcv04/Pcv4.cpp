//============================================================================
// Name        : Pcv4.cpp
// Author      : Ronny Haensch, Andreas Ley
// Version     : 2.0
// Copyright   : -
// Description : Estimation of Fundamental Matrix
//============================================================================

#include "Pcv4.h"

#include <random>
#include <opencv2/features2d.hpp>

using namespace cv;
using namespace std;


namespace pcv4 {
    
/**
 * @brief Applies a 2D transformation to an array of points or lines
 * @param H Matrix representing the transformation
 * @param geomObjects Array of input objects, each in homogeneous coordinates
 * @param type The type of the geometric objects, point or line. All are the same type.
 * @returns Array of transformed objects.
 */
std::vector<cv::Vec3f> applyH_2D(const std::vector<cv::Vec3f>& geomObjects, const cv::Matx33f &H, GeometryType type)
{
    // TO DO !!! Done
    std::vector<cv::Vec3f> result;
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
 * @brief Get the conditioning matrix of given points
 * @param p The points as matrix
 * @returns The condition matrix
 */
cv::Matx33f getCondition2D(const std::vector<cv::Vec3f>& points2D)
{
    // TO DO !!! Done
    cv::Scalar t = cv::mean(points2D);
    float tx = (float) t[0];
    float ty = (float) t[1];
    float Sx = 0;
    float Sy = 0;
    size_t N = points2D.size();
    for (const auto& point : points2D){
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
 * @brief Define the design matrix as needed to compute fundamental matrix
 * @param p1 first set of points
 * @param p2 second set of points
 * @returns The design matrix to be computed
 */
cv::Mat_<float> getDesignMatrix_fundamental(const std::vector<cv::Vec3f>& p1_conditioned, const std::vector<cv::Vec3f>& p2_conditioned)
{
    // TO DO !!! Done
    int N = (int) p1_conditioned.size();
    cv::Mat_<float> A = cv::Mat_<float>::zeros(N, 9);
    for (int i = 0; i < N; i++) {
        cv::Vec3f p1 = p1_conditioned[i];
        cv::Vec3f p2 = p2_conditioned[i];
        A(i, 0) = p1[0] * p2[0];
        A(i, 1) = p1[1] * p2[0];
        A(i, 2) = p1[2] * p2[0];
        A(i, 3) = p1[0] * p2[1];
        A(i, 4) = p1[1] * p2[1];
        A(i, 5) = p1[2] * p2[1];
        A(i, 6) = p1[0] * p2[2];
        A(i, 7) = p1[1] * p2[2];
        A(i, 8) = p1[2] * p2[2];
    }
    return A;
}


/**
 * @brief Solve homogeneous equation system by usage of SVD
 * @param A The design matrix
 * @returns The estimated fundamental matrix
 */
cv::Matx33f solve_dlt_fundamental(const cv::Mat_<float>& A)
{
    // TO DO !!! Done
    cv::SVD svd(A, cv::SVD::FULL_UV);
    cv::Matx33f F = cv::Matx33f::zeros();
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            F(i, j) = svd.vt.at<float>(8, 3*i+j);
        }
    }
    return F;
}


/**
 * @brief Enforce rank of 2 on fundamental matrix
 * @param F The matrix to be changed
 * @return The modified fundamental matrix
 */
cv::Matx33f forceSingularity(const cv::Matx33f& F)
{
    // TO DO !!! Done
    cv::SVD svd(F);
    cv::Vec3f d = svd.w; 
    cv::Matx33f D(d[0], 0, 0,
                  0, d[1], 0,
                  0, 0, 0);
    cv::Matx33f U(svd.u);
    cv::Matx33f Vt(svd.vt);
    return U*D*Vt;
}


/**
 * @brief Decondition a fundamental matrix that was estimated from conditioned points
 * @param T1 Conditioning matrix of set of 2D image points
 * @param T2 Conditioning matrix of set of 2D image points
 * @param F Conditioned fundamental matrix that has to be un-conditioned
 * @return Un-conditioned fundamental matrix
 */
cv::Matx33f decondition_fundamental(const cv::Matx33f& T1, const cv::Matx33f& T2, const cv::Matx33f& F)
{
    // TO DO !!! Done
    return T2.t() * forceSingularity(F) * T1;
}


/**
 * @brief Compute the fundamental matrix
 * @param p1 first set of points
 * @param p2 second set of points
 * @return The estimated fundamental matrix
 */
cv::Matx33f getFundamentalMatrix(const std::vector<cv::Vec3f>& p1, const std::vector<cv::Vec3f>& p2)
{
    // TO DO !!! Done
    cv::Matx33f T_fst = getCondition2D(p1);
    cv::Matx33f T_snd = getCondition2D(p2);

    std::vector<cv::Vec3f> cond2D_fst = applyH_2D(p1, T_fst, GEOM_TYPE_POINT);
    std::vector<cv::Vec3f> cond2D_snd = applyH_2D(p2, T_snd, GEOM_TYPE_POINT);

    cv::Mat_<float> A = getDesignMatrix_fundamental(cond2D_fst, cond2D_snd);
    cv::Matx33f F = solve_dlt_fundamental(A);
    cv::Matx33f F_singular = forceSingularity(F);

    return decondition_fundamental(T_fst, T_snd, F_singular);
}



/**
 * @brief Calculate geometric error of estimated fundamental matrix for a single point pair
 * @details Implement the "Sampson distance"
 * @param p1		first point
 * @param p2		second point
 * @param F		fundamental matrix
 * @returns geometric error
 */
float getError(const cv::Vec3f& p1, const cv::Vec3f& p2, const cv::Matx33f& F)
{
    // TO DO !!! Done
    cv::Vec3f line1 = F * p1; 
    cv::Vec3f line2 = F.t() * p2; 
    float numerator = std::pow(p2.dot(line1), 2);
    float denominator = std::pow(line1[0], 2) + std::pow(line1[1], 2) + std::pow(line2[0], 2) + std::pow(line2[1], 2);
    return numerator / denominator;
}

/**
 * @brief Calculate geometric error of estimated fundamental matrix for a set of point pairs
 * @details Implement the mean "Sampson distance"
 * @param p1		first set of points
 * @param p2		second set of points
 * @param F		fundamental matrix
 * @returns geometric error
 */
float getError(const std::vector<cv::Vec3f>& p1, const std::vector<cv::Vec3f>& p2, const cv::Matx33f& F)
{
    // TO DO !!! Done
    float sum = 0;
    int N = (int) p1.size();
    for (size_t i = 0; i < N; ++i) {
        sum += getError(p1[i], p2[i], F);
    }
    return sum / N; 
}

/**
 * @brief Count the number of inliers of an estimated fundamental matrix
 * @param p1		first set of points
 * @param p2		second set of points
 * @param F		fundamental matrix
 * @param threshold Maximal "Sampson distance" to still be counted as an inlier
 * @returns		Number of inliers
 */
unsigned countInliers(const std::vector<cv::Vec3f>& p1, const std::vector<cv::Vec3f>& p2, const cv::Matx33f& F, float threshold)
{
    // TO DO !!! Done
    unsigned count = 0;
    int N = (int) p1.size();
    for (int i = 0; i < N; i++) {
        float error = getError(p1[i], p2[i], F);
        if (error < threshold) {
            count++;
        }
    }
    return count;
}




/**
 * @brief Estimate the fundamental matrix robustly using RANSAC
 * @details Use the number of inliers as the score
 * @param p1 first set of points
 * @param p2 second set of points
 * @param numIterations How many subsets are to be evaluated
 * @param threshold Maximal "Sampson distance" to still be counted as an inlier
 * @returns The fundamental matrix
 */
cv::Matx33f estimateFundamentalRANSAC(const std::vector<cv::Vec3f>& p1, const std::vector<cv::Vec3f>& p2, unsigned numIterations, float threshold)
{
    const unsigned subsetSize = 8;
    
    std::mt19937 rng;
    std::uniform_int_distribution<unsigned> uniformDist(0, p1.size()-1);
    // Draw a random point index with unsigned index = uniformDist(rng);
    
    // TO DO !!!
    unsigned mostInlier = 0;
    cv::Matx33f bestF = cv::Matx33f::eye();
    for (unsigned i = 0; i < numIterations; i++) {
        std::vector<cv::Vec3f> subset1, subset2;        
        for (int i=0; i  < subsetSize; i++) {
            unsigned idx = uniformDist(rng);
            subset1.push_back(p1[idx]);
            subset2.push_back(p2[idx]);
        }
        cv::Matx33f F = getFundamentalMatrix(subset1, subset2);
        unsigned count = countInliers (p1, p2, F, threshold);
        if (count > mostInlier) {
            mostInlier = count;
            bestF = F;
        }
    }
    return bestF;
}




/**
 * @brief Draw points and corresponding epipolar lines into both images
 * @param img1 Structure containing first image
 * @param img2 Structure containing second image
 * @param p1 First point set (points in first image)
 * @param p2 First point set (points in second image)
 * @param F Fundamental matrix (mapping from point in img1 to lines in img2)
 */
void visualize(const cv::Mat& img1, const cv::Mat& img2, const std::vector<cv::Vec3f>& p1, const std::vector<cv::Vec3f>& p2, const cv::Matx33f& F)
{
    // make a copy to not draw into the original images and destroy them
    cv::Mat img1_copy = img1.clone();
    cv::Mat img2_copy = img2.clone();
        
    // TO DO !!!
    // Compute epilines for both images and draw them with drawEpiLine() into img1_copy and img2_copy respectively
    // Use cv::circle(image, cv::Point2f(x, y), 2, cv::Scalar(0, 255, 0), 2); to draw the points.
    for (int i = 0; i < (int) p1.size(); i++) {
        circle(img1_copy, cv::Point2f(p1[i][0], p1[i][1]), 2, cv::Scalar(0, 255, 0), 2);
        circle(img2_copy, cv::Point2f(p2[i][0], p2[i][1]), 2, cv::Scalar(255, 0, 0), 2);
    }
    for (int i = 0; i < (int) p1.size(); i++) {
        cv::Vec3f line1 = F.t() * p2[i];
        cv::Vec3f line2 = F * p1[i];
        drawEpiLine(img1_copy, line2[0], line2[1], line2[2]);
        drawEpiLine(img2_copy, line1[0], line1[1], line1[2]);
    }
    
    // show images
    cv::imshow("Epilines img1", img1_copy);
    cv::imshow("Epilines img2", img2_copy);
    
    cv::waitKey(0);
}



/**
 * @brief Filters the raw matches
 * @details Applies cross consistency check and ratio test (ratio of 0.75) and returns the point pairs that pass both.
 * @param rawOrbMatches Structure containing keypoints and raw matches obtained from comparing feature descriptors (see Helper.h)
 * @param p1 Points within the first image (returned in the array by this method)
 * @param p2 Points within the second image (returned in the array by this method)
 */
void filterMatches(const RawOrbMatches &rawOrbMatches, std::vector<cv::Vec3f>& p1, std::vector<cv::Vec3f>& p2)
{
    
/******* Small std::map cheat sheet ************************************

// This std::map stores pairs of ints and floats (key value pairs). Each float (value) can quickly be looked up with it's corresponding int (key).
std::map<int, float> exampleMap;
 
// Looking up an element:
int key = 5;
auto it = exampleMap.find(key);
if (it == exampleMap.end()) {
    // no entry with key 5 in the map
} else {
    float value = it->second;
    // do s.th. with the value
}

// Iteration over all elements: 
for (const auto &pair : exampleMap) {
    int key = pair.first;
    float value = pair.second;
}

**************************************************************************/

    p1.clear();
    p2.clear();
    
    const float ratio = 0.75f;

    for (const auto &pair : rawOrbMatches.matches_1_2) {
        
        // TO DO !!!
        // Skip those pairs that don't fulfill the ratio test or cross consistency check
        auto it2 = rawOrbMatches.matches_2_1.find(pair.first);
        if (it2 != rawOrbMatches.matches_2_1.end()) {
            auto &match1 = pair.second;
            auto &match2 = it2->second;
            bool isRatioTestPassed = match1.closestDistance / match1.secondClosestDistance < ratio;
            bool isCrossConsistencyCheckPassed = match2.closest == pair.first;
            if (isRatioTestPassed && isCrossConsistencyCheckPassed) {
                p1.push_back(rawOrbMatches.keypoints1[pair.first]);
                p2.push_back(rawOrbMatches.keypoints2[match1.closest]);
            }
        }
    }
}

/**
 * @brief Computes matches automatically.
 * @details Points will be in homogeneous coordinates.
 * @param img1 The first image
 * @param img2 The second image
 * @param p1 Points within the first image (returned in the array by this method)
 * @param p2 Points within the second image (returned in the array by this method)
 */
void getPointsAutomatic(const cv::Mat &img1, const cv::Mat &img2, std::vector<cv::Vec3f>& p1, std::vector<cv::Vec3f>& p2)
{
    // TO DO !!!
    RawOrbMatches rawOrbMatches = extractRawOrbMatches(img1, img2);
    filterMatches(rawOrbMatches, p1, p2);
}


}
