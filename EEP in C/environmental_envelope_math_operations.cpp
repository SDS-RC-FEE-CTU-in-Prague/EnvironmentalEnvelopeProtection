//
// Created by efrem on 10/30/2023.
//

#include <iostream>
#include <vector>
#include <casadi/casadi.hpp>

// Function to multiply two matrices
std::vector<std::vector<double>>
multiplyMatrices(const std::vector<std::vector<double>> &A, const std::vector<std::vector<double>> &B) {
    unsigned long rows_A = A.size();
    unsigned long cols_A = A[0].size();
    unsigned long cols_B = B[0].size();

    // Create a new matrix to store the result
    std::vector<std::vector<double>> result(rows_A, std::vector<double>(cols_B, 0));

    for (unsigned long i = 0; i < rows_A; i++) {
        for (unsigned long j = 0; j < cols_B; j++) {
            for (unsigned long k = 0; k < cols_A; k++) {
                result[i][j] += A[i][k] * B[k][j];
            }
        }
    }

    return result;
}

// Function to transpose a matrix
std::vector<std::vector<double>> transposeMatrix(const std::vector<std::vector<double>> &matrix) {
    unsigned long rows = matrix.size();
    unsigned long cols = matrix[0].size();

    // Create a new matrix with dimensions swapped
    std::vector<std::vector<double>> result(cols, std::vector<double>(rows));

    for (unsigned long i = 0; i < rows; i++) {
        for (unsigned long j = 0; j < cols; j++) {
            result[j][i] = matrix[i][j];
        }
    }

    return result;
}

// Function to calculate the inverse of a 4x4 matrix
std::vector<std::vector<double>> inverse4x4(const std::vector<std::vector<double>> &matrix) {
    // Check if the matrix is square and 4x4
    if (matrix.size() != 4 || matrix[0].size() != 4) {
        std::cerr << "Error: Input matrix must be a 4x4 matrix." << std::endl;
        return {{0.0, 0.0, 0.0, 0.0},
                {0.0, 0.0, 0.0, 0.0},
                {0.0, 0.0, 0.0, 0.0},
                {0.0, 0.0, 0.0, 0.0}}; // Return an empty matrix to indicate an error
    }

    // Augment the matrix with the identity matrix
    std::vector<std::vector<double>> augmentedMatrix(4, std::vector<double>(8, 0.0));
    for (int i = 0; i < 4; ++i) {
        for (int j = 0; j < 4; ++j) {
            augmentedMatrix[i][j] = matrix[i][j];
            augmentedMatrix[i][j + 4] = (i == j) ? 1.0 : 0.0;
        }
    }

    // Perform Gauss-Jordan elimination
    for (int i = 0; i < 4; ++i) {
        // Make the diagonal element 1
        double diagElement = augmentedMatrix[i][i];
        if (diagElement == 0.0) {
            std::cerr << "Error: Matrix is singular, inverse does not exist." << std::endl;
            return {{0.0, 0.0, 0.0, 0.0},
                    {0.0, 0.0, 0.0, 0.0},
                    {0.0, 0.0, 0.0, 0.0},
                    {0.0, 0.0, 0.0, 0.0}}; // Return an empty matrix to indicate an error
        }
        for (int j = 0; j < 8; ++j) {
            augmentedMatrix[i][j] /= diagElement;
        }

        // Make the other elements in the column 0
        for (int k = 0; k < 4; ++k) {
            if (k != i) {
                double factor = augmentedMatrix[k][i];
                for (int j = 0; j < 8; ++j) {
                    augmentedMatrix[k][j] -= factor * augmentedMatrix[i][j];
                }
            }
        }
    }

    std::vector<std::vector<double>> result(4, std::vector<double>(4));
    // Extract the inverse matrix from the augmented matrix
    for (int i = 0; i < 4; ++i) {
        for (int j = 0; j < 4; ++j) {
            result[i][j] = augmentedMatrix[i][j + 4];
        }
    }

    return result;
}


// Function to calculate the inverse of a 3x3 matrix
std::vector<std::vector<double>> inverse3x3(const std::vector<std::vector<double>> &matrix) {
    if (matrix.size() != 3 || matrix[0].size() != 3 || matrix[1].size() != 3 || matrix[2].size() != 3) {
        std::cerr << "Error: Input matrix is not 3x3." << std::endl;
        return {{0.0, 0.0, 0.0},
                {0.0, 0.0, 0.0},
                {0.0, 0.0, 0.0}}; // Return an empty matrix to indicate an error
    }

    double a = matrix[0][0];
    double b = matrix[0][1];
    double c = matrix[0][2];
    double d = matrix[1][0];
    double e = matrix[1][1];
    double f = matrix[1][2];
    double g = matrix[2][0];
    double h = matrix[2][1];
    double i = matrix[2][2];

    double det = a * (e * i - f * h) - b * (d * i - f * g) + c * (d * h - e * g);

    if (det == 0.0) {
        std::cerr << "Error: The determinant is zero, and the matrix is not invertible." << std::endl;
        return {{0.0, 0.0, 0.0},
                {0.0, 0.0, 0.0},
                {0.0, 0.0, 0.0}}; // Return an empty matrix to indicate an error
    }

    double invDet = 1.0 / det;

    std::vector<std::vector<double>> result(3, std::vector<double>(3));

    result[0][0] = (e * i - f * h) * invDet;
    result[0][1] = (c * h - b * i) * invDet;
    result[0][2] = (b * f - c * e) * invDet;
    result[1][0] = (f * g - d * i) * invDet;
    result[1][1] = (a * i - c * g) * invDet;
    result[1][2] = (c * d - a * f) * invDet;
    result[2][0] = (d * h - e * g) * invDet;
    result[2][1] = (b * g - a * h) * invDet;
    result[2][2] = (a * e - b * d) * invDet;

    return result;
}

// Function to square each element of a vector
std::vector<double> squareElements(const std::vector<double> &input) {
    std::vector<double> result;

    for (const double &value: input) {
        result.push_back(value * value);
    }

    return result;
}

// Function to square each element of a vector
std::vector<double> cubeElements(const std::vector<double> &input) {
    std::vector<double> result;

    for (const double &value: input) {
        result.push_back(value * value * value);
    }

    return result;
}

// Function to multiply a matrix by a vector
std::vector<double>
matrixVectorMultiply(const std::vector<std::vector<double>> &matrix, const std::vector<double> &vector) {
    unsigned long numRows = matrix.size();
    unsigned long numCols = matrix[0].size();

    if (numCols != vector.size()) {
        std::cerr << "Error: Matrix columns must match vector size for multiplication." << std::endl;
        return std::vector<double>(); // Return an empty vector to indicate an error
    }

    std::vector<double> result(numRows, 0.0);

    for (unsigned long i = 0; i < numRows; i++) {
        for (unsigned long j = 0; j < numCols; j++) {
            result[i] += matrix[i][j] * vector[j];
        }
    }

    return result;
}

// fit points with a parabolic curve
void calculateParabolicCoefficients(const std::vector<double> &x, const std::vector<double> &y, double &a0, double &a1,
                                    double &a2) {

    const std::vector<std::vector<double>> &tmp = {{1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0}, x,
                                                   squareElements(x)};
    const std::vector<std::vector<double>> &X = transposeMatrix(tmp);
    const std::vector<std::vector<double>> &M = multiplyMatrices(inverse3x3(multiplyMatrices(transposeMatrix(X), X)),
                                                                 transposeMatrix(X));
    const std::vector<double> &a = matrixVectorMultiply(M, y);

    a0 = a[0];
    a1 = a[1];
    a2 = a[2];

}

void
calculateCubicParabolicCoefficients(const std::vector<double> &x, const std::vector<double> &y, double &a0, double &a1,
                                    double &a2, double &a3) {

    const std::vector<std::vector<double>> &tmp = {{1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0}, x,
                                                   squareElements(x), cubeElements(x)};
    const std::vector<std::vector<double>> &X = transposeMatrix(tmp);
    const std::vector<std::vector<double>> &M = multiplyMatrices(inverse4x4(multiplyMatrices(transposeMatrix(X), X)),
                                                                 transposeMatrix(X));
    const std::vector<double> &a = matrixVectorMultiply(M, y);

    a0 = a[0];
    a1 = a[1];
    a2 = a[2];
    a3 = a[3];

}