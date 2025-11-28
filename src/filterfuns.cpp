#include <iostream>
#include <vector>
#include <algorithm>
#include <cmath>
#include <stdexcept>

#include "ImageFilter.h"


bool ImageFilter::isClose(double a, double b, double eps) {
    return std::fabs(a - b) < eps;
}

// 生成1D高斯核
std::vector<double> ImageFilter::gaussian_kernel1d(double sigma, int radius) {

    double sigma2 = sigma * sigma;
    std::vector<double> kernel(2 * radius + 1);
    double sum = 0.0;

    for (size_t i = 0; i < kernel.size(); ++i) {
        int x = static_cast<int>(i) - radius;
        kernel[i] = std::exp(-0.5 / sigma2 * x * x);
        sum += kernel[i];
    }

    // 归一化
    for (auto& val : kernel) {
        val /= sum;
    }

    return kernel;
}


Mat2D ImageFilter::extractSlice(const Mat3D& mat3d, int z) {
    if (mat3d.empty()) throw std::invalid_argument("Input 3D matrix is empty");
    if (z < 0 || z >= static_cast<int>(mat3d.size())) {
        throw std::out_of_range("Z index out of range");
    }
    return mat3d[z];
}

// 2D边界填充
Mat2D ImageFilter::pad2D(const Mat2D& input, int pad_row, int pad_col,
                    int borderType, double cval) {
    if (input.empty()) return {};
    int rows = input.size();
    int cols = input[0].size();

    int new_rows = rows + 2 * pad_row;
    int new_cols = cols + 2 * pad_col;
    Mat2D result(new_rows, std::vector<double>(new_cols, 0.0));

    // 填充中间原始数据
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            result[pad_row + i][pad_col + j] = input[i][j];
        }
    }

    // 填充上下边界
    for (int i = 0; i < pad_row; ++i) {
        int row_idx = i;                  // 上边界行索引
        int mirror_row = pad_row - i - 1; // 镜像行索引

        for (int j = 0; j < new_cols; ++j) {
            // 计算原始列索引
            int col_idx = j - pad_col;
            int original_col = getMirrorIndex(col_idx, cols, borderType);

            double val;
            if (col_idx >= 0 && col_idx < cols) {
                val = getBorderValue(input, mirror_row - pad_row, original_col, 
                                    rows, cols, borderType, cval);
            } else {
                val = getBorderValue(input, mirror_row - pad_row, original_col, 
                                    rows, cols, borderType, cval);
            }
            result[row_idx][j] = val;
            result[new_rows - 1 - i][j] = val; // 下边界对称位置
        }
    }

    // 填充左右边界（不包括已填充的四角）
    for (int i = pad_row; i < new_rows - pad_row; ++i) {
        for (int j = 0; j < pad_col; ++j) {
            int col_idx = j;                  // 左边界列索引
            int mirror_col = pad_col - j - 1; // 镜像列索引
            int original_row = i - pad_row;

            double val = getBorderValue(input, original_row, mirror_col - pad_col, 
                                        rows, cols, borderType, cval);
            result[i][col_idx] = val;
            result[i][new_cols - 1 - j] = val; // 右边界对称位置
        }
    }

    return result;
}

// 3D边界填充
Mat3D ImageFilter::pad3D(const Mat3D& input, const std::vector<int>& pads,
                    int borderType, double cval) {
    if (input.empty()) return {};
    if (pads.size() < 2) throw std::invalid_argument("Pads must have at least 2 elements");
    
    int pad_row = pads[0];
    int pad_col = pads[1];
    int pad_depth = pads.size() > 2 ? pads[2] : 0;
    int depth = input.size();

    // 先处理每个2D切片的行和列填充
    std::vector<Mat2D> sliced_padded(depth);
    for (int z = 0; z < depth; ++z) {
        sliced_padded[z] = pad2D(input[z], pad_row, pad_col, borderType, cval);
    }

    // 处理深度方向填充
    int new_depth = depth + 2 * pad_depth;
    int new_rows = sliced_padded[0].size();
    int new_cols = sliced_padded[0][0].size();
    Mat3D result(new_depth, Mat2D(new_rows, std::vector<double>(new_cols, 0.0)));

    // 复制中间原始深度切片
    for (int z = 0; z < depth; ++z) {
        result[pad_depth + z] = sliced_padded[z];
    }

    // 填充深度方向前半部分
    for (int z = 0; z < pad_depth; ++z) {
        int src_z = getMirrorIndex(z - pad_depth, depth, borderType);
        src_z = std::clamp(src_z, 0, depth - 1);
        
        for (int i = 0; i < new_rows; ++i) {
            for (int j = 0; j < new_cols; ++j) {
                if (borderType == 0) { // CONSTANT
                    result[z][i][j] = cval;
                } else {
                    result[z][i][j] = sliced_padded[src_z][i][j];
                }
            }
        }
    }

    // 填充深度方向后半部分
    for (int z = 0; z < pad_depth; ++z) {
        int dst_z = pad_depth + depth + z;
        int src_z = getMirrorIndex(depth + z, depth, borderType);
        src_z = std::clamp(src_z, 0, depth - 1);
        
        for (int i = 0; i < new_rows; ++i) {
            for (int j = 0; j < new_cols; ++j) {
                if (borderType == 0) { // CONSTANT
                    result[dst_z][i][j] = cval;
                } else {
                    result[dst_z][i][j] = sliced_padded[src_z][i][j];
                }
            }
        }
    }

    return result;
}




    // 1D相关运算
void ImageFilter::correlate1d(const Mat3D& input, const std::vector<double>& weights,
                        int axis, Mat3D& output, int borderType, double cval) {
    if (input.empty() || weights.empty()) return;
    if (axis < 0 || axis > 2) throw std::invalid_argument("Invalid axis (0-2)");

    int ksize = weights.size();
    int k_half = ksize / 2;
    std::vector<int> pads(3, 0);
    pads[axis] = k_half;

    // 边界填充
    Mat3D padded = pad3D(input, pads, borderType, cval);

    // 初始化输出
    int depth = input.size();
    int rows = input[0].size();
    int cols = input[0][0].size();
    output = Mat3D(depth, Mat2D(rows, std::vector<double>(cols, 0.0)));

    // 判断核对称性
    bool symmetric = true;
    bool anti_symmetric = true;
    for (int i = 1; i <= k_half; ++i) {
        if (!isClose(weights[k_half + i], weights[k_half - i])) symmetric = false;
        if (!isClose(weights[k_half + i], -weights[k_half - i])) anti_symmetric = false;
    }

    if (axis == 0) {  // 行方向
        for (int z = 0; z < depth; ++z) {
            for (int c = 0; c < cols; ++c) {
                for (int r = 0; r < rows; ++r) {
                    int center_y = k_half + r;
                    double sum = 0.0;

                    if (symmetric) {
                        sum = padded[z + pads[2]][center_y][c + pads[1]] * weights[k_half];
                        for (int i = 1; i <= k_half; ++i) {
                            sum += (padded[z + pads[2]][center_y - i][c + pads[1]] +
                                    padded[z + pads[2]][center_y + i][c + pads[1]]) * weights[k_half + i];
                        }
                    } else if (anti_symmetric) {
                        sum = padded[z + pads[2]][center_y][c + pads[1]] * weights[k_half];
                        for (int i = 1; i <= k_half; ++i) {
                            sum += (padded[z + pads[2]][center_y - i][c + pads[1]] -
                                    padded[z + pads[2]][center_y + i][c + pads[1]]) * weights[k_half + i];
                        }
                    } else {
                        for (int i = 0; i < ksize; ++i) {
                            sum += padded[z + pads[2]][center_y - k_half + i][c + pads[1]] * weights[i];
                        }
                    }

                    output[z][r][c] = sum;
                }
            }
        }
    } else if (axis == 1) {  // 列方向
        for (int z = 0; z < depth; ++z) {
            for (int r = 0; r < rows; ++r) {
                for (int c = 0; c < cols; ++c) {
                    int center_x = k_half + c;
                    double sum = 0.0;

                    if (symmetric) {
                        sum = padded[z + pads[2]][r + pads[0]][center_x] * weights[k_half];
                        for (int i = 1; i <= k_half; ++i) {
                            sum += (padded[z + pads[2]][r + pads[0]][center_x - i] +
                                    padded[z + pads[2]][r + pads[0]][center_x + i]) * weights[k_half + i];
                        }
                    } else if (anti_symmetric) {
                        sum = padded[z + pads[2]][r + pads[0]][center_x] * weights[k_half];
                        for (int i = 1; i <= k_half; ++i) {
                            sum += (padded[z + pads[2]][r + pads[0]][center_x - i] -
                                    padded[z + pads[2]][r + pads[0]][center_x + i]) * weights[k_half + i];
                        }
                    } else {
                        for (int i = 0; i < ksize; ++i) {
                            sum += padded[z + pads[2]][r + pads[0]][center_x - k_half + i] * weights[i];
                        }
                    }

                    output[z][r][c] = sum;
                }
            }
        }
    } else if (axis == 2) {  // 深度方向
        for (int r = 0; r < rows; ++r) {
            for (int c = 0; c < cols; ++c) {
                for (int z = 0; z < depth; ++z) {
                    int center_z = k_half + z;
                    double sum = 0.0;

                    if (symmetric) {
                        sum = padded[center_z][r + pads[0]][c + pads[1]] * weights[k_half];
                        for (int i = 1; i <= k_half; ++i) {
                            sum += (padded[center_z - i][r + pads[0]][c + pads[1]] +
                                    padded[center_z + i][r + pads[0]][c + pads[1]]) * weights[k_half + i];
                        }
                    } else if (anti_symmetric) {
                        sum = padded[center_z][r + pads[0]][c + pads[1]] * weights[k_half];
                        for (int i = 1; i <= k_half; ++i) {
                            sum += (padded[center_z - i][r + pads[0]][c + pads[1]] -
                                    padded[center_z + i][r + pads[0]][c + pads[1]]) * weights[k_half + i];
                        }
                    } else {
                        for (int i = 0; i < ksize; ++i) {
                            sum += padded[center_z - k_half + i][r + pads[0]][c + pads[1]] * weights[i];
                        }
                    }

                    output[z][r][c] = sum;
                }
            }
        }
    }
}











// 1D高斯滤波
void ImageFilter::gaussian_filter1d(const Mat3D& input, double sigma, int axis,
                                Mat3D& output, int borderType, double cval) {
    int radius = static_cast<int>(4 * sigma + 0.5);
    auto kernel = gaussian_kernel1d(sigma, radius);
    std::reverse(kernel.begin(), kernel.end()); // 卷积需要核反转
    correlate1d(input, kernel, axis, output, borderType, cval);
}

// 高斯滤波
Mat3D ImageFilter::gaussian_filter(const Mat3D& input, double sigma, 
                            int borderType, double cval) {
    if (input.empty()) return {};
    Mat3D result = input;
    int dims = 3; // 假设3D数据

    for (int axis = 0; axis < dims; ++axis) {
        Mat3D temp;
        gaussian_filter1d(result, sigma, axis, temp, borderType, cval);
        result = temp;
    }

    return result;
}










// Sobel滤波
Mat3D ImageFilter::sobel(const Mat3D& input, int axis, 
                    int borderType, double cval) {
    if (input.empty()) return {};
    if (axis < 0 || axis > 2) throw std::invalid_argument("Invalid axis (0-2)");

    // 梯度核
    std::vector<double> grad_kernel = {-1, 0, 1};
    Mat3D result;
    correlate1d(input, grad_kernel, axis, result, borderType, cval);

    // 平滑核
    std::vector<double> smooth_kernel = {1, 2, 1};
    for (int ax = 0; ax < 3; ++ax) {
        if (ax != axis) {
            Mat3D temp;
            correlate1d(result, smooth_kernel, ax, temp, borderType, cval);
            result = temp;
        }
    }

    return result;
}




// 获取边界填充的镜像索引
int ImageFilter::getMirrorIndex(int idx, int size, int borderType) {
    if (size <= 0) return 0;
    
    switch (borderType) {
        case 1:  // REPLICATE
            return std::clamp(idx, 0, size - 1);
        case 2:  // REFLECT
            if (idx < 0) return -idx - 1;
            if (idx >= size) return 2 * size - idx - 1;
            return idx;
        case 3:  // REFLECT_101
            if (idx < 0) return -idx;
            if (idx >= size) return 2 * size - idx - 2;
            return idx;
        default: // CONSTANT
            return 0;
    }
}

// 获取边界填充值
double ImageFilter::getBorderValue(const Mat2D& input, int row, int col,
                            int rows, int cols, int borderType, double cval) {
    if (row >= 0 && row < rows && col >= 0 && col < cols) {
        return input[row][col];
    }

    if (borderType == 0) { // CONSTANT
        return cval;
    }

    int r = getMirrorIndex(row, rows, borderType);
    int c = getMirrorIndex(col, cols, borderType);
    return input[r][c];
}




/* 
int main() {
  try {
    
    Mat3D input(1, Mat2D(3, std::vector<double>(3, 0.0)));
    input[0][0] = {1, 2, 3};
    input[0][1] = {4, 5, 6};
    input[0][2] = {7, 8, 9};

    Mat3D gaussian_out = ImageFilter::gaussian_filter(input, 4.0, 1);  // borderType=REPLICATE

    std::cout << "高斯滤波后的结果：" << std::endl;
    for (int i = 0; i < gaussian_out[0].size(); ++i) {
      for (int j = 0; j < gaussian_out[0][i].size(); ++j) {
        std::cout << gaussian_out[0][i][j] << " ";
      }
      std::cout << std::endl;
    }

    Mat3D sobel_out = ImageFilter::sobel(input, 1);  // axis=1（列方向）
    std::cout << "\nSobel 滤波后的结果：" << std::endl;
    for (int i = 0; i < sobel_out[0].size(); ++i) {
      for (int j = 0; j < sobel_out[0][i].size(); ++j) {
        std::cout << sobel_out[0][i][j] << " ";
      }
      std::cout << std::endl;
    }

  } catch (const std::exception& e) {
    std::cerr << "错误：" << e.what() << std::endl;
    return 1;
  }

  return 0;
} */