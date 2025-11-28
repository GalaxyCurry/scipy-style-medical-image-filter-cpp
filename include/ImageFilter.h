#ifndef IMAGE_FILTER_H
#define IMAGE_FILTER_H

#include <vector>
#include <cmath>

using Mat3D = std::vector<std::vector<std::vector<double>>>;
using Mat2D = std::vector<std::vector<double>>;

class ImageFilter {
public:
    /**
     * @brief 判断两个浮点数是否接近（在误差范围内）
     * @param a 第一个浮点数
     * @param b 第二个浮点数
     * @param eps 误差范围，默认1e-6
     * @return 若接近则返回true，否则返回false
     */
    static bool isClose(double a, double b, double eps = 1e-6);

    /**
     * @brief 生成1D高斯核
     * @param sigma 高斯标准差
     * @param radius 核半径（核大小为2*radius+1）
     * @return 归一化后的1D高斯核
     */
    static std::vector<double> gaussian_kernel1d(double sigma, int radius);

    /**
     * @brief 提取3D矩阵中的第z个2D切片
     * @param mat3d 输入3D矩阵
     * @param z 切片索引
     * @return 提取的2D切片
     * @throws std::invalid_argument 若输入矩阵为空
     * @throws std::out_of_range 若z索引越界
     */
    static Mat2D extractSlice(const Mat3D& mat3d, int z);

    /**
     * @brief 对2D矩阵进行边界填充
     * @param input 输入2D矩阵
     * @param pad_row 行方向填充量
     * @param pad_col 列方向填充量
     * @param borderType 边界填充类型（0:CONSTANT, 1:REPLICATE, 2:REFLECT, 3:REFLECT_101）
     * @param cval 当borderType为CONSTANT时的填充值，默认0.0
     * @return 填充后的2D矩阵
     */
    static Mat2D pad2D(const Mat2D& input, int pad_row, int pad_col,
                      int borderType, double cval = 0.0);

    /**
     * @brief 对3D矩阵进行边界填充
     * @param input 输入3D矩阵
     * @param pads 填充量数组（[pad_row, pad_col, pad_depth]）
     * @param borderType 边界填充类型（0:CONSTANT, 1:REPLICATE, 2:REFLECT, 3:REFLECT_101）
     * @param cval 当borderType为CONSTANT时的填充值，默认0.0
     * @return 填充后的3D矩阵
     * @throws std::invalid_argument 若输入矩阵为空或pads元素不足
     */
    static Mat3D pad3D(const Mat3D& input, const std::vector<int>& pads,
                      int borderType, double cval = 0.0);

    /**
     * @brief 对3D矩阵进行1D相关运算
     * @param input 输入3D矩阵
     * @param weights 1D卷积核
     * @param axis 运算轴（0:行, 1:列, 2:深度）
     * @param output 输出结果矩阵
     * @param borderType 边界填充类型，默认1(REPLICATE)
     * @param cval 当borderType为CONSTANT时的填充值，默认0.0
     * @throws std::invalid_argument 若输入矩阵为空、核为空或轴无效
     */
    static void correlate1d(const Mat3D& input, const std::vector<double>& weights,
                           int axis, Mat3D& output, int borderType = 1, double cval = 0.0);

    /**
     * @brief 对3D矩阵进行1D高斯滤波
     * @param input 输入3D矩阵
     * @param sigma 高斯标准差
     * @param axis 滤波轴（0:行, 1:列, 2:深度）
     * @param output 输出结果矩阵
     * @param borderType 边界填充类型，默认1(REPLICATE)
     * @param cval 当borderType为CONSTANT时的填充值，默认0.0
     */
    static void gaussian_filter1d(const Mat3D& input, double sigma, int axis,
                                 Mat3D& output, int borderType = 1, double cval = 0.0);

    /**
     * @brief 对3D矩阵进行3D高斯滤波（在所有轴上依次应用1D高斯滤波）
     * @param input 输入3D矩阵
     * @param sigma 高斯标准差
     * @param borderType 边界填充类型，默认1(REPLICATE)
     * @param cval 当borderType为CONSTANT时的填充值，默认0.0
     * @return 滤波后的3D矩阵
     */
    static Mat3D gaussian_filter(const Mat3D& input, double sigma,
                                int borderType = 1, double cval = 0.0);

    /**
     * @brief 对3D矩阵进行Sobel滤波（计算指定轴方向的梯度）
     * @param input 输入3D矩阵
     * @param axis 梯度方向轴（0:行, 1:列, 2:深度），默认0
     * @param borderType 边界填充类型，默认1(REPLICATE)
     * @param cval 当borderType为CONSTANT时的填充值，默认0.0
     * @return 梯度计算结果
     * @throws std::invalid_argument 若输入矩阵为空或轴无效
     */
    static Mat3D sobel(const Mat3D& input, int axis = 0,
                      int borderType = 1, double cval = 0.0);

private:
    /**
     * @brief 计算边界填充的镜像索引
     * @param idx 原始索引
     * @param size 维度大小
     * @param borderType 边界填充类型
     * @return 处理后的索引
     */
    static int getMirrorIndex(int idx, int size, int borderType);

    /**
     * @brief 获取边界填充值
     * @param input 输入2D矩阵
     * @param row 行索引
     * @param col 列索引
     * @param rows 矩阵行数
     * @param cols 矩阵列数
     * @param borderType 边界填充类型
     * @param cval 当borderType为CONSTANT时的填充值
     * @return 填充值
     */
    static double getBorderValue(const Mat2D& input, int row, int col,
                               int rows, int cols, int borderType, double cval);
};

#endif 