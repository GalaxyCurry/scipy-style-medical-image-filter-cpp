#include <iostream>
#include <vector>
#include <string>
#include <algorithm>
#include <filesystem>
#include <stdexcept>
#include <cstdint>
#include <cmath>
#include <cstddef>

#include "ImageFilter.h"

#include <dcmtk/dcmdata/dcfilefo.h>    // DcmFileFormat
#include <dcmtk/dcmdata/dcdatset.h>   // DcmDataset
#include <dcmtk/dcmdata/dcdeftag.h>    // 预定义DICOM标签
#include <dcmtk/ofstd/ofcond.h>      // OFCondition
#include <dcmtk/dcmdata/dctypes.h>     // 数据类型（Uint16等）
#include <dcmtk/dcmdata/dcxfer.h>      // 传输语法
#include <dcmtk/dcmdata/dcpixseq.h>    // 像素序列
#include <dcmtk/dcmimgle/dcmimage.h>   // DicomImage
#include <dcmtk/dcmdata/dcostrmz.h>    // 压缩保存
#include <dcmtk/ofstd/ofstring.h>      // OFString



// -------------------------- 工具函数 --------------------------
// 1. 遍历文件夹，获取所有.dcm文件路径
std::vector<std::string> get_all_dcm_files(const std::string& folder_path) {
    std::vector<std::string> dcm_paths;
    namespace fs = std::filesystem;

    if (!fs::exists(folder_path) || !fs::is_directory(folder_path)) {
        throw std::runtime_error("文件夹不存在或不是目录：" + folder_path);
    }

    for (const auto& entry : fs::directory_iterator(folder_path)) {
        if (entry.is_regular_file() && entry.path().extension() == ".dcm") {
            dcm_paths.push_back(entry.path().string());
        }
    }

    if (dcm_paths.empty()) {
        throw std::runtime_error("文件夹中未找到DCM文件：" + folder_path);
    }

    return dcm_paths;
}

// 2. 读取单个DCM文件的像素数据和z轴位置（用于排序）
struct DcmSlice {
    Mat2D pixel_data;  // 2D切片数据（double类型）
    double z_pos;      // z轴位置（用于按深度排序）
    int height;        // 切片高度
    int width;         // 切片宽度
    DcmFileFormat file_format; // 保存原始文件格式（用于复制元数据）
};

DcmSlice read_dcm_slice(const std::string& dcm_path) {
    DcmSlice slice;

    // 读取DCM文件
    OFCondition status = slice.file_format.loadFile(dcm_path.c_str());
    if (!status.good()) {
        throw std::runtime_error("读取DCM文件失败：" + dcm_path + "，错误：" + status.text());
    }

    // 关键修复1：const DcmFileFormat转换为非const，调用getDataset()
    DcmDataset* dataset = const_cast<DcmFileFormat*>(&slice.file_format)->getDataset();
    if (!dataset) {
        throw std::runtime_error("获取DCM数据集失败：" + dcm_path);
    }

    // 获取图像尺寸（Rows/Columns）- 非const数据集调用接口
    Uint16 rows = 0, cols = 0;
    // 去掉多余参数，使用默认值，避免匹配错误重载
    if (!dataset->findAndGetUint16(DCM_Rows, rows).good() ||
        !dataset->findAndGetUint16(DCM_Columns, cols).good()) {
        throw std::runtime_error("获取DCM图像尺寸失败：" + dcm_path);
    }
    slice.height = static_cast<int>(rows);
    slice.width = static_cast<int>(cols);

    // 获取z轴位置（SliceLocation）- 非const调用
    if (!dataset->findAndGetFloat64(DCM_SliceLocation, slice.z_pos).good()) {
        // 兼容无SliceLocation的情况，读取ImagePositionPatient的z分量
        const Float64* pos_ptr = nullptr; // 非const指针（适配接口）
        unsigned long pos_count = 0;
        if (dataset->findAndGetFloat64Array(DCM_ImagePositionPatient, pos_ptr, &pos_count).good() && pos_count >= 3) {
            slice.z_pos = pos_ptr[2];
            std::cout << "警告：DCM文件" << dcm_path << "无SliceLocation，使用ImagePositionPatient的z分量：" << slice.z_pos << std::endl;
        } else {
            throw std::runtime_error("获取DCM切片位置失败：" + dcm_path + "（无SliceLocation和ImagePositionPatient标签）");
        }
    }

    // 读取像素数据：直接使用非const dataset构造DicomImage
    DicomImage dcm_image(dataset, dataset->getOriginalXfer());
    if (dcm_image.getStatus() != EIS_Normal) {
        throw std::runtime_error("解析DCM像素数据失败：" + dcm_path + "，错误状态：" + std::to_string(dcm_image.getStatus()));
    }

    // 获取像素数据指针（16位无符号整数）
    const void* pixel_ptr = dcm_image.getOutputData(16);
    if (!pixel_ptr) {
        throw std::runtime_error("获取DCM像素指针失败：" + dcm_path);
    }
    const uint16_t* pixel_data_16 = static_cast<const uint16_t*>(pixel_ptr);

    // 填充到Mat2D（double类型）
    const size_t height = static_cast<size_t>(slice.height);
    const size_t width = static_cast<size_t>(slice.width);
    slice.pixel_data.resize(height, std::vector<double>(width));
    for (size_t i = 0; i < height; ++i) {
        for (size_t j = 0; j < width; ++j) {
            slice.pixel_data[i][j] = static_cast<double>(pixel_data_16[i * width + j]);
        }
    }

    return slice;
}

// 3. 将DcmSlice数组转换为Mat3D（[depth, height, width]）
Mat3D slices_to_mat3d(const std::vector<DcmSlice>& slices) {
    if (slices.empty()) {
        throw std::runtime_error("切片数组为空，无法转换为3D数组");
    }

    const size_t depth = slices.size();
    const size_t height = static_cast<size_t>(slices[0].height);
    const size_t width = static_cast<size_t>(slices[0].width);

    Mat3D mat3d(depth, Mat2D(height, std::vector<double>(width)));

    for (size_t z = 0; z < depth; ++z) {
        if (static_cast<size_t>(slices[z].height) != height || static_cast<size_t>(slices[z].width) != width) {
            throw std::runtime_error("第" + std::to_string(z) + "个切片尺寸不一致（预期：" + 
                                   std::to_string(height) + "x" + std::to_string(width) + 
                                   "，实际：" + std::to_string(slices[z].height) + "x" + std::to_string(slices[z].width) + "）");
        }
        for (size_t y = 0; y < height; ++y) {
            std::copy(slices[z].pixel_data[y].begin(), slices[z].pixel_data[y].end(), mat3d[z][y].begin());
        }
    }

    return mat3d;
}

// 4. 将Mat3D保存为一系列DCM文件（保留原始DICOM元数据）
void save_mat3d_to_dcm(const Mat3D& mat3d, const std::string& output_folder,
                      const std::vector<DcmSlice>& original_slices,
                      const std::vector<double>& spacing) {
    if (mat3d.empty() || mat3d[0].empty() || mat3d[0][0].empty()) {
        throw std::runtime_error("3D数组为空，无法保存DCM文件");
    }

    if (mat3d.size() != original_slices.size()) {
        throw std::runtime_error("处理后的数据与原始切片数量不匹配（处理后：" + 
                               std::to_string(mat3d.size()) + "，原始：" + std::to_string(original_slices.size()) + "）");
    }

    namespace fs = std::filesystem;
    if (!fs::exists(output_folder)) {
        if (!fs::create_directories(output_folder)) {
            throw std::runtime_error("创建输出文件夹失败：" + output_folder);
        }
    }

    const size_t depth = mat3d.size();
    const size_t height = mat3d[0].size();
    const size_t width = mat3d[0][0].size();

    for (size_t z = 0; z < depth; ++z) {
        if (mat3d[z].size() != height || mat3d[z][0].size() != width) {
            throw std::runtime_error("第" + std::to_string(z) + "个处理后切片尺寸异常");
        }

        // 复制原始元数据（深拷贝）
        DcmFileFormat file_format;
        OFCondition status = file_format.copyFrom(original_slices[z].file_format);
        if (!status.good()) {
            throw std::runtime_error("复制DCM元数据失败（切片" + std::to_string(z) + "）：" + status.text());
        }
        DcmDataset* dataset = file_format.getDataset(); // 保存时用非const（需修改）
        if (!dataset) {
            throw std::runtime_error("获取输出DCM数据集失败（切片" + std::to_string(z) + "）");
        }

        // 准备像素数据（16位无符号整数）
        std::vector<uint16_t> pixel_data_16(height * width);
        for (size_t y = 0; y < height; ++y) {
            for (size_t x = 0; x < width; ++x) {
                double val = std::clamp(mat3d[z][y][x], 0.0, 65535.0);
                pixel_data_16[y * width + x] = static_cast<uint16_t>(std::round(val));
            }
        }

        // 更新像素数据
        status = dataset->putAndInsertUint16Array(DCM_PixelData, pixel_data_16.data(), static_cast<unsigned long>(pixel_data_16.size()));
        if (!status.good()) {
            throw std::runtime_error("更新DCM像素数据失败（切片" + std::to_string(z) + "）：" + status.text());
        }

        // 更新像素间距
        if (spacing.size() >= 2) {
            dataset->putAndInsertFloat64(DCM_PixelSpacing, spacing[0], 0);
            dataset->putAndInsertFloat64(DCM_PixelSpacing, spacing[1], 1);
        }
        if (spacing.size() >= 3) {
            dataset->putAndInsertFloat64(DCM_SliceThickness, spacing[2]);
        }

        // 更新切片位置
        dataset->putAndInsertFloat64(DCM_SliceLocation, original_slices[z].z_pos);

        // 关键修复2：const DcmFileFormat转换为非const，获取传输语法
        DcmDataset* original_dataset = const_cast<DcmFileFormat*>(&original_slices[z].file_format)->getDataset();
        E_TransferSyntax xfer = original_dataset->getOriginalXfer();
        
        // 保存文件
        std::string output_path = output_folder + "/filtered_slice_" + std::to_string(z + 1) + ".dcm";
        status = file_format.saveFile(output_path.c_str(), xfer);
        if (!status.good()) {
            throw std::runtime_error("保存DCM文件失败：" + output_path + "，错误：" + status.text());
        }

        if ((z + 1) % 50 == 0) {
            std::cout << "已保存 " << (z + 1) << "/" << depth << " 个DCM文件" << std::endl;
        }
    }

    std::cout << "所有DCM文件保存完成！共 " << depth << " 个文件，输出路径：" << output_folder << std::endl;
}

// 5. 从DCM获取像素间距（x/y/z轴）
std::vector<double> get_dcm_spacing(const std::vector<DcmSlice>& slices) {
    std::vector<double> spacing = {1.0, 1.0, 1.0};
    if (slices.empty()) {
        return spacing;
    }

    // 关键修复3：const DcmFileFormat转换为非const，获取数据集
    DcmDataset* dataset = const_cast<DcmFileFormat*>(&slices[0].file_format)->getDataset();
    if (dataset) {
        double spacing_x = 1.0, spacing_y = 1.0;
        // 非const调用findAndGetFloat64
        if (dataset->findAndGetFloat64(DCM_PixelSpacing, spacing_x, 0).good() &&
            dataset->findAndGetFloat64(DCM_PixelSpacing, spacing_y, 1).good()) {
            spacing[0] = spacing_x;
            spacing[1] = spacing_y;
        } else {
            std::cout << "警告：无法获取像素间距，使用默认值（1.0, 1.0）" << std::endl;
        }

        // 读取切片厚度
        double slice_thickness = 1.0;
        if (dataset->findAndGetFloat64(DCM_SliceThickness, slice_thickness).good()) {
            spacing[2] = slice_thickness;
        }
    }

    // 计算z轴间距（如果未获取到切片厚度）
    if (slices.size() >= 2 && spacing[2] <= 0.01) {
        double z_diff = std::fabs(slices[1].z_pos - slices[0].z_pos);
        if (z_diff > 0.01) {
            spacing[2] = z_diff;
        }
    }

    return spacing;
}




// -------------------------- 主函数 --------------------------
int main(int argc, char* argv[]) {
    // 标记未使用参数，避免警告
    (void)argc;
    (void)argv;

    try {
        // 配置参数（修改为实际路径）
        const std::string dcm_folder = "D:/tasks/Smart/Smart_Screw_Inspection/Data/testdata_dcm/haidian_S009_REFVOL00012";
        const std::string output_folder = "D:/tasks/Smart/Smart_Screw_Inspection/Data/testdata_dcm/haidian_s009_12";
        const double gaussian_sigma = 4;
        const int border_type = 1;
        const bool use_gaussian_filter = true;
        const bool use_sobel_filter = false;

        // 读取DCM文件路径
        std::cout << "===== 开始读取DCM文件 =====" << std::endl;
        std::vector<std::string> dcm_paths = get_all_dcm_files(dcm_folder);
        std::cout << "成功找到 " << dcm_paths.size() << " 个DCM文件" << std::endl;

        // 解析DCM切片
        std::cout << "\n===== 开始解析DCM切片 =====" << std::endl;
        std::vector<DcmSlice> slices;
        slices.reserve(dcm_paths.size());
        for (size_t i = 0; i < dcm_paths.size(); ++i) {
            std::cout << "解析第 " << (i + 1) << "/" << dcm_paths.size() << " 个文件：" << dcm_paths[i] << std::endl;
            slices.push_back(read_dcm_slice(dcm_paths[i]));
        }

        // 按z轴排序
        std::cout << "\n===== 按切片位置排序 =====" << std::endl;
        std::sort(slices.begin(), slices.end(), [](const DcmSlice& a, const DcmSlice& b) {
            return a.z_pos < b.z_pos;
        });
        std::cout << "排序完成，共 " << slices.size() << " 个切片" << std::endl;

        // 构建3D体数据
        std::cout << "\n===== 构建3D体数据 =====" << std::endl;
        Mat3D input_vol = slices_to_mat3d(slices);
        const size_t depth = input_vol.size();
        const size_t height = input_vol[0].size();
        const size_t width = input_vol[0][0].size();
        std::cout << "3D体数据尺寸：z=" << depth << " × y=" << height << " × x=" << width << std::endl;

        // 获取像素间距
        std::vector<double> spacing = get_dcm_spacing(slices);
        std::cout << "像素间距：x=" << spacing[0] << "mm, y=" << spacing[1] << "mm, z=" << spacing[2] << "mm" << std::endl;

        // 滤波处理
        std::cout << "\n===== 开始滤波处理 =====" << std::endl;
        Mat3D filtered_vol;
        if (use_gaussian_filter) {
            std::cout << "执行高斯滤波（sigma=" << gaussian_sigma << "）..." << std::endl;
            filtered_vol = ImageFilter::gaussian_filter(input_vol, gaussian_sigma, border_type);
        } else if (use_sobel_filter) {
            const int sobel_axis = 2;
            std::cout << "执行Sobel滤波（轴=" << sobel_axis << "）..." << std::endl;
            filtered_vol = ImageFilter::sobel(input_vol, sobel_axis, border_type);
        } else {
            throw std::runtime_error("未选择任何滤波方式！请设置use_gaussian_filter或use_sobel_filter为true");
        }
        std::cout << "滤波处理完成" << std::endl;

        // 保存DCM文件
        std::cout << "\n===== 开始保存DCM文件 =====" << std::endl;
        save_mat3d_to_dcm(filtered_vol, output_folder, slices, spacing);

        std::cout << "\n===== 所有流程执行完成！=====" << std::endl;

    } catch (const std::exception& e) {
        std::cerr << "\n===== 错误：" << e.what() << " =====" << std::endl;
        return 1;
    }

    return 0;
}