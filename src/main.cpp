#include <iostream>
#include <vector>
#include <string>
#include <algorithm>
#include <filesystem>
#include <stdexcept>
#include <cstdint>
#include <cmath>
#include <cstddef>

#include "ImageFilter.h"  // 你的滤波类头文件

// ITK库包含（严格适配ITK 5.4）
#include "itkImage.h"
#include "itkGDCMImageIO.h"
#include "itkGDCMSeriesFileNames.h"
#include "itkImageSeriesReader.h"
#include "itkImageSeriesWriter.h"
#include "itkRescaleIntensityImageFilter.h"
#include "itkMetaDataDictionary.h"
#include "itkMetaDataObject.h"

// 定义ITK图像类型（保持不变）
using PixelType = uint16_t;
constexpr unsigned int Dimension = 3;
using ImageType = itk::Image<PixelType, Dimension>;
using ReaderType = itk::ImageSeriesReader<ImageType>;
using ImageIOType = itk::GDCMImageIO;
using SeriesFileNamesType = itk::GDCMSeriesFileNames;
using WriterType = itk::ImageSeriesWriter<ImageType, ImageType>;

// 定义矩阵类型（保持不变）
using Mat2D = std::vector<std::vector<double>>;
using Mat3D = std::vector<Mat2D>;

// -------------------------- 工具函数（完全适配ITK 5.4） --------------------------
// 1. 获取DICOM序列文件路径（ITK 5.4正确用法）
std::vector<std::string> get_all_dcm_files(const std::string& folder_path) {
    namespace fs = std::filesystem;
    if (!fs::exists(folder_path) || !fs::is_directory(folder_path)) {
        throw std::runtime_error("文件夹不存在或不是目录：" + folder_path);
    }

    auto seriesFileNames = SeriesFileNamesType::New();
    seriesFileNames->SetDirectory(folder_path);  // 必须先设置目录
    
    // ITK 5.4：GetSeriesUIDs() 获取所有序列UID
    std::vector<std::string> seriesUIDs = seriesFileNames->GetSeriesUIDs();
    if (seriesUIDs.empty()) {
        throw std::runtime_error("文件夹中未找到DCM序列：" + folder_path);
    }

    // ITK 5.4：GetFileNames() 必须传入序列UID作为参数（无SetSeriesUID方法）
    const std::string targetSeriesUID = seriesUIDs[0];
    std::vector<std::string> dcm_paths = seriesFileNames->GetFileNames(targetSeriesUID);
    
    if (dcm_paths.empty()) {
        throw std::runtime_error("获取DCM文件路径失败（UID：" + targetSeriesUID + "）");
    }

    return dcm_paths;
}

// 2. 读取DICOM序列并转换为Mat3D（ITK 5.4正确用法）
void read_dcm_series(const std::string& folder_path, Mat3D& volume, 
                    std::vector<double>& spacing, 
                    // 元数据类型：适配ITK 5.4的返回值（const std::vector<MetaDataDictionary*>*）
                    std::vector<itk::MetaDataDictionary*>& metaDictionaries) {
    // 创建读取器和IO对象
    auto reader = ReaderType::New();
    auto dicomIO = ImageIOType::New();
    reader->SetImageIO(dicomIO);

    // 获取序列文件名称（ITK 5.4正确用法）
    auto seriesFileNames = SeriesFileNamesType::New();
    seriesFileNames->SetDirectory(folder_path);
    
    std::vector<std::string> seriesUIDs = seriesFileNames->GetSeriesUIDs();
    if (seriesUIDs.empty()) {
        throw std::runtime_error("未找到DICOM序列：" + folder_path);
    }

    // 直接传入UID获取文件列表（无SetSeriesUID）
    std::vector<std::string> fileNames = seriesFileNames->GetFileNames(seriesUIDs[0]);
    reader->SetFileNames(fileNames);

    try {
        reader->Update();
    } catch (const itk::ExceptionObject& e) {
        throw std::runtime_error("读取DICOM序列失败：" + std::string(e.GetDescription()));
    }

    // 获取图像信息（保持不变）
    ImageType::Pointer image = reader->GetOutput();
    ImageType::RegionType region = image->GetLargestPossibleRegion();
    ImageType::SizeType size = region.GetSize();
    
    // 获取间距信息（保持不变）
    ImageType::SpacingType itkSpacing = image->GetSpacing();
    spacing = {itkSpacing[0], itkSpacing[1], itkSpacing[2]};

    // ITK 5.4：处理元数据字典（返回值为 const std::vector<MetaDataDictionary*>*）
    auto metaDictArray = reader->GetMetaDataDictionaryArray();
    if (!metaDictArray || metaDictArray->empty()) {
        throw std::runtime_error("获取DICOM元数据失败");
    }
    // 复制到非const向量（适配写入器接口）
    metaDictionaries.clear();
    for (auto dictPtr : *metaDictArray) {
        metaDictionaries.push_back(const_cast<itk::MetaDataDictionary*>(dictPtr));
    }

    // 转换为Mat3D格式（保持不变）
    volume.resize(size[2]);
    ImageType::IndexType index;

    for (size_t z = 0; z < size[2]; ++z) {
        volume[z].resize(size[1], std::vector<double>(size[0]));
        index[2] = z;
        
        for (size_t y = 0; y < size[1]; ++y) {
            index[1] = y;
            for (size_t x = 0; x < size[0]; ++x) {
                index[0] = x;
                volume[z][y][x] = static_cast<double>(image->GetPixel(index));
            }
        }
    }

    std::cout << "3D体数据尺寸：z=" << size[2] << " × y=" << size[1] << " × x=" << size[0] << std::endl;
}

// 3. 保存Mat3D为DICOM序列（ITK 5.4正确用法）
void save_mat3d_to_dcm(const Mat3D& mat3d, const std::string& output_folder,
                      // 元数据类型：ITK 5.4写入器需要 std::vector<MetaDataDictionary*>*
                      const std::vector<itk::MetaDataDictionary*>& originalMetaDictionaries,
                      const std::vector<double>& spacing) {
    if (mat3d.empty() || mat3d[0].empty() || mat3d[0][0].empty()) {
        throw std::runtime_error("3D数组为空，无法保存DCM文件");
    }

    if (mat3d.size() != originalMetaDictionaries.size()) {
        throw std::runtime_error("处理后的数据与原始切片数量不匹配（处理后：" + 
                               std::to_string(mat3d.size()) + "，原始：" + std::to_string(originalMetaDictionaries.size()) + "）");
    }

    namespace fs = std::filesystem;
    if (!fs::exists(output_folder)) {
        if (!fs::create_directories(output_folder)) {
            throw std::runtime_error("创建输出文件夹失败：" + output_folder);
        }
    }

    // 创建图像并填充数据（保持不变）
    const size_t depth = mat3d.size();
    const size_t height = mat3d[0].size();
    const size_t width = mat3d[0][0].size();

    ImageType::Pointer image = ImageType::New();
    ImageType::RegionType region;
    ImageType::IndexType start;
    ImageType::SizeType size;

    start[0] = 0;
    start[1] = 0;
    start[2] = 0;

    size[0] = width;
    size[1] = height;
    size[2] = depth;

    region.SetSize(size);
    region.SetIndex(start);
    image->SetRegions(region);

    // 设置间距（保持不变）
    ImageType::SpacingType itkSpacing;
    itkSpacing[0] = spacing[0];
    itkSpacing[1] = spacing[1];
    itkSpacing[2] = spacing[2];
    image->SetSpacing(itkSpacing);

    // 分配内存并填充数据（保持不变）
    image->Allocate();
    ImageType::IndexType index;

    for (size_t z = 0; z < depth; ++z) {
        if (mat3d[z].size() != height || mat3d[z][0].size() != width) {
            throw std::runtime_error("第" + std::to_string(z) + "个处理后切片尺寸异常");
        }

        index[2] = z;
        for (size_t y = 0; y < height; ++y) {
            index[1] = y;
            for (size_t x = 0; x < width; ++x) {
                index[0] = x;
                double val = std::clamp(mat3d[z][y][x], 0.0, 65535.0);
                image->SetPixel(index, static_cast<PixelType>(std::round(val)));
            }
        }
    }

    // 设置写入器（保持不变）
    auto writer = WriterType::New();
    auto dicomIO = ImageIOType::New();
    writer->SetImageIO(dicomIO);
    writer->SetInput(image);

    // 生成输出文件名（保持不变）
    std::vector<std::string> outputFileNames;
    for (size_t i = 0; i < depth; ++i) {
        std::string fileName = output_folder + "/filtered_slice_" + std::to_string(i + 1) + ".dcm";
        outputFileNames.push_back(fileName);
    }
    writer->SetFileNames(outputFileNames);

    // ITK 5.4：SetMetaDataDictionaryArray 需要传入 std::vector<MetaDataDictionary*>*
    writer->SetMetaDataDictionaryArray(&originalMetaDictionaries);

    // 写入DICOM文件（保持不变）
    try {
        writer->Update();
    } catch (const itk::ExceptionObject& e) {
        throw std::runtime_error("保存DICOM文件失败：" + std::string(e.GetDescription()));
    }

    std::cout << "所有DCM文件保存完成！共 " << depth << " 个文件，输出路径：" << output_folder << std::endl;
}

// -------------------------- 主函数（适配ITK 5.4） --------------------------
int main(int argc, char* argv[]) {
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

        // 读取DCM文件路径（保持不变）
        std::cout << "===== 开始读取DCM文件 =====" << std::endl;
        std::vector<std::string> dcm_paths = get_all_dcm_files(dcm_folder);
        std::cout << "成功找到 " << dcm_paths.size() << " 个DCM文件" << std::endl;

        // 解析DCM序列并构建3D体数据（元数据类型改为非const指针向量）
        std::cout << "\n===== 开始解析DCM序列 =====" << std::endl;
        Mat3D input_vol;
        std::vector<double> spacing;
        std::vector<itk::MetaDataDictionary*> metaDictionaries;  // ITK 5.4适配类型
        
        read_dcm_series(dcm_folder, input_vol, spacing, metaDictionaries);
        
        const size_t depth = input_vol.size();
        const size_t height = input_vol[0].size();
        const size_t width = input_vol[0][0].size();
        std::cout << "3D体数据尺寸：z=" << depth << " × y=" << height << " × x=" << width << std::endl;
        std::cout << "像素间距：x=" << spacing[0] << "mm, y=" << spacing[1] << "mm, z=" << spacing[2] << "mm" << std::endl;

        // 滤波处理（保持不变）
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

        // 保存DCM文件（保持不变）
        std::cout << "\n===== 开始保存DCM文件 =====" << std::endl;
        save_mat3d_to_dcm(filtered_vol, output_folder, metaDictionaries, spacing);

        std::cout << "\n===== 所有流程执行完成！=====" << std::endl;

    } catch (const std::exception& e) {
        std::cerr << "\n===== 错误：" << e.what() << " =====" << std::endl;
        return 1;
    }

    return 0;
}