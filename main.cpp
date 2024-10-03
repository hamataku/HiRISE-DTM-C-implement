#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
#include <gdal.h>
#include <gdal_priv.h>
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include <stb_image_write.h>
#include <cstring>

#define ALTITUDE_MULTIPLIER 2
#define MESH_STEP 1
#define NORMAL_STEP 1

struct Vec3 {
    float x, y, z;
};

void compute_mesh(const std::vector<std::vector<float>>& data, int rows, int cols, int step, std::vector<Vec3>& positions)
{
    std::cout << "compute_mesh" << std::endl;
    float h = ALTITUDE_MULTIPLIER;
    for (int i0 = 0; i0 < rows - 1; i0 += step) {
        for (int j0 = 0; j0 < cols - 1; j0 += step) {
            int i1 = i0 + step;
            int j1 = j0 + step;
            float y00 = data[i0][j0];
            float y01 = data[i0][j1];
            float y10 = data[i1][j0];
            float y11 = data[i1][j1];
            int z0 = i0 * step, z1 = i1 * step;
            int x0 = j0 * step, x1 = j1 * step;

            if (!std::isnan(y00) && !std::isnan(y01) && !std::isnan(y10)) {
                positions.push_back({(float)x0, y10 * h, (float)z1});
                positions.push_back({(float)x1, y01 * h, (float)z0});
                positions.push_back({(float)x0, y00 * h, (float)z0});
            }
            if (!std::isnan(y11) && !std::isnan(y10) && !std::isnan(y01)) {
                positions.push_back({(float)x1, y11 * h, (float)z1});
                positions.push_back({(float)x1, y01 * h, (float)z0});
                positions.push_back({(float)x0, y10 * h, (float)z1});
            }
        }
    }
    std::cout << positions.size() / 3 << std::endl;
}

Vec3 normal_from_points(Vec3 a, Vec3 b, Vec3 c)
{
    Vec3 ab = {b.x - a.x, b.y - a.y, b.z - a.z};
    Vec3 ac = {c.x - a.x, c.y - a.y, c.z - a.z};
    Vec3 normal;
    normal.x = ab.y * ac.z - ab.z * ac.y;
    normal.y = ab.z * ac.x - ab.x * ac.z;
    normal.z = ab.x * ac.y - ab.y * ac.x;
    float d = std::sqrt(normal.x * normal.x + normal.y * normal.y + normal.z * normal.z);
    return {normal.x / d, normal.y / d, normal.z / d};
}

void save_binary_stl(const std::vector<Vec3>& positions, const std::string& path)
{
    std::cout << "save_binary_stl" << std::endl;
    std::ofstream file(path, std::ios::binary);

    char header[80] = {0};
    file.write(header, 80);

    uint32_t numTriangles = positions.size() / 3;
    file.write(reinterpret_cast<const char*>(&numTriangles), sizeof(uint32_t));

    for (size_t i = 0; i < positions.size(); i += 3) {
        Vec3 normal = normal_from_points(positions[i], positions[i + 1], positions[i + 2]);

        file.write(reinterpret_cast<const char*>(&normal), sizeof(Vec3));

        file.write(reinterpret_cast<const char*>(&positions[i]), sizeof(Vec3));
        file.write(reinterpret_cast<const char*>(&positions[i + 1]), sizeof(Vec3));
        file.write(reinterpret_cast<const char*>(&positions[i + 2]), sizeof(Vec3));

        uint16_t attr = 0;
        file.write(reinterpret_cast<const char*>(&attr), sizeof(uint16_t));
    }

    file.close();
}

void compute_normals(const std::vector<std::vector<float>>& data, int rows, int cols, int step, std::vector<std::vector<Vec3>>& normals)
{
    std::cout << "compute_normals" << std::endl;
    normals.resize(rows, std::vector<Vec3>(cols));

    for (int i = 1; i < rows - 1; i += step) {
        for (int j = 1; j < cols - 1; j += step) {
            float dzdx = (data[i][j + 1] - data[i][j - 1]) / 2.0f;
            float dzdy = (data[i + 1][j] - data[i - 1][j]) / 2.0f;
            Vec3 normal = {-dzdx, ALTITUDE_MULTIPLIER, -dzdy};
            float length = std::sqrt(normal.x * normal.x + normal.y * normal.y + normal.z * normal.z);
            normal = {normal.x / length, normal.y / length, normal.z / length};
            normals[i][j] = normal;
        }
    }
}

void save_normal_map(const std::vector<std::vector<Vec3>>& normals, int rows, int cols, const std::string& path)
{
    std::cout << "save_normal_map" << std::endl;
    uint8_t* image = new uint8_t[rows * cols * 3];
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            float nx = (normals[i][j].x + 1.0f) / 2.0f * 255.0f;
            float ny = (normals[i][j].y + 1.0f) / 2.0f * 255.0f;
            float nz = (normals[i][j].z + 1.0f) / 2.0f * 255.0f;
            image[(i * cols + j) * 3 + 0] = (uint8_t)nx;
            image[(i * cols + j) * 3 + 1] = (uint8_t)ny;
            image[(i * cols + j) * 3 + 2] = (uint8_t)nz;
        }
    }
    stbi_write_png(path.c_str(), cols, rows, 3, image, cols * 3);
    delete[] image;
}

std::vector<std::vector<float>> load(const std::string& path, int& rows, int& cols)
{
    std::cout << "load" << std::endl;
    GDALAllRegister();
    GDALDataset* dataset = (GDALDataset*)GDALOpen(path.c_str(), GA_ReadOnly);
    if (!dataset) {
        std::cerr << "Failed to load dataset" << std::endl;
        return {};
    }

    GDALRasterBand* band = dataset->GetRasterBand(1);
    rows = band->GetYSize();
    cols = band->GetXSize();

    std::vector<std::vector<float>> data(rows, std::vector<float>(cols));
    for (int i = 0; i < rows; ++i) {
        CPLErr err = band->RasterIO(GF_Read, 0, i, cols, 1, data[i].data(), cols, 1, GDT_Float32, 0, 0);
        if (err != CE_None) {
            std::cerr << "Error reading raster data" << std::endl;
            GDALClose(dataset);
            return {};
        }
    }

    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            if (data[i][j] < -1e9) {
                data[i][j] = NAN;
            }
        }
    }

    GDALClose(dataset);
    std::cout << "Loaded data with dimensions: " << rows << "x" << cols << std::endl;
    return data;
}

int main(int argc, char* argv[])
{
    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " <input_file>" << std::endl;
        return 1;
    }

    std::string path = argv[1];
    std::string name = path.substr(0, path.find_last_of('.'));

    int rows, cols;
    std::vector<std::vector<float>> data = load(path, rows, cols);
    if (data.empty())
        return 1;

    if (MESH_STEP) {
        std::vector<Vec3> positions;
        compute_mesh(data, rows, cols, MESH_STEP, positions);
        save_binary_stl(positions, name + ".stl");
    }

    if (NORMAL_STEP) {
        std::vector<std::vector<Vec3>> normals;
        compute_normals(data, rows, cols, NORMAL_STEP, normals);
        save_normal_map(normals, rows, cols, name + ".png");
    }

    return 0;
}
