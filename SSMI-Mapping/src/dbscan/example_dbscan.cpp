#include "dbscan_clustering.hpp"  // DBSCAN
#include "dbscan_point_cloud.hpp" // PointCloud

#include <chrono>   // std::chrono
#include <csignal>  // std::signal
#include <iostream> // std::cout
#include <memory>   // std::unique_ptr
#include <random>   // std::random_device
#include <variant>  // std::variant

using namespace clustering;

// Constants
constexpr static int NUMBER_OF_DIMENSIONS = 3;
constexpr static int NUMBER_OF_POINTS = 10000;
constexpr static int NUMBER_OF_ITERATIONS = 100;
constexpr static double NEAREST_NEIGHBOR_PROXIMITY = 0.5;
constexpr static std::int32_t MINIMUM_POINTS_TO_FORM_CLUSTER = 5;

using CoordinateType = float;
using PointType = DBSCANPoint<CoordinateType, NUMBER_OF_DIMENSIONS>;
using PointCloudType = DBSCANPointCloud<CoordinateType, NUMBER_OF_DIMENSIONS>;
using DBSCAN_Type = DBSCANClustering<CoordinateType, NUMBER_OF_DIMENSIONS>;
using IndicesT = std::vector<std::uint32_t>;
using ClusterT = std::unordered_map<std::int32_t, IndicesT>;
/// @brief Generate and return random 3D points
static inline void generateRandomData(PointCloudType &point_cloud)
{
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<CoordinateType> dist(-10.0, 10.0);

    point_cloud.clear();
    point_cloud.reserve(NUMBER_OF_POINTS);
    for (auto i = 0; i < NUMBER_OF_POINTS; ++i)
    {
        PointType point;
        for (auto j = 0; j < NUMBER_OF_DIMENSIONS; ++j)
        {
            point[j] = dist(gen);
        }
        point_cloud.push_back(point);
    }
}

namespace
{
volatile std::sig_atomic_t signal_status;
}

int main(int argc, const char **const argv)
{
    // Install signal handler
    std::signal(SIGINT, [](int signal) {
        signal_status = signal;
        std::cout << "Keyboard interrupt\n";
        std::exit(signal);
    });

    // Generate point cloud data
    PointCloudType point_cloud;
    generateRandomData(point_cloud);

    // Run several iterations of DBSCAN
    double average_time_dbscan = 0.0;
    // for (std::int32_t i = 0; i < NUMBER_OF_ITERATIONS; ++i)
    // {
        auto start_time = std::chrono::high_resolution_clock::now();

        DBSCAN_Type dbscan(point_cloud, NEAREST_NEIGHBOR_PROXIMITY, MINIMUM_POINTS_TO_FORM_CLUSTER);
        
        auto stop_time = std::chrono::high_resolution_clock::now();

        auto elapsed_seconds = (stop_time - start_time).count() / 1e9;
        average_time_dbscan += elapsed_seconds;
        std::cout << "Testing" << std::endl;
        std::cout << "Elapsed time (s): " << elapsed_seconds << std::endl;
    // }
    dbscan.formClusters();
    ClusterT cluster_indices = dbscan.getClusterIndices();

        for (const auto& cluster : cluster_indices) {
            std::cout << "Cluster ID: " << cluster.first << std::endl;
            std::cout << "Indices: ";
            for (const auto& index : cluster.second) {
                std::cout << index << " ";
            }
        std::cout << std::endl;
        }
    average_time_dbscan /= NUMBER_OF_ITERATIONS;
    std::cout << "Average time per loop (DBSCAN): " << average_time_dbscan << " seconds\n";

    return EXIT_SUCCESS;
}