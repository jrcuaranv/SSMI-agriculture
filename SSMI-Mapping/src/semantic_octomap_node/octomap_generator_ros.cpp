#include <semantic_octomap_node/octomap_generator_ros.h>
#include <pcl_ros/transforms.h>
#include <pcl_ros/impl/transforms.hpp>
#include <octomap_msgs/conversions.h>
#include <nav_msgs/OccupancyGrid.h>
#include <pcl/conversions.h>
#include <pcl/filters/voxel_grid.h>
#include <cmath>
#include <sstream>
#include <cstring> // For std::memcpy
#include <chrono>

OctomapGeneratorNode::OctomapGeneratorNode(ros::NodeHandle& nh): nh_(nh)
{
    // Initiate octree
    ROS_INFO("Semantic octomap generated!");
    toggle_color_service_ = nh_.advertiseService("toggle_use_semantic_color", &OctomapGeneratorNode::toggleUseSemanticColor, this);
    RLE_service_ = nh_.advertiseService("querry_RLE", &OctomapGeneratorNode::querry_RLE, this);

    reset();
    fullmap_pub_ = nh_.advertise<octomap_msgs::Octomap>("octomap_full", 1, true);
    occ_map_pub_ = nh_.advertise<nav_msgs::OccupancyGrid>("occupancy_map_2D", 1, true);
    pointcloud_sub_ = new message_filters::Subscriber<sensor_msgs::PointCloud2> (nh_, pointcloud_topic_, 5);
    tf_pointcloud_sub_ = new tf::MessageFilter<sensor_msgs::PointCloud2> (*pointcloud_sub_, tf_listener_, world_frame_id_, 5);
    tf_pointcloud_sub_->registerCallback(boost::bind(&OctomapGeneratorNode::insertCloudCallback, this, _1));
    reset_octomap_sub_ = nh_.subscribe("/reset_octomap", 1000, &OctomapGeneratorNode::resetCallback, this); //jrcv
    semantic_centroids_pub_ = nh.advertise<sensor_msgs::PointCloud2>("semantic_centroids", 10); //jrcv
    octomap_status_pub_ = nh.advertise<std_msgs::Float32>("octomap_status", 10); //jrcv
    coarse_pointcloud_pub_ = nh.advertise<sensor_msgs::PointCloud2> ("coarse_pointcloud", 1); //jrcv
}

OctomapGeneratorNode::~OctomapGeneratorNode() {}
/// Clear octomap and reset values to paramters from parameter server
void OctomapGeneratorNode::resetCallback(const std_msgs::Float32::ConstPtr& msg)
{
    float data = msg->data;
    if (data == 1){
        delete octomap_generator_;
        reset();
    }
}
void OctomapGeneratorNode::reset()
{
    octomap_generator_ = new OctomapGenerator<PCLSemantics, SemanticOctree>();
    nh_.getParam("/octomap/pointcloud_topic", pointcloud_topic_);
    nh_.getParam("/octomap/world_frame_id", world_frame_id_);
    nh_.getParam("/octomap/resolution", resolution_);
    nh_.getParam("/octomap/max_range", max_range_);
    nh_.getParam("/octomap/raycast_range", raycast_range_);
    nh_.getParam("/octomap/clamping_thres_min", clamping_thres_min_);
    nh_.getParam("/octomap/clamping_thres_max", clamping_thres_max_);
    nh_.getParam("/octomap/occupancy_thres", occupancy_thres_);
    nh_.getParam("/octomap/prob_hit", prob_hit_);
    nh_.getParam("/octomap/prob_miss", prob_miss_);
    nh_.getParam("/octomap/psi", psi_);
    nh_.getParam("/octomap/phi", phi_);
    nh_.getParam("/octomap/publish_2d_map", publish_2d_map);
    nh_.getParam("/octomap/min_ground_z", min_ground_z);
    nh_.getParam("/octomap/max_ground_z", max_ground_z);
    nh_.getParam("/octomap/save_path", output_dir_);

    octomap_generator_->setClampingThresMin(clamping_thres_min_);
    octomap_generator_->setClampingThresMax(clamping_thres_max_);
    octomap_generator_->setResolution(resolution_);
    octomap_generator_->setOccupancyThres(occupancy_thres_);
    octomap_generator_->setProbHit(prob_hit_);
    octomap_generator_->setProbMiss(prob_miss_);
    octomap_generator_->setPsi(psi_);
    octomap_generator_->setPhi(phi_);
    octomap_generator_->setRayCastRange(raycast_range_);
    octomap_generator_->setMaxRange(max_range_);
}

bool OctomapGeneratorNode::toggleUseSemanticColor(std_srvs::Empty::Request& request, std_srvs::Empty::Response& response)
{
    octomap_generator_->setUseSemanticColor(!octomap_generator_->isUseSemanticColor());
    if(octomap_generator_->isUseSemanticColor())
        ROS_INFO("Using semantic color");
    else
        ROS_INFO("Using rgb color");
    if (octomap_msgs::fullMapToMsg(*octomap_generator_->getOctree(), map_msg_))
        fullmap_pub_.publish(map_msg_);
    else
        ROS_ERROR("Error serializing OctoMap");
    return true;
}

bool OctomapGeneratorNode::querry_RLE(semantic_octomap::GetRLE::Request& request, semantic_octomap::GetRLE::Response& response)
{
    const octomap::point3d origin(request.origin.x, request.origin.y, request.origin.z);

    for (int i = 0; i < (int)request.endPoints.size(); ++i)
    {   
        const octomap::point3d endPoint(request.endPoints[i].x, request.endPoints[i].y, request.endPoints[i].z);
        semantic_octomap::RayRLE rayRLE_msg;
        if (octomap_generator_->get_ray_RLE(origin, endPoint, rayRLE_msg))
        {          
            response.RLE_list.push_back(rayRLE_msg);
        } 
    }
    
    return true;
}

void OctomapGeneratorNode::insertCloudCallback(const sensor_msgs::PointCloud2::ConstPtr& cloud_msg)
{
    // Voxel filter to down sample the point cloud
    // Create the filtering object
    pcl::PCLPointCloud2::Ptr cloud (new pcl::PCLPointCloud2 ());
    pcl_conversions::toPCL(*cloud_msg, *cloud);
    // Get tf transform
    tf::StampedTransform sensorToWorldTf;
    try
    {
        tf_listener_.lookupTransform(world_frame_id_, cloud_msg->header.frame_id, cloud_msg->header.stamp, sensorToWorldTf);
    }
    catch(tf::TransformException& ex)
    {
        ROS_ERROR_STREAM( "Transform error of sensor data: " << ex.what() << ", quitting callback");
        return;
    }
    // Transform coordinate
    Eigen::Matrix4f sensorToWorld;
    pcl_ros::transformAsMatrix(sensorToWorldTf, sensorToWorld);
    auto start_time = std::chrono::high_resolution_clock::now();
    octomap_generator_->insertPointCloud(cloud, sensorToWorld);
    auto stop_time = std::chrono::high_resolution_clock::now();
    auto elapsed_seconds = (stop_time - start_time).count() / 1e9;
    std::cout << "Time insertPointCloud: " << elapsed_seconds << std::endl;

    // publish semantic centroids
    std::vector<octomap::point3d> sem_centroids;
    octomap_generator_->compute_centroids(sem_centroids, output_dir_);
    

    // Publish full octomap
    start_time = std::chrono::high_resolution_clock::now();
    map_msg_.header.frame_id = world_frame_id_;
    map_msg_.header.stamp = cloud_msg->header.stamp;
    if (octomap_msgs::fullMapToMsg(*octomap_generator_->getOctree(), map_msg_))
        fullmap_pub_.publish(map_msg_);
    else
        ROS_ERROR("Error serializing full OctoMap");
    
    stop_time = std::chrono::high_resolution_clock::now();
    elapsed_seconds = (stop_time - start_time).count() / 1e9;
    std::cout << "Time publish full octomap: " << elapsed_seconds << std::endl;

    // Publish 2D occupancy map
    if (publish_2d_map)
        publish2DOccupancyMap(octomap_generator_->getOctree(), cloud_msg->header.stamp, world_frame_id_);
    
    // publish Semantic Centroids
    publishSemanticCentroids(octomap_generator_->getOctree(), cloud_msg->header.stamp, world_frame_id_);
    // publish octomap status done
    std_msgs::Float32 oct_status_msg;
    oct_status_msg.data = 1.0;
    octomap_status_pub_.publish(oct_status_msg);

}

void OctomapGeneratorNode::publishSemanticCentroids(const SemanticOctree* octomap,
                                                 const ros::Time& stamp,
                                                 const std::string& frame_id)
{
    std::vector<octomap::point3d> points;
    points = octomap->getSemanticCentroids();
    // Create a PointCloud2 message
    sensor_msgs::PointCloud2::Ptr point_cloud (new sensor_msgs::PointCloud2());
    
    point_cloud->header.frame_id = frame_id; // Set the frame ID appropriately
    point_cloud->header.stamp = stamp;
    point_cloud->height = 1;
    point_cloud->width = points.size();
    point_cloud->fields.resize(3);
    point_cloud->fields[0].name = "x";
    point_cloud->fields[0].offset = 0;
    point_cloud->fields[0].datatype = sensor_msgs::PointField::FLOAT32;
    point_cloud->fields[0].count = 1;
    point_cloud->fields[1].name = "y";
    point_cloud->fields[1].offset = 4;
    point_cloud->fields[1].datatype = sensor_msgs::PointField::FLOAT32;
    point_cloud->fields[1].count = 1;
    point_cloud->fields[2].name = "z";
    point_cloud->fields[2].offset = 8;
    point_cloud->fields[2].datatype = sensor_msgs::PointField::FLOAT32;
    point_cloud->fields[2].count = 1;
    point_cloud->is_bigendian = false;
    point_cloud->point_step = 12;
    point_cloud->row_step = 12 * points.size();
    point_cloud->is_dense = true;
    point_cloud->data.resize(point_cloud->row_step);
    auto float_data_iter = reinterpret_cast<float*>(point_cloud->data.data());
    // for (const auto& point : points) {
    //     for (float value : point) {
    //         *float_data_iter = value;
    //         ++float_data_iter;
    //     }
    // }
    for (const auto& point : points) {
        *float_data_iter = point.x();
        ++float_data_iter;
        *float_data_iter = point.y();
        ++float_data_iter;
        *float_data_iter = point.z();
        ++float_data_iter;
    }


    // Publish the PointCloud2 message
    semantic_centroids_pub_.publish(*point_cloud);
    
}
void OctomapGeneratorNode::publish2DOccupancyMap(const SemanticOctree* octomap,
                                                 const ros::Time& stamp,
                                                 const std::string& frame_id)
{
  // get dimensions of octree
  double minX, minY, minZ, maxX, maxY, maxZ;
  octomap->getMetricMin(minX, minY, minZ);
  octomap->getMetricMax(maxX, maxY, maxZ);
  octomap::point3d minPt = octomap::point3d(minX, minY, minZ);

  unsigned int tree_depth = octomap->getTreeDepth();

  octomap::OcTreeKey paddedMinKey = octomap->coordToKey(minPt);

  nav_msgs::OccupancyGrid::Ptr occupancy_map (new nav_msgs::OccupancyGrid());

  unsigned int width, height;
  double res;

  unsigned int ds_shift = tree_depth-16;

  occupancy_map->header.stamp = stamp;
  occupancy_map->header.frame_id = frame_id;
  occupancy_map->info.resolution = res = octomap->getNodeSize(16);
  occupancy_map->info.width = width = (maxX-minX) / res + 1;
  occupancy_map->info.height = height = (maxY-minY) / res + 1;
  occupancy_map->info.origin.position.x = minX  - (res / (float)(1<<ds_shift) ) + res;
  occupancy_map->info.origin.position.y = minY  - (res / (float)(1<<ds_shift) );

  occupancy_map->data.clear();
  occupancy_map->data.resize(width*height, -1);

    // traverse all leafs in the tree:
  unsigned int treeDepth = std::min<unsigned int>(16, octomap->getTreeDepth());
  for (typename SemanticOctree::iterator it = octomap->begin(treeDepth), end = octomap->end(); it != end; ++it)
  {
  
    double node_z = it.getZ();
    double node_half_side = pow(it.getSize(), 1/3) / 2;
    double top_side = node_z + node_half_side;
    double bottom_side = node_z - node_half_side;
    
    if((bottom_side >= min_ground_z && bottom_side <= max_ground_z) ||
       (top_side >= min_ground_z && top_side <= max_ground_z) ||
       (bottom_side <= min_ground_z && top_side >= max_ground_z))
    {
      bool occupied = octomap->isNodeOccupied(*it);
      int intSize = 1 << (16 - it.getDepth());

      octomap::OcTreeKey minKey=it.getIndexKey();

      for (int dx = 0; dx < intSize; dx++)
      {
        for (int dy = 0; dy < intSize; dy++)
        {
          int posX = std::max<int>(0, minKey[0] + dx - paddedMinKey[0]);
          posX>>=ds_shift;

          int posY = std::max<int>(0, minKey[1] + dy - paddedMinKey[1]);
          posY>>=ds_shift;

          int idx = width * posY + posX;

          if (occupied)
            occupancy_map->data[idx] = 100;
          else if (occupancy_map->data[idx] == -1)
          {
            occupancy_map->data[idx] = 0;
          }

        }
      }
    }
  }

  occ_map_pub_.publish(*occupancy_map);
}

bool OctomapGeneratorNode::save(const char* filename) const
{
    return octomap_generator_->save(filename);
}

int main(int argc, char** argv)
{
    ros::init(argc, argv, "octomap_generator");
    ros::NodeHandle nh;
    OctomapGeneratorNode octomapGeneratorNode(nh);
    ros::spin();
    std::string save_path;
    nh.getParam("/octomap/save_path", save_path);
    std::string map_file_path = save_path + "map.ot";

    octomapGeneratorNode.save(map_file_path.c_str());
    ROS_INFO("OctoMap saved.");
    return 0;
}
