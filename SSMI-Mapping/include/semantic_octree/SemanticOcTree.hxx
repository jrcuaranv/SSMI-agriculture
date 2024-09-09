#ifndef SEMANTIC_OCTOMAP_SEMANTICOCTREE_HXX
#define SEMANTIC_OCTOMAP_SEMANTICOCTREE_HXX



namespace octomap {

    static float computeDistance(const octomap::point3d& p1, const octomap::point3d& p2) {
    float dx = p2.x() - p1.x();
    float dy = p2.y() - p1.y();
    float dz = p2.z() - p1.z();
    return std::sqrt(dx * dx + dy * dy + dz * dz);
    }

    // Tree implementation  --------------------------------------
    template<class SEMANTICS>
    SemanticOcTree<SEMANTICS>::SemanticOcTree(double resolution)
            : OccupancyOcTreeBase<SemanticOcTreeNode<SEMANTICS> >(this->resolution),
                    phiTree(), psiTree(), maxLogOddsTree(), minLogOddsTree(), minOccupancyLogOdds(), semanticCentroids()
    {
        semanticOcTreeMemberInit.ensureLinking();
    };

    template<class SEMANTICS>
    bool SemanticOcTree<SEMANTICS>::pruneNode(SemanticOcTreeNode<SEMANTICS>* node) {
        // Same as ColorOcTree
        if (!isNodeCollapsible(node))
            return false;

        // Set value to children's values (all assumed equal)
        node->copyData(*(this->getNodeChild(node, 0)));
        
        if (node->isColorSet()) // TODO check
            node->setColor(node->getAverageChildColor());
            

        // Delete children
        for (unsigned int i=0;i<8;i++) {
            this->deleteNodeChild(node, i);
        }
        delete[] node->children;
        node->children = NULL;

        return true;
    }

    template<class SEMANTICS>
    bool SemanticOcTree<SEMANTICS>::isNodeCollapsible(const SemanticOcTreeNode<SEMANTICS>* node) const
    {
        // All children must exist, must not have children of
        // their own and have same occupancy probability and same log-odds vector
        if(!this->nodeChildExists(node, 0))
            return false;
        const SemanticOcTreeNode<SEMANTICS>* firstChild = this->getNodeChild(node, 0);
        if(this->nodeHasChildren(firstChild))
            return false;
        bool firstChildFree = firstChild->getValue() <= this->minOccupancyLogOdds;
        for (unsigned int i = 1; i<8; i++) {
            // Compare nodes using their occupancy and log-odds vector, ignoring color for pruning
            if (!this->nodeChildExists(node, i) || this->nodeHasChildren(this->getNodeChild(node, i))
                || !(this->getNodeChild(node, i)->getValue() == firstChild->getValue())
                || !(this->getNodeChild(node, i)->getSemantics() == firstChild->getSemantics() || firstChildFree))
                return false;
        }
        return true;
    }

    template<class SEMANTICS>
    void SemanticOcTree<SEMANTICS>::setUseSemanticColor(bool use)
    {
        // Traverse all tree nodes
        for(typename SemanticOcTree<SEMANTICS>::tree_iterator it = this->begin_tree(), end=this->end_tree(); it!= end; ++it)
            it->use_semantic_color = use;
    }

    template<class SEMANTICS>
    SemanticOcTreeNode<SEMANTICS>* SemanticOcTree<SEMANTICS>::setNodeColor(const OcTreeKey& key,
                                                                           uint8_t r,
                                                                           uint8_t g,
                                                                           uint8_t b)
    {
        SemanticOcTreeNode<SEMANTICS>* n = this->search (key);
        if (n != 0) {
            n->setColor(r, g, b);
        }
        return n;
    }

    template<class SEMANTICS>
    SemanticOcTreeNode<SEMANTICS>* SemanticOcTree<SEMANTICS>::averageNodeColor(SemanticOcTreeNode<SEMANTICS>* node, uint8_t r,
                                                                               uint8_t g, uint8_t b)
    {
        if (node != 0) {
            if (node->isColorSet()) {
                ColorOcTreeNode::Color prev_color = node->getColor();
                node->setColor((prev_color.r + r)/2, (prev_color.g + g)/2, (prev_color.b + b)/2);
            }
            else {
                node->setColor(r, g, b);
            }
        }
        return node;
    }

    template<class SEMANTICS>
    SemanticOcTreeNode<SEMANTICS>* SemanticOcTree<SEMANTICS>::updateNodeSemantics(SemanticOcTreeNode<SEMANTICS>* node,
                                                                                  ColorOcTreeNode::Color obs)
    {
        if (node != 0)
        {
            
            SEMANTICS sem;
            
            if (node->isSemanticsSet())
            {
                sem = SEMANTICS::fuseObs(node->semantics, obs, phiTree, psiTree, maxLogOddsTree, minLogOddsTree);
            }
            else
            {
                sem = SEMANTICS::initSemantics(obs, node->getLogOdds(), phiTree, psiTree, maxLogOddsTree, minLogOddsTree);
            }
            
            float logOddsValue = SEMANTICS::getOccFromSem(sem);
            
            node->setSemantics(sem);
            node->setLogOdds(logOddsValue);
        }
        
        return node;
    }
    
    template<class SEMANTICS>
    SemanticOcTreeNode<SEMANTICS>* SemanticOcTree<SEMANTICS>::updateFreeNode(SemanticOcTreeNode<SEMANTICS>* node)
    {
        if(node->isSemanticsSet()) // if node is semantic, and observation is free
        {
            SEMANTICS sem = SEMANTICS::fuseObsFree(node->semantics, phiTree, minLogOddsTree);
            float logOddsValue = SEMANTICS::getOccFromSem(sem);
            
            node->setSemantics(sem);
            node->setLogOdds(logOddsValue);
        } else // if node is free and observation is free
        {
            float nodeLogOdds = node->getLogOdds() + phiTree;
            if(nodeLogOdds < minOccupancyLogOdds)
                nodeLogOdds = minOccupancyLogOdds;
            
            node->setLogOdds(nodeLogOdds);
        }
        
        return node;
    }

    template<class SEMANTICS>
    SemanticOcTreeNode<SEMANTICS>* SemanticOcTree<SEMANTICS>::updateNode(const OcTreeKey& key, bool occupied,
                                                                         const ColorOcTreeNode::Color& class_obs,
                                                                         const ColorOcTreeNode::Color& color_obs,
                                                                         bool lazy_eval)
    {
        // early abort (no change will happen).
        // may cause an overhead in some configuration, but more often helps
        SemanticOcTreeNode<SEMANTICS>* leaf = this->search(key);
        // no change: node already at threshold
        if (leaf && !checkNeedsUpdate(leaf, occupied, class_obs))
        {
          return leaf;
        }

        bool createdRoot = false;
        if (this->root == NULL){
          this->root = new SemanticOcTreeNode<SEMANTICS>();
          this->tree_size++;
          createdRoot = true;
        }

        return updateNodeRecurs(this->root, createdRoot, key, 0, occupied, class_obs, color_obs, lazy_eval);
    }

    template<class SEMANTICS>
    SemanticOcTreeNode<SEMANTICS>* SemanticOcTree<SEMANTICS>::updateNode(float x, float y, float z, bool occupied,
                                                                         const ColorOcTreeNode::Color& class_obs,
                                                                         const ColorOcTreeNode::Color& color_obs,
                                                                         bool lazy_eval)
    {
        OcTreeKey key;
        if (!this->coordToKeyChecked(x, y, z, key))
            return NULL;
        return updateNode(key, occupied, class_obs, color_obs, lazy_eval);
    }

    template<class SEMANTICS>
    SemanticOcTreeNode<SEMANTICS>* SemanticOcTree<SEMANTICS>::updateNodeRecurs(SemanticOcTreeNode<SEMANTICS>* node, bool node_just_created, const OcTreeKey& key,
                                                                               unsigned int depth, bool occupied, const ColorOcTreeNode::Color& class_obs,
                                                                               const ColorOcTreeNode::Color& color_obs, bool lazy_eval)
    {
        bool created_node = false;

        assert(node);

        // follow down to last level
        if (depth < this->tree_depth) {
          unsigned int pos = computeChildIdx(key, this->tree_depth -1 - depth);
          if (!this->nodeChildExists(node, pos)) {
            // child does not exist, but maybe it's a pruned node?
            if (!this->nodeHasChildren(node) && !node_just_created ) {
              // current node does not have children AND it is not a new node 
              // -> expand pruned node
              this->expandNode(node);
            }
            else {
              // not a pruned node, create requested child
              this->createNodeChild(node, pos);
              created_node = true;
            }
          }

          if (lazy_eval)
            return updateNodeRecurs(this->getNodeChild(node, pos), created_node, key, depth+1, occupied, class_obs, color_obs, lazy_eval);
          else {
            SemanticOcTreeNode<SEMANTICS>* retval = updateNodeRecurs(this->getNodeChild(node, pos), created_node, key, depth+1, occupied, class_obs, color_obs, lazy_eval);
            // prune node if possible, otherwise set own probability
            // note: combining both did not lead to a speedup!
            if (this->pruneNode(node)){
              // return pointer to current parent (pruned), the just updated node no longer exists
              retval = node;
            } else {
              node->updateSemanticsChildren();
              node->updateColorChildren();
            }

            return retval;
          }
        }

        // at last level, update node, end of recursion
        else {
          if (this->use_change_detection) {
            bool occBefore = this->isNodeOccupied(node);
            updateNodeLogOdds(node, occupied, class_obs, color_obs);

            if (node_just_created){  // new node
              this->changed_keys.insert(std::pair<OcTreeKey,bool>(key, true));
            } else if (occBefore != this->isNodeOccupied(node)) {  // occupancy changed, track it
              KeyBoolMap::iterator it = this->changed_keys.find(key);
              if (it == this->changed_keys.end())
                this->changed_keys.insert(std::pair<OcTreeKey,bool>(key, false));
              else if (it->second == false)
                this->changed_keys.erase(it);
            }
          } else {
            updateNodeLogOdds(node, occupied, class_obs, color_obs);
          }
          return node;
        }
    }

    template<class SEMANTICS>
    void SemanticOcTree<SEMANTICS>::updateNodeLogOdds(SemanticOcTreeNode<SEMANTICS>* node, bool occupied,
                                                      const ColorOcTreeNode::Color& class_obs,
                                                      const ColorOcTreeNode::Color& color_obs)
    {
        // node color update
        averageNodeColor(node, color_obs.r, color_obs.g, color_obs.b);
        // node semantics and occupancy update
        if(occupied)
        {
            updateNodeSemantics(node, class_obs);
        } else {
            updateFreeNode(node);
        }
    }

    template<class SEMANTICS>
    void SemanticOcTree<SEMANTICS>::insertPointCloud(const Pointcloud& scan, const octomap::point3d& sensor_origin,
                                                     double maxrange, bool discretize)
    {
        KeySet free_cells, occupied_cells;
        if (discretize)
          this->computeDiscreteUpdate(scan, sensor_origin, free_cells, occupied_cells, maxrange);
        else
          this->computeUpdate(scan, sensor_origin, free_cells, occupied_cells, maxrange);

        // Insert data into tree only for the free cells, occupied cells will be taken care of separately  -----------------------
        for (KeySet::iterator it = free_cells.begin(); it != free_cells.end(); ++it) {
          updateNode(*it, false);
        }
        // std::vector<octomap::point3d> sem_centroids;
        // if (this->compute_centroids(sem_centroids)){
        //     this->setSemanticCentroids(sem_centroids);
        // }
    }

    template<class SEMANTICS>
    void SemanticOcTree<SEMANTICS>::updateInnerOccupancy() {
        this->updateInnerOccupancyRecurs(this->root, 0);
    }

    template<class SEMANTICS>
    void SemanticOcTree<SEMANTICS>::updateInnerOccupancyRecurs(SemanticOcTreeNode<SEMANTICS>* node, unsigned int depth) {
        // Only recurse and update for inner nodes:
        if (this->nodeHasChildren(node)){
            // Return early for last level:
            if (depth < this->tree_depth){
                for (unsigned int i=0; i<8; i++) {
                    if (this->nodeChildExists(node, i)) {
                        updateInnerOccupancyRecurs(this->getNodeChild(node, i), depth+1);
                    }
                }
            }
            // Update occupancy, semantics and color for inner nodes
            node->updateSemanticsChildren();
            node->updateColorChildren();
        }
    }
    
    template<class SEMANTICS>
    bool SemanticOcTree<SEMANTICS>::checkNeedsUpdate(const SemanticOcTreeNode<SEMANTICS>* node, bool occupied, const ColorOcTreeNode::Color& class_obs)
    {   
        SEMANTICS sem = node->getSemantics();
        
        if((!occupied && node->getLogOdds() <= minOccupancyLogOdds) ||
           (sem.data[0].color == class_obs &&
            sem.data[0].logOdds >= maxLogOddsTree && sem.data[1].logOdds <= minLogOddsTree && 
            sem.data[2].logOdds <= minLogOddsTree && sem.others <= minLogOddsTree))
        {
            return false;
        } else {
            return true;
        }
    }



    template<class SEMANTICS>
    bool SemanticOcTree<SEMANTICS>::compute_centroids(std::vector<octomap::point3d>& sem_centroids, std::string& output_dir)
    {
        auto start_time = std::chrono::high_resolution_clock::now();

        // Vectors to store coordinates and occupancy values
        std::vector<octomap::point3d> coordinates;
        // Iterate through all the nodes of the OctoMap
        for(typename SemanticOcTree<SEMANTICS>::leaf_iterator it = this->begin_leafs(), end=this->end_leafs(); it!= end; ++it)
        {
            
            SEMANTICS leaf_sem = it->getSemantics();
            ColorOcTreeNode::Color leaf_color =  ColorOcTreeNode::Color(leaf_sem.data[0].color.r,
                                                                    leaf_sem.data[0].color.g,
                                                                    leaf_sem.data[0].color.b);
            if (leaf_color == ColorOcTreeNode::Color(255,0,0) && leaf_sem.data[0].logOdds > 1.0)
            {
                octomap::point3d coord = it.getCoordinate();
                coordinates.push_back(coord);
            }
            

        }
        if (coordinates.size()  == 0)
        {
            std::cout << "no fruit semantic points were found" << std::endl;
            return false;
        }

        // saving 3D points in a .txt. file
        std::string pointcloud_path = output_dir + "semantic_points.txt";
        std::ofstream outfile(pointcloud_path);

        if (outfile.is_open()) {
            // Write each point's coordinates to the file
            for (const auto& point : coordinates) {
                outfile << point.x() << " " << point.y() << " " << point.z() << std::endl;
            }

            // Close the file
            outfile.close();
            std::cout << "Points written to points.txt successfully." << std::endl;
        } else {
            std::cerr << "Error opening file for writing." << std::endl;
        }
       
        clustering::DBSCANPointCloud<float, 3> dbscan_pointcloud;

        // Iterate over the vector of 3D points and add them to the point cloud
        for (const auto& coord : coordinates) {
            std::array<float, 3> dbscan_point;
            dbscan_point[0] = coord.x();
            dbscan_point[1] = coord.y();
            dbscan_point[2] = coord.z();
            dbscan_pointcloud.push_back(dbscan_point);

        }

        std::cout << "Point cloud created with " << dbscan_pointcloud.size() << " points." << std::endl;
        // Set up the parameters for clustering
        // float cluster_tolerance = 0.03;   // Cluster tolerance (in meters)
        float epsilon = 0.05; // neighborhood threshold for dbscan
        int min_cluster_size = 3;//7;       // Minimum number of points in a cluster
        int max_cluster_size = 10000;    // Maximum number of points in a cluster


        clustering::DBSCANClustering<float, 3> dbscan(dbscan_pointcloud, epsilon, min_cluster_size);
        dbscan.formClusters();
        std::unordered_map<std::int32_t, std::vector<std::uint32_t>> cluster_indices = dbscan.getClusterIndices();
        
        if (cluster_indices.size()  == 0)
        {
            std::cout << "no clusters were found" << std::endl;
            return false;
        }
        // Output the clusters
        
        std::cout << "Number of clusters: " << cluster_indices.size() << std::endl;
        double total_clusters_volume = 0;

        for (const auto& cluster : cluster_indices) {
            std::cout << "Cluster ID: " << cluster.first << std::endl;
            // std::array<float, 3> db_point;
            float x_sum = 0.0;
            float y_sum = 0.0;
            float z_sum = 0.0;
            std::vector<float> centroid;
            pcl::PointCloud<pcl::PointXYZ>::Ptr cluster_cloud(new pcl::PointCloud<pcl::PointXYZ>);
            
            
            // Computing centroid of clusters
            for (const auto& index : cluster.second) {
                x_sum = x_sum + dbscan_pointcloud[index][0];
                y_sum = y_sum + dbscan_pointcloud[index][1];
                z_sum = z_sum + dbscan_pointcloud[index][2];
                pcl::PointXYZ pcl_point;
                pcl_point.x = dbscan_pointcloud[index][0];
                pcl_point.y = dbscan_pointcloud[index][1];
                pcl_point.z = dbscan_pointcloud[index][2];
                cluster_cloud->points.push_back(pcl_point);
            }

            centroid.push_back(x_sum/cluster.second.size());
            centroid.push_back(y_sum/cluster.second.size());
            centroid.push_back(z_sum/cluster.second.size());
            octomap::point3d cent(centroid[0],centroid[1],centroid[2]);
            sem_centroids.push_back(cent);

            // Compute convex hull of the cluster

            pcl::ConvexHull<pcl::PointXYZ> hull;
            hull.setInputCloud(cluster_cloud);

            pcl::PointCloud<pcl::PointXYZ> points;
            std::vector<pcl::Vertices > polygons;
            hull.setComputeAreaVolume (true);
            hull.reconstruct(points, polygons);
            double volume = hull.getTotalVolume();
            total_clusters_volume = total_clusters_volume + volume;
            std::cout << "Cluster "<< cluster.first << " Volume: "<< volume << std::endl;

        }

        std::cout << "Total clusters volume: " << total_clusters_volume << std::endl;

        auto stop_time = std::chrono::high_resolution_clock::now();
        auto elapsed_seconds = (stop_time - start_time).count() / 1e9;
        std::cout << "Time clustering dbscan (s): " << elapsed_seconds << std::endl;

        this->setSemanticCentroids(sem_centroids);
        return true;
    }    

    template<class SEMANTICS>
    bool SemanticOcTree<SEMANTICS>::get_ray_RLE(const octomap::point3d& origin,
                                                const octomap::point3d& end,
                                                semantic_octomap::RayRLE& rayRLE_msg)
    {
        std::vector<octomap::point3d> sem_centroids;
        sem_centroids = this->getSemanticCentroids();
        KeyRay* keyray = &(this->keyrays.at(0));
        float offset_irrelevant_semantics = 5.0;
        float offset_main_semantics = 2.0;
        double node_occupancy_prob; //= currNode->getOccupancy();
        float node_occupancy_logodds; //= currNode->getLogOdds();
        if (this->computeRayKeys(origin, end, *keyray))
        {   
            semantic_octomap::LE le_msg_prev;
            le_msg_prev.le.push_back(1);
            SemanticOcTreeNode<SEMANTICS>* prevNode = this->search(*(keyray->begin()));
            if (!prevNode) //if previous node is unknown
            {
                for (int i = 0; i < 4; i++)
                {
                    le_msg_prev.le.push_back(-0.1);
                }
            } else if (!prevNode->isSemanticsSet()) { // if previous node is free
                float l = prevNode->getLogOdds() - log(4.);
                for (int i = 0; i < 4; i++)
                {
                    le_msg_prev.le.push_back(l);
                }
            } else { // when previous node is semantic
                SEMANTICS nodeSem = prevNode->getSemantics();
                for (int i = 0; i < 3; i++)
                {
                    le_msg_prev.le.push_back(nodeSem.data[i].logOdds);
                }
                le_msg_prev.le.push_back(nodeSem.others);
            }

            for (KeyRay::iterator it = keyray->begin()+1; it != keyray->end(); ++it)
            {
                octomap::point3d node_coord = this->keyToCoord(*it);
                
                SemanticOcTreeNode<SEMANTICS>* currNode = this->search(*it);
                
                
                if (!currNode && !prevNode) // When both previous and current nodes are unknown
                {
                    le_msg_prev.le[0] += 1;
                    if (it == keyray->end()-1){
                            le_msg_prev.le.push_back(127); // just a random color for unknown nodes, jrcv
                            le_msg_prev.le.push_back(127);
                            le_msg_prev.le.push_back(127);
                            le_msg_prev.le.push_back(node_coord.x());
                            le_msg_prev.le.push_back(node_coord.y());
                            le_msg_prev.le.push_back(node_coord.z());
                            le_msg_prev.le.push_back(0.5);
                            
                            rayRLE_msg.le_list.push_back(le_msg_prev); // pushing unknown nodes data

                        }
                } else if (currNode && prevNode) { // if previous and current nodes are known (free or semantic)
                    if (*currNode == *prevNode)
                    {
                        le_msg_prev.le[0] += 1;
                    } else { // if current node is a different category
                        
                        SEMANTICS nodeSem = prevNode->getSemantics();
                        ColorOcTreeNode::Color nodeColor = prevNode->getColor();
                        
                        // getting additional data for debugging
                        le_msg_prev.le.push_back((float)nodeSem.data[0].color.r);
                        le_msg_prev.le.push_back((float)nodeSem.data[0].color.g);
                        le_msg_prev.le.push_back((float)nodeSem.data[0].color.b);
                        le_msg_prev.le.push_back(node_coord.x());
                        le_msg_prev.le.push_back(node_coord.y());
                        le_msg_prev.le.push_back(node_coord.z());
                        le_msg_prev.le.push_back(prevNode->getOccupancy());
                        
                        if (nodeSem.data[0].logOdds > 0.0){ //skipping free nodes and semantic nodes with negative logodds as they are actually  probabilistically free
                            rayRLE_msg.le_list.push_back(le_msg_prev);    
                        }
                        le_msg_prev.le.clear();
                        le_msg_prev.le.push_back(1);
                        
                        if (!currNode->isSemanticsSet()) {  //if current node is free
                            float l = currNode->getLogOdds(); // - log(4.); //what is log(4) doing here?
                            for (int i = 0; i < 4; i++)
                            {
                                le_msg_prev.le.push_back(l);
                            }
                        } else {
                            SEMANTICS nodeSem = currNode->getSemantics();
                            for (int i = 0; i < 3; i++)
                            {
                                le_msg_prev.le.push_back(nodeSem.data[i].logOdds);
                            }
                            le_msg_prev.le.push_back(nodeSem.others);
                        }
                        prevNode = currNode;
                    }
                } else {// when current node OR previous node are unknown, jrcv
                    
                    if(prevNode){
                        SEMANTICS nodeSem = prevNode->getSemantics();
                        ColorOcTreeNode::Color nodeColor = prevNode->getColor();
                        // getting additional data for debugging
                        le_msg_prev.le.push_back((float)nodeSem.data[0].color.r);
                        le_msg_prev.le.push_back((float)nodeSem.data[0].color.g);
                        le_msg_prev.le.push_back((float)nodeSem.data[0].color.b);
                        le_msg_prev.le.push_back(node_coord.x());
                        le_msg_prev.le.push_back(node_coord.y());
                        le_msg_prev.le.push_back(node_coord.z());
                        le_msg_prev.le.push_back(prevNode->getOccupancy());
                        
                        if (nodeSem.data[0].logOdds > 0.0){ //skipping free nodes and semantic nodes with negative logodds as they are actually  probabilistically free
                            rayRLE_msg.le_list.push_back(le_msg_prev);    
                        }
                        le_msg_prev.le.clear();
                    }
                    else { // if previous node is unknown and current node is known
                        le_msg_prev.le.push_back(127); // just a random color for unknown nodes for debugging, jrcv
                        le_msg_prev.le.push_back(127);
                        le_msg_prev.le.push_back(127);
                        le_msg_prev.le.push_back(node_coord.x());
                        le_msg_prev.le.push_back(node_coord.y());
                        le_msg_prev.le.push_back(node_coord.z());
                        le_msg_prev.le.push_back(0.5);
                        
                        rayRLE_msg.le_list.push_back(le_msg_prev); // pushing unknown nodes data
                    }
                    
                    le_msg_prev.le.clear();
                    le_msg_prev.le.push_back(1);
                    
                    if (!currNode) // if current node is the unknown
                    {
                        for (int i = 0; i < 4; i++)
                        {
                            le_msg_prev.le.push_back(-0.1);
                        }
                    } else { // if current node is known and privious unknown
                        
                        if (!currNode->isSemanticsSet()) { // if current node is free?
                            float l = currNode->getLogOdds(); // - log(4.);
                            for (int i = 0; i < 4; i++)
                            {
                                le_msg_prev.le.push_back(l);
                            }
                        } else {
                            SEMANTICS nodeSem = currNode->getSemantics();
                            for (int i = 0; i < 3; i++)
                            {
                                le_msg_prev.le.push_back(nodeSem.data[i].logOdds);
                            }
                            le_msg_prev.le.push_back(nodeSem.others);
                        }
                    }
                    
                    prevNode = currNode;
                }
                
            }
            //rayRLE_msg.le_list.push_back(le_msg_prev);
            
            return true;
        } else {
            return false;
        }
    }

    template<class SEMANTICS>
    typename SemanticOcTree<SEMANTICS>::StaticMemberInitializer SemanticOcTree<SEMANTICS>::semanticOcTreeMemberInit;

} // end namespace

#endif //SEMANTIC_OCTOMAP_SEMANTICOCTREE_HXX
