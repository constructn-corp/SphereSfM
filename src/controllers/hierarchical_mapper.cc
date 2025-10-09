// Copyright (c) 2023, ETH Zurich and UNC Chapel Hill.
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
//
//     * Redistributions in binary form must reproduce the above copyright
//       notice, this list of conditions and the following disclaimer in the
//       documentation and/or other materials provided with the distribution.
//
//     * Neither the name of ETH Zurich and UNC Chapel Hill nor the names of
//       its contributors may be used to endorse or promote products derived
//       from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDERS OR CONTRIBUTORS BE
// LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
// CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
// SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
// INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
// CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
// ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
// POSSIBILITY OF SUCH DAMAGE.
//
// Author: Johannes L. Schoenberger (jsch-at-demuc-dot-de)

#include "controllers/hierarchical_mapper.h"

#include "base/scene_clustering.h"
#include "util/misc.h"
#include <sys/stat.h>  
#include <sys/types.h> 
#include <fstream>

namespace colmap {
namespace {

// Helper function to read cluster files from a list file
std::vector<std::vector<std::string>> ReadCustomClusterFiles(const std::string& cluster_list_path) {
    std::vector<std::vector<std::string>> clusters;
    
    std::ifstream cluster_list_file(cluster_list_path);
    if (!cluster_list_file.is_open()) {
        std::cerr << "ERROR: Could not open cluster list file: " << cluster_list_path << std::endl;
        return clusters;
    }
    
    std::string cluster_file_path;
    int cluster_count = 0;
    while (std::getline(cluster_list_file, cluster_file_path)) {
        if (cluster_file_path.empty() || cluster_file_path[0] == '#') {
            continue;  // Skip empty lines and comments
        }
        
        std::ifstream cluster_file(cluster_file_path);
        if (!cluster_file.is_open()) {
            std::cerr << "WARNING: Could not open cluster file: " << cluster_file_path << std::endl;
            continue;
        }
        
        std::vector<std::string> cluster_images;
        std::string image_name;
        while (std::getline(cluster_file, image_name)) {
            if (!image_name.empty() && image_name[0] != '#') {
                // Remove any trailing whitespace or path separators
                image_name.erase(image_name.find_last_not_of(" \t\r\n") + 1);
                cluster_images.push_back(image_name);
            }
        }
        
        if (!cluster_images.empty()) {
            clusters.push_back(cluster_images);
            cluster_count++;
            std::cout << "Loaded cluster " << cluster_count << " with " << cluster_images.size() 
                      << " images from " << cluster_file_path << std::endl;
        }
        cluster_file.close();
    }
    
    cluster_list_file.close();
    return clusters;
}

// Helper function to convert image names to image IDs
std::vector<image_t> ConvertNamesToIds(const std::vector<std::string>& image_names, 
                                       const std::unordered_map<std::string, image_t>& name_to_id_map) {
    std::vector<image_t> image_ids;
    image_ids.reserve(image_names.size());
    
    for (const auto& name : image_names) {
        auto it = name_to_id_map.find(name);
        if (it != name_to_id_map.end()) {
            image_ids.push_back(it->second);
        } else {
            std::cerr << "WARNING: Image not found in database: " << name << std::endl;
        }
    }
    
    return image_ids;
}


void MergeClusters(
    const SceneClustering::Cluster& cluster,
    std::unordered_map<const SceneClustering::Cluster*, ReconstructionManager>*
        reconstruction_managers) {
  // Extract all reconstructions from all child clusters.
  std::vector<Reconstruction*> reconstructions;
  for (const auto& child_cluster : cluster.child_clusters) {
    if (!child_cluster.child_clusters.empty()) {
      MergeClusters(child_cluster, reconstruction_managers);
    }

    auto& reconstruction_manager = reconstruction_managers->at(&child_cluster);
    for (size_t i = 0; i < reconstruction_manager.Size(); ++i) {
      reconstructions.push_back(&reconstruction_manager.Get(i));
    }
  }

  // Try to merge all child cluster reconstruction.
  while (reconstructions.size() > 1) {
    bool merge_success = false;
    for (size_t i = 0; i < reconstructions.size(); ++i) {
      for (size_t j = 0; j < i; ++j) {
        const double kMaxReprojError = 8.0;
        if (reconstructions[i]->Merge(*reconstructions[j], kMaxReprojError)) {
          reconstructions.erase(reconstructions.begin() + j);
          merge_success = true;
          break;
        }
      }

      if (merge_success) {
        break;
      }
    }

    if (!merge_success) {
      break;
    }
  }

  // Insert a new reconstruction manager for merged cluster.
  auto& reconstruction_manager = (*reconstruction_managers)[&cluster];
  for (const auto& reconstruction : reconstructions) {
    reconstruction_manager.Add();
    reconstruction_manager.Get(reconstruction_manager.Size() - 1) =
        *reconstruction;
  }

  // Delete all merged child cluster reconstruction managers.
  for (const auto& child_cluster : cluster.child_clusters) {
    reconstruction_managers->erase(&child_cluster);
  }
}

}  // namespace

bool HierarchicalMapperController::Options::Check() const {
  CHECK_OPTION_GT(init_num_trials, -1);
  CHECK_OPTION_GE(num_workers, -1);
  
  if (cluster_outpath.empty()) {
    std::cerr << "ERROR: cluster_outpath is empty" << std::endl;
    return false;
  }
  // Check custom cluster options
  if (use_custom_clusters) {
    if (custom_cluster_list_path.empty()) {
      std::cerr << "ERROR: Custom cluster list path is empty" << std::endl;
      return false;
    }
    std::ifstream test_file(custom_cluster_list_path);
    if (!test_file.is_open()) {
      std::cerr << "ERROR: Custom cluster list file does not exist or cannot be opened: " 
                << custom_cluster_list_path << std::endl;
      return false;
    }
    test_file.close();
  }
  
  return true;
}

HierarchicalMapperController::HierarchicalMapperController(
    const Options& options, const SceneClustering::Options& clustering_options,
    const IncrementalMapperOptions& mapper_options,
    ReconstructionManager* reconstruction_manager)
    : options_(options),
      clustering_options_(clustering_options),
      mapper_options_(mapper_options),
      reconstruction_manager_(reconstruction_manager) {
  CHECK(options_.Check());
  CHECK(clustering_options_.Check());
  CHECK(mapper_options_.Check());
  CHECK_EQ(clustering_options_.branching, 2);
}

void HierarchicalMapperController::Run() {
  PrintHeading1("Partitioning the scene");

  //////////////////////////////////////////////////////////////////////////////
  // Cluster scene
  //////////////////////////////////////////////////////////////////////////////

  std::unordered_map<image_t, std::string> image_id_to_name;
  std::unordered_map<std::string, image_t> image_name_to_id;

  Database database(options_.database_path);

  std::cout << "Reading images..." << std::endl;
  const auto images = database.ReadAllImages();
  for (const auto& image : images) {
    image_id_to_name.emplace(image.ImageId(), image.Name());
    image_name_to_id.emplace(image.Name(), image.ImageId());
  }


  //////////////////////////////////////////////////////////////////////////////
  // Create clusters (custom or automatic)
  //////////////////////////////////////////////////////////////////////////////
  
  std::vector<SceneClustering::Cluster> custom_leaf_clusters;
  std::vector<const SceneClustering::Cluster*> leaf_clusters;
  std::unique_ptr<SceneClustering> scene_clustering_ptr;
  
  std::cout << "Using cluster_outpath as : " << options_.cluster_outpath << std::endl;
  
  if (options_.use_custom_clusters) {
    std::cout << "Using custom clusters from: " << options_.custom_cluster_list_path << std::endl;
    
    // Read custom cluster files
    auto custom_clusters = ReadCustomClusterFiles(options_.custom_cluster_list_path);
    
    if (custom_clusters.empty()) {
      std::cerr << "ERROR: No valid clusters found in custom cluster files!" << std::endl;
      return;
    }
    
    // Convert custom clusters to SceneClustering::Cluster format
    custom_leaf_clusters.reserve(custom_clusters.size());
    for (size_t i = 0; i < custom_clusters.size(); ++i) {
      auto image_ids = ConvertNamesToIds(custom_clusters[i], image_name_to_id);
      if (!image_ids.empty()) {
        SceneClustering::Cluster cluster;
        cluster.image_ids = image_ids;
        custom_leaf_clusters.push_back(std::move(cluster));
      } else {
        std::cerr << "WARNING: Cluster " << (i + 1) << " has no valid images, skipping." << std::endl;
      }
    }
    
    // Create pointers to custom clusters
    for (auto& cluster : custom_leaf_clusters) {
      leaf_clusters.push_back(&cluster);
    }
    
    std::cout << "Successfully created " << leaf_clusters.size() << " custom clusters." << std::endl;
  } else {
    std::cout << "Using automatic clustering..." << std::endl;
    
    // Use original automatic clustering - create object with unique_ptr for proper memory management
    scene_clustering_ptr = std::make_unique<SceneClustering>(
        SceneClustering::Create(clustering_options_, database));
    leaf_clusters = scene_clustering_ptr->GetLeafClusters();
  }

  size_t total_num_images = 0;
  
  // Create clusters directory
  const std::string clusters_dir = options_.cluster_outpath;
  
  auto CreateDir = [](const std::string& path) {
#if defined(_WIN32)
    return _mkdir(path.c_str()) == 0;
#else 
    return mkdir(path.c_str(), 0733) == 0;
#endif
  };
  
  auto ExistsDir = [](const std::string& path) {
    struct stat info;
    return stat(path.c_str(), &info) == 0 && (info.st_mode & S_IFDIR);
  };
  
  if (!ExistsDir(clusters_dir)) {
    CreateDir(clusters_dir);
  }


  // Print cluster information with meaningful names
  std::cout << "    ######################################" << std::endl;  
  for (size_t i = 0; i < leaf_clusters.size(); ++i) {
    const auto& cluster = leaf_clusters[i];
    total_num_images += cluster->image_ids.size();


    std::cout << StringPrintf("  %s %d with %d images:", 
                              options_.use_custom_clusters ? "Custom Cluster" : "Cluster",
                              i + 1, static_cast<int>(cluster->image_ids.size()))
              << std::endl;


    for (const image_t image_id : cluster->image_ids) {
      std::cout << "    Image ID: " << image_id;


      auto it = image_id_to_name.find(image_id);
      if (it != image_id_to_name.end()) {
        std::cout << " (" << it->second << ")";
      }
      std::cout << std::endl;
    }
  }

  std::cout << StringPrintf("Clusters have %d images", static_cast<int>(total_num_images))
            << std::endl;

  //////////////////////////////////////////////////////////////////////////////
  // Reconstruct clusters
  //////////////////////////////////////////////////////////////////////////////

  PrintHeading1("Reconstructing clusters");

  // Determine the number of workers and threads per worker.
  const int kMaxNumThreads = -1;
  const int num_eff_threads = GetEffectiveNumThreads(kMaxNumThreads);
  const int kDefaultNumWorkers = 8;
  const int num_eff_workers =
      options_.num_workers < 1
          ? std::min(static_cast<int>(leaf_clusters.size()),
                     std::min(kDefaultNumWorkers, num_eff_threads))
          : options_.num_workers;
  const int num_threads_per_worker =
      std::max(1, num_eff_threads / num_eff_workers);

  // Function to reconstruct one cluster using incremental mapping.
  auto ReconstructCluster = [&, this](
                                const SceneClustering::Cluster& cluster,
                                ReconstructionManager* reconstruction_manager,
                                const std::string& cluster_name) {
    if (cluster.image_ids.empty()) {
      return;
    }

    IncrementalMapperOptions custom_options = mapper_options_;
    custom_options.max_model_overlap = 3;
    custom_options.init_num_trials = options_.init_num_trials;
    if (custom_options.num_threads < 0) {
      custom_options.num_threads = num_threads_per_worker;
    }

    for (const auto image_id : cluster.image_ids) {
      custom_options.image_names.insert(image_id_to_name.at(image_id));
    }


    // Create cluster-specific directory
    const std::string cluster_dir = clusters_dir + "/" + cluster_name;
    if (!ExistsDir(cluster_dir)) {
      CreateDir(cluster_dir);
    }


    // Create a temporary reconstruction manager for this cluster
    ReconstructionManager cluster_reconstruction_manager;
    
    IncrementalMapperController mapper(&custom_options, options_.image_path,
                                       options_.database_path,
                                       &cluster_reconstruction_manager);
    mapper.Start();
    mapper.Wait();


    // Save the cluster reconstruction
    if (cluster_reconstruction_manager.Size() > 0) {
      // Save the reconstruction in the cluster directory
      cluster_reconstruction_manager.Get(0).WriteBinary(cluster_dir);
      cluster_reconstruction_manager.Get(0).WriteText(cluster_dir);
      cluster_reconstruction_manager.Get(0).ExportPLY(cluster_dir + "/points.ply");
      
      // Save statistics
      std::ofstream stats_file(cluster_dir + "/cluster_stats.txt");
      if (stats_file.is_open()) {
        const Reconstruction& reconstruction = cluster_reconstruction_manager.Get(0);
        stats_file << "Cluster: " << cluster_name << "\n";
        stats_file << "Number of images: " << reconstruction.NumImages() << "\n";
        stats_file << "Registered images: " << reconstruction.NumRegImages() << "\n";
        stats_file << "Number of points: " << reconstruction.NumPoints3D() << "\n";
        stats_file << "Number of cameras: " << reconstruction.NumCameras() << "\n";
        stats_file.close();
      }
      
      LOG(INFO) << "Saved cluster reconstruction to: " << cluster_dir;
      
      // Copy the reconstruction to the main reconstruction manager
      reconstruction_manager->Add();
      reconstruction_manager->Get(reconstruction_manager->Size() - 1) = 
          cluster_reconstruction_manager.Get(0);
    }
  };

  // Start reconstructing the bigger clusters first for resource usage.
  std::sort(leaf_clusters.begin(), leaf_clusters.end(),
            [](const SceneClustering::Cluster* cluster1,
               const SceneClustering::Cluster* cluster2) {
              return cluster1->image_ids.size() > cluster2->image_ids.size();
            });

  // Start the reconstruction workers.

  std::unordered_map<const SceneClustering::Cluster*, ReconstructionManager>
      reconstruction_managers;
  reconstruction_managers.reserve(leaf_clusters.size());

  ThreadPool thread_pool(num_eff_workers);
  for (size_t i = 0; i < leaf_clusters.size(); ++i) {
    const auto& cluster = leaf_clusters[i];
    const std::string cluster_name = options_.use_custom_clusters ? 
                                     ("custom_cluster_" + std::to_string(i + 1)) :
                                     ("cluster_" + std::to_string(i + 1));
    
    thread_pool.AddTask(ReconstructCluster, *cluster,
                        &reconstruction_managers[cluster], cluster_name);
  }
  thread_pool.Wait();

  //////////////////////////////////////////////////////////////////////////////
  // Merge clusters
  //////////////////////////////////////////////////////////////////////////////

if (leaf_clusters.size() > 1) {
  PrintHeading1("Merging clusters");
  
  if (options_.use_custom_clusters) {
    std::cout << "Merging " << reconstruction_managers.size() << " custom clusters..." << std::endl;
    
    // Create a new reconstruction manager for the final result
    ReconstructionManager final_reconstruction_manager;
    
    // Find the largest reconstruction as the base
    ReconstructionManager* base_manager = nullptr;
    size_t max_points = 0;
    
    for (auto& pair : reconstruction_managers) {
      if (pair.second.Size() > 0) {
        size_t num_points = pair.second.Get(0).NumPoints3D();
        if (num_points > max_points) {
          max_points = num_points;
          base_manager = &pair.second;
        }
      }
    }
    
    if (base_manager && base_manager->Size() > 0) {
      // Copy the base reconstruction
      final_reconstruction_manager.Add();
      final_reconstruction_manager.Get(0) = base_manager->Get(0);
      
      std::cout << "Base reconstruction has " << max_points << " points" << std::endl;
      
      // Try to merge other reconstructions into the base
      for (auto& pair : reconstruction_managers) {
        if (&pair.second != base_manager && pair.second.Size() > 0) {
          const double kMaxReprojError = 8.0;
          size_t points_before = final_reconstruction_manager.Get(0).NumPoints3D();
          
          if (final_reconstruction_manager.Get(0).Merge(pair.second.Get(0), kMaxReprojError)) {
            size_t points_after = final_reconstruction_manager.Get(0).NumPoints3D();
            std::cout << "Successfully merged cluster, points: " << points_before 
                      << " -> " << points_after << std::endl;
          } else {
            std::cout << "Failed to merge cluster with " 
                      << pair.second.Get(0).NumPoints3D() << " points" << std::endl;
          }
        }
      }
      
      *reconstruction_manager_ = std::move(final_reconstruction_manager);
    }
    
  } else {
    // Use original merging for automatic clustering
    if (scene_clustering_ptr) {
      MergeClusters(*scene_clustering_ptr->GetRootCluster(), &reconstruction_managers);
    }
    
    CHECK_EQ(reconstruction_managers.size(), 1);
    *reconstruction_manager_ = std::move(reconstruction_managers.begin()->second);
  }
} else {
  // Only one cluster, no merging needed
  if (!reconstruction_managers.empty()) {
    *reconstruction_manager_ = std::move(reconstruction_managers.begin()->second);
  }
}
}
}  // namespace colmap
