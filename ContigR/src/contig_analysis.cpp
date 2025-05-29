#include <Rcpp.h>
#include "contigs.hpp"
#include "overlaps.hpp"
#include "assign_multiplicity.hpp"
#include "output.hpp"

using namespace Rcpp;

// [[Rcpp::export]]
List analyze_contigs_cpp(std::string filepath, int maxk, int mink, 
                        std::string output_dir, int num_iterations) {
    // Create output directory if it doesn't exist
    std::filesystem::path output_path(output_dir);
    if (!std::filesystem::exists(output_path)) {
        std::filesystem::create_directory(output_path);
    }
    
    // Create CSV file for execution times
    std::ofstream csv_file(output_path / "execution_times.csv");
    csv_file << "Iteration,Contig Collection (ms),Overlap Detection (ms),"
             << "Multiplicity Assignment (ms),File Writing (ms),Total Time (ms)\n";
    
    std::vector<ExecutionTimes> all_execution_times;
    
    for(int iteration = 0; iteration < num_iterations; ++iteration) {
        Rcout << "\nStarting iteration " << iteration + 1 << " of " 
              << num_iterations << std::endl;
        
        ExecutionTimes times;
        auto total_start_time = std::chrono::high_resolution_clock::now();
        
        // Get contig collection
        auto start_time = std::chrono::high_resolution_clock::now();
        ContigCollection contig_collection = get_contig_collection(filepath, maxk);
        auto end_time = std::chrono::high_resolution_clock::now();
        times.contig_collection_time = std::chrono::duration_cast<std::chrono::milliseconds>
            (end_time - start_time).count();
        
        // Detect overlaps
        start_time = std::chrono::high_resolution_clock::now();
        OverlapCollection overlap_collection = detect_adjacent_contigs(contig_collection, mink, maxk);
        end_time = std::chrono::high_resolution_clock::now();
        times.overlap_detection_time = std::chrono::duration_cast<std::chrono::milliseconds>
            (end_time - start_time).count();
        
        // Assign multiplicity
        start_time = std::chrono::high_resolution_clock::now();
        assign_multiplicity(contig_collection, overlap_collection);
        end_time = std::chrono::high_resolution_clock::now();
        times.multiplicity_assignment_time = std::chrono::duration_cast<std::chrono::milliseconds>
            (end_time - start_time).count();
        
        // Write output files
        start_time = std::chrono::high_resolution_clock::now();
        std::filesystem::path iteration_dir = output_path / 
            ("iteration_" + std::to_string(iteration + 1));
        if (!std::filesystem::exists(iteration_dir)) {
            std::filesystem::create_directory(iteration_dir);
        }
        std::string outdpath = iteration_dir.string();
        
        write_summary(contig_collection, overlap_collection, filepath, outdpath);
        write_adjacency_table_and_full_log(contig_collection, overlap_collection, outdpath);
        write_genbank(contig_collection, overlap_collection, outdpath);
        end_time = std::chrono::high_resolution_clock::now();
        times.file_writing_time = std::chrono::duration_cast<std::chrono::milliseconds>
            (end_time - start_time).count();
        
        // Calculate total time
        auto total_end_time = std::chrono::high_resolution_clock::now();
        times.total_time = std::chrono::duration_cast<std::chrono::milliseconds>
            (total_end_time - total_start_time).count();
        
        all_execution_times.push_back(times);
        
        // Write to CSV
        csv_file << iteration + 1 << ","
                << times.contig_collection_time << ","
                << times.overlap_detection_time << ","
                << times.multiplicity_assignment_time << ","
                << times.file_writing_time << ","
                << times.total_time << "\n";
        
        Rcout << "Iteration " << iteration + 1 << " execution times:" << std::endl;
        Rcout << "  Contig Collection: " << times.contig_collection_time << " ms" << std::endl;
        Rcout << "  Overlap Detection: " << times.overlap_detection_time << " ms" << std::endl;
        Rcout << "  Multiplicity Assignment: " << times.multiplicity_assignment_time << " ms" << std::endl;
        Rcout << "  File Writing: " << times.file_writing_time << " ms" << std::endl;
        Rcout << "  Total Time: " << times.total_time << " ms" << std::endl;
    }
    
    // Calculate averages
    ExecutionTimes avg_times = {0, 0, 0, 0, 0};
    for(const auto& times : all_execution_times) {
        avg_times.contig_collection_time += times.contig_collection_time;
        avg_times.overlap_detection_time += times.overlap_detection_time;
        avg_times.multiplicity_assignment_time += times.multiplicity_assignment_time;
        avg_times.file_writing_time += times.file_writing_time;
        avg_times.total_time += times.total_time;
    }
    
    avg_times.contig_collection_time /= num_iterations;
    avg_times.overlap_detection_time /= num_iterations;
    avg_times.multiplicity_assignment_time /= num_iterations;
    avg_times.file_writing_time /= num_iterations;
    avg_times.total_time /= num_iterations;
    
    // Write averages to CSV
    csv_file << "\nAverage,"
            << avg_times.contig_collection_time << ","
            << avg_times.overlap_detection_time << ","
            << avg_times.multiplicity_assignment_time << ","
            << avg_times.file_writing_time << ","
            << avg_times.total_time << "\n";
    
    csv_file.close();
    
    // Write final statistics
    std::ofstream stats_file(output_path / "final_statistics.txt");
    stats_file << "Final Statistics\n";
    stats_file << "================\n\n";
    stats_file << "Average execution times:\n";
    stats_file << "  Contig Collection: " << avg_times.contig_collection_time << " ms\n";
    stats_file << "  Overlap Detection: " << avg_times.overlap_detection_time << " ms\n";
    stats_file << "  Multiplicity Assignment: " << avg_times.multiplicity_assignment_time << " ms\n";
    stats_file << "  File Writing: " << avg_times.file_writing_time << " ms\n";
    stats_file << "  Total Time: " << avg_times.total_time << " ms\n";
    stats_file.close();
    
    // Return results as a list
    return List::create(
        Named("adjacency_table_path") = (output_path / "iteration_1__adjacent_contigs.tsv").string(),
        Named("execution_times") = DataFrame::create(
            Named("iteration") = seq_len(num_iterations),
            Named("contig_collection") = sapply(all_execution_times, [](const ExecutionTimes& t) { return t.contig_collection_time; }),
            Named("overlap_detection") = sapply(all_execution_times, [](const ExecutionTimes& t) { return t.overlap_detection_time; }),
            Named("multiplicity_assignment") = sapply(all_execution_times, [](const ExecutionTimes& t) { return t.multiplicity_assignment_time; }),
            Named("file_writing") = sapply(all_execution_times, [](const ExecutionTimes& t) { return t.file_writing_time; }),
            Named("total_time") = sapply(all_execution_times, [](const ExecutionTimes& t) { return t.total_time; })
        ),
        Named("average_times") = List::create(
            Named("contig_collection") = avg_times.contig_collection_time,
            Named("overlap_detection") = avg_times.overlap_detection_time,
            Named("multiplicity_assignment") = avg_times.multiplicity_assignment_time,
            Named("file_writing") = avg_times.file_writing_time,
            Named("total_time") = avg_times.total_time
        )
    );
} 