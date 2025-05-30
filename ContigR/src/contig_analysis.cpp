#include <Rcpp.h>
#include <filesystem>
#include <fstream>
#include <vector>
#include <chrono>
#include "contigs.hpp"
#include "overlaps.hpp"
#include "assign_multiplicity.hpp"
#include "output.hpp"

using namespace Rcpp;

// [[Rcpp::export]]
List analyze_contigs_cpp(std::string filepath, int maxk, int mink, 
                         std::string output_dir, int num_iterations) {

    std::filesystem::path output_path(output_dir);
    if (!std::filesystem::exists(output_path)) {
        std::filesystem::create_directory(output_path);
    }

    std::ofstream csv_file(output_path / "execution_times.csv");
    csv_file << "Iteration,Contig Collection (ms),Overlap Detection (ms),"
             << "Multiplicity Assignment (ms),File Writing (ms),Total Time (ms)\n";

    std::vector<long> contig_collection_times;
    std::vector<long> overlap_detection_times;
    std::vector<long> multiplicity_assignment_times;
    std::vector<long> file_writing_times;
    std::vector<long> total_times;

    for (int iteration = 0; iteration < num_iterations; ++iteration) {
        Rcout << "\nStarting iteration " << iteration + 1 << " of " << num_iterations << std::endl;

        auto total_start_time = std::chrono::high_resolution_clock::now();

        // Contig Collection
        auto start_time = std::chrono::high_resolution_clock::now();
        ContigCollection contig_collection = get_contig_collection(filepath, maxk);
        auto end_time = std::chrono::high_resolution_clock::now();
        long contig_time = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time).count();

        // Overlap Detection
        start_time = std::chrono::high_resolution_clock::now();
        OverlapCollection overlap_collection = detect_adjacent_contigs(contig_collection, mink, maxk);
        end_time = std::chrono::high_resolution_clock::now();
        long overlap_time = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time).count();

        // Multiplicity Assignment
        start_time = std::chrono::high_resolution_clock::now();
        assign_multiplicity(contig_collection, overlap_collection);
        end_time = std::chrono::high_resolution_clock::now();
        long multiplicity_time = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time).count();

        // File Writing
        start_time = std::chrono::high_resolution_clock::now();
        std::filesystem::path iteration_dir = output_path / ("iteration_" + std::to_string(iteration + 1));
        if (!std::filesystem::exists(iteration_dir)) {
            std::filesystem::create_directory(iteration_dir);
        }
        std::string outdpath = iteration_dir.string();

        write_summary(contig_collection, overlap_collection, filepath, outdpath);
        write_adjacency_table_and_full_log(contig_collection, overlap_collection, outdpath);
        write_genbank(contig_collection, overlap_collection, outdpath);
        end_time = std::chrono::high_resolution_clock::now();
        long file_time = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time).count();

        // Total time
        auto total_end_time = std::chrono::high_resolution_clock::now();
        long total_time = std::chrono::duration_cast<std::chrono::milliseconds>(total_end_time - total_start_time).count();

        contig_collection_times.push_back(contig_time);
        overlap_detection_times.push_back(overlap_time);
        multiplicity_assignment_times.push_back(multiplicity_time);
        file_writing_times.push_back(file_time);
        total_times.push_back(total_time);

        csv_file << iteration + 1 << "," << contig_time << "," << overlap_time << ","
                 << multiplicity_time << "," << file_time << "," << total_time << "\n";

        Rcout << "Iteration " << iteration + 1 << " execution times:\n"
              << "  Contig Collection: " << contig_time << " ms\n"
              << "  Overlap Detection: " << overlap_time << " ms\n"
              << "  Multiplicity Assignment: " << multiplicity_time << " ms\n"
              << "  File Writing: " << file_time << " ms\n"
              << "  Total Time: " << total_time << " ms\n";
    }

    // Calculate averages
    long avg_contig = std::accumulate(contig_collection_times.begin(), contig_collection_times.end(), 0L) / num_iterations;
    long avg_overlap = std::accumulate(overlap_detection_times.begin(), overlap_detection_times.end(), 0L) / num_iterations;
    long avg_multiplicity = std::accumulate(multiplicity_assignment_times.begin(), multiplicity_assignment_times.end(), 0L) / num_iterations;
    long avg_file = std::accumulate(file_writing_times.begin(), file_writing_times.end(), 0L) / num_iterations;
    long avg_total = std::accumulate(total_times.begin(), total_times.end(), 0L) / num_iterations;

    csv_file << "\nAverage," << avg_contig << "," << avg_overlap << ","
             << avg_multiplicity << "," << avg_file << "," << avg_total << "\n";

    csv_file.close();

    std::ofstream stats_file(output_path / "final_statistics.txt");
    stats_file << "Final Statistics\n================\n\n"
               << "Average execution times:\n"
               << "  Contig Collection: " << avg_contig << " ms\n"
               << "  Overlap Detection: " << avg_overlap << " ms\n"
               << "  Multiplicity Assignment: " << avg_multiplicity << " ms\n"
               << "  File Writing: " << avg_file << " ms\n"
               << "  Total Time: " << avg_total << " ms\n";
    stats_file.close();

    return List::create(
        Named("adjacency_table_path") = (output_path / "iteration_1__adjacent_contigs.tsv").string(),
        Named("execution_times") = DataFrame::create(
            Named("iteration") = seq_len(num_iterations),
            Named("contig_collection") = contig_collection_times,
            Named("overlap_detection") = overlap_detection_times,
            Named("multiplicity_assignment") = multiplicity_assignment_times,
            Named("file_writing") = file_writing_times,
            Named("total_time") = total_times
        ),
        Named("average_times") = List::create(
            Named("contig_collection") = avg_contig,
            Named("overlap_detection") = avg_overlap,
            Named("multiplicity_assignment") = avg_multiplicity,
            Named("file_writing") = avg_file,
            Named("total_time") = avg_total
        )
    );
}
