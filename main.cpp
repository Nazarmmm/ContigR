#include <iostream> //fasta bondage
#include <string>
#include <chrono>
#include <vector>
#include <fstream>
#include <iomanip>
#include <filesystem>

#include "contigs.hpp"                                //2 files cra, c,сборщик spades оставляет концы равное к-меру, fastg bondage
#include "overlaps.hpp"
#include "assign_multiplicity.hpp"
#include "output.hpp"

using namespace std;
namespace fs = std::filesystem;

// Структура для хранения времени выполнения каждого этапа
struct ExecutionTimes {
    long long contig_collection_time;
    long long overlap_detection_time;
    long long multiplicity_assignment_time;
    long long file_writing_time;
    long long total_time;
};

int main() {
    const int NUM_ITERATIONS = 100;
    std::vector<ExecutionTimes> all_execution_times;
    
    std::string filepath = "C:/Users/admin/Downloads/testg.fasta";
    int maxk=50;
    int mink=5;
    
    // Создаем директорию Output, если она не существует
    fs::path output_dir = "Output";
    if (!fs::exists(output_dir)) {
        fs::create_directory(output_dir);
    }
    
    // Создаем CSV файл для записи результатов в папке Output
    std::ofstream csv_file(output_dir / "execution_times.csv");
    csv_file << "Iteration,Contig Collection (ms),Overlap Detection (ms),Multiplicity Assignment (ms),File Writing (ms),Total Time (ms)\n";
    
    for(int iteration = 0; iteration < NUM_ITERATIONS; ++iteration) {
        std::cout << "\nStarting iteration " << iteration + 1 << " of " << NUM_ITERATIONS << std::endl;
        
        ExecutionTimes times;
        auto total_start_time = std::chrono::high_resolution_clock::now();
        
        // Замер времени получения коллекции контигов
        auto start_time = std::chrono::high_resolution_clock::now();
        ContigCollection contig_collection = get_contig_collection(filepath, maxk);
        auto end_time = std::chrono::high_resolution_clock::now();
        times.contig_collection_time = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time).count();
        
        // Замер времени определения перекрытий
        start_time = std::chrono::high_resolution_clock::now();
        OverlapCollection overlap_collection = detect_adjacent_contigs(contig_collection, mink, maxk);
        end_time = std::chrono::high_resolution_clock::now();
        times.overlap_detection_time = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time).count();
        
        // Замер времени назначения множественности
        start_time = std::chrono::high_resolution_clock::now();
        assign_multiplicity(contig_collection, overlap_collection);
        end_time = std::chrono::high_resolution_clock::now();
        times.multiplicity_assignment_time = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time).count();
        
        // Замер времени записи файлов
        start_time = std::chrono::high_resolution_clock::now();
        // Создаем поддиректорию для текущей итерации
        fs::path iteration_dir = output_dir / ("iteration_" + std::to_string(iteration + 1));
        if (!fs::exists(iteration_dir)) {
            fs::create_directory(iteration_dir);
        }
        std::string outdpath = iteration_dir.string();
        
        write_summary(contig_collection, overlap_collection, filepath, outdpath);
        write_adjacency_table_and_full_log(contig_collection, overlap_collection, outdpath);
        write_genbank(contig_collection, overlap_collection, outdpath);
        end_time = std::chrono::high_resolution_clock::now();
        times.file_writing_time = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time).count();
        
        // Подсчет общего времени
        auto total_end_time = std::chrono::high_resolution_clock::now();
        times.total_time = std::chrono::duration_cast<std::chrono::milliseconds>(total_end_time - total_start_time).count();
        
        all_execution_times.push_back(times);
        
        // Запись результатов в CSV
        csv_file << iteration + 1 << ","
                << times.contig_collection_time << ","
                << times.overlap_detection_time << ","
                << times.multiplicity_assignment_time << ","
                << times.file_writing_time << ","
                << times.total_time << "\n";
        
        std::cout << "Iteration " << iteration + 1 << " execution times:" << std::endl;
        std::cout << "  Contig Collection: " << times.contig_collection_time << " ms" << std::endl;
        std::cout << "  Overlap Detection: " << times.overlap_detection_time << " ms" << std::endl;
        std::cout << "  Multiplicity Assignment: " << times.multiplicity_assignment_time << " ms" << std::endl;
        std::cout << "  File Writing: " << times.file_writing_time << " ms" << std::endl;
        std::cout << "  Total Time: " << times.total_time << " ms" << std::endl;
    }
    
    // Вычисляем и записываем средние значения
    ExecutionTimes avg_times = {0, 0, 0, 0, 0};
    for(const auto& times : all_execution_times) {
        avg_times.contig_collection_time += times.contig_collection_time;
        avg_times.overlap_detection_time += times.overlap_detection_time;
        avg_times.multiplicity_assignment_time += times.multiplicity_assignment_time;
        avg_times.file_writing_time += times.file_writing_time;
        avg_times.total_time += times.total_time;
    }
    
    avg_times.contig_collection_time /= NUM_ITERATIONS;
    avg_times.overlap_detection_time /= NUM_ITERATIONS;
    avg_times.multiplicity_assignment_time /= NUM_ITERATIONS;
    avg_times.file_writing_time /= NUM_ITERATIONS;
    avg_times.total_time /= NUM_ITERATIONS;
    
    csv_file << "\nAverage,"
            << avg_times.contig_collection_time << ","
            << avg_times.overlap_detection_time << ","
            << avg_times.multiplicity_assignment_time << ","
            << avg_times.file_writing_time << ","
            << avg_times.total_time << "\n";
    
    csv_file.close();
    
    // Записываем итоговую статистику в отдельный файл
    std::ofstream stats_file(output_dir / "final_statistics.txt");
    stats_file << "Final Statistics\n";
    stats_file << "================\n\n";
    stats_file << "Average execution times:\n";
    stats_file << "  Contig Collection: " << avg_times.contig_collection_time << " ms\n";
    stats_file << "  Overlap Detection: " << avg_times.overlap_detection_time << " ms\n";
    stats_file << "  Multiplicity Assignment: " << avg_times.multiplicity_assignment_time << " ms\n";
    stats_file << "  File Writing: " << avg_times.file_writing_time << " ms\n";
    stats_file << "  Total Time: " << avg_times.total_time << " ms\n";
    stats_file.close();
    
    std::cout << "\nAverage execution times:" << std::endl;
    std::cout << "  Contig Collection: " << avg_times.contig_collection_time << " ms" << std::endl;
    std::cout << "  Overlap Detection: " << avg_times.overlap_detection_time << " ms" << std::endl;
    std::cout << "  Multiplicity Assignment: " << avg_times.multiplicity_assignment_time << " ms" << std::endl;
    std::cout << "  File Writing: " << avg_times.file_writing_time << " ms" << std::endl;
    std::cout << "  Total Time: " << avg_times.total_time << " ms" << std::endl;
    
    std::cout << "\nResults have been saved to the Output directory" << std::endl;
    
    return 0;
}