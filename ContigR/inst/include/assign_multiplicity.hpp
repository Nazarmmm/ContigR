#pragma once

#include <algorithm>

#include "contigs.hpp"
#include "overlaps.hpp"
#include "output.hpp"

using namespace std;

float _calc_multiplty_by_coverage(float curr_coverage, float first_contig_coverage) {
    return curr_coverage / first_contig_coverage;
}

float _calc_multiply_by_overlaps(const std::vector<Overlap>& ovl_list) {
    // Function for calculating multiplicity of a given contig
    // based on the number of overlaps of this contig.
    // :param ovl_list: list of overlaps of current contig;

    // Count overlaps associated with start
    int num_start_matches = std::count_if(ovl_list.begin(), ovl_list.end(), is_start_match);
    // Count overlaps associated with end
    int num_end_matches = std::count_if(ovl_list.begin(), ovl_list.end(), is_end_match);

    // Obtain multiplicity based on the number of overlaps
    float multiplicity = std::max(1, std::min(num_start_matches, num_end_matches));

    return multiplicity;
}

void assign_multiplicity(ContigCollection& contig_collection, const OverlapCollection& overlap_collection) {
    if (contig_collection.empty()) {
        std::cerr << "[ERROR] Contig collection is empty. Cannot assign multiplicity." << std::endl;
        return;
    }
    // Function assigns multiplicity (copies of this contig in the genome) to contigs.
    // :param contig_collection: instance of `ContigCollection`
    //   returned by function `get_contig_collection`;
    // Function modifies `contig_collection` argument passed to it.

    const float COVERAGE_THRESHOLD = 1e-6;
    const float first_contig_coverage = contig_collection[0].cov;
    const bool first_cov_is_valid = first_contig_coverage > COVERAGE_THRESHOLD;

    if (!first_cov_is_valid) {
        std::cout << "\nFirst contig has insufficient coverage (less than " << COVERAGE_THRESHOLD << ").\n"
                  << "Multiplicity of contigs will be calculated based on overlaps instead of coverage.\n" << std::endl;
    }

    // Calculate multiplicity of contigs:
    #pragma omp parallel for
    for (size_t i = 0; i < contig_collection.size(); ++i) {
        if (first_cov_is_valid && contig_collection[i].cov > COVERAGE_THRESHOLD) {
            contig_collection[i].multplty = _calc_multiplty_by_coverage(
                contig_collection[i].cov,
                first_contig_coverage
            );
        } else {
            contig_collection[i].multplty = _calc_multiply_by_overlaps(
                overlap_collection[i]
            );
        }
    }
} 