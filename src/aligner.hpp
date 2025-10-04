#ifndef STROBEALIGN_ALIGNER_HPP
#define STROBEALIGN_ALIGNER_HPP

#include <string>
#include <tuple>
#include <optional>
#include "ssw/ssw_cpp.h"
#include "cigar.hpp"


struct AlignmentParameters {
    // match is a score, the others are penalties (all are nonnegative)
    int match;
    int mismatch;
    int gap_open;
    int gap_extend;
    int end_bonus;

    // Ancient DNA support: separate penalties for transitions vs transversions
    bool ancient_dna;
    int transition_penalty;    // C<->T, A<->G
    int transversion_penalty;  // All other mismatches
};

std::ostream& operator<<(std::ostream& os, const AlignmentParameters& params);

struct AlignmentInfo {
    Cigar cigar;
    unsigned int edit_distance{0};
    unsigned int ref_start{0};
    unsigned int ref_end{0};
    unsigned int query_start{0};
    unsigned int query_end{0};
    int sw_score{0};

    int ref_span() const { return ref_end - ref_start; }
};

struct Aligner {
public:
    Aligner(AlignmentParameters parameters)
        : parameters(parameters)
        , ssw_aligner(StripedSmithWaterman::Aligner(
            parameters.match,
            // Use average penalty for ancient DNA mode since SSW doesn't support per-base penalties
            // Assuming ~50% transitions in aDNA: (transition_penalty + transversion_penalty) / 2
            parameters.ancient_dna ? (parameters.transition_penalty + parameters.transversion_penalty) / 2 : parameters.mismatch,
            parameters.gap_open,
            parameters.gap_extend))
    { }

    std::optional<AlignmentInfo> align(const std::string &query, const std::string &ref) const;

    AlignmentParameters parameters;

    unsigned calls_count() {
        return m_align_calls;
    }

private:
    const StripedSmithWaterman::Aligner ssw_aligner;
    const StripedSmithWaterman::Filter filter;
    mutable unsigned m_align_calls{0};  // no. of calls to the align() method
};

inline int hamming_distance(const std::string &s, const std::string &t) {
    if (s.length() != t.length()){
        return -1;
    }

    int mismatches = 0;
    for (size_t i = 0; i < s.length(); i++) {
        if (s[i] != t[i]) {
            mismatches++;
        }
    }

    return mismatches;
}

std::tuple<size_t, size_t, int> highest_scoring_segment(
    const std::string& query, const std::string& ref, const AlignmentParameters& params
);

AlignmentInfo hamming_align(
    const std::string &query, const std::string &ref, const AlignmentParameters& params
);

#endif
