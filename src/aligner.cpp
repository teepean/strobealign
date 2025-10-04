/*
 * Low-level alignment functions
 *
 * This is for anything that returns an aln_info object, currently
 * Aligner::align and hamming_align.
 */
#include <tuple>
#include <algorithm>
#include <cassert>
#include <ostream>
#include "aligner.hpp"

// Helper function to determine if a mismatch is a transition (C<->T or A<->G)
// Used for ancient DNA damage pattern support
static inline bool is_transition(char query_base, char ref_base) {
    // Transitions: C<->T and A<->G (purines and pyrimidines)
    if ((query_base == 'C' || query_base == 'c') && (ref_base == 'T' || ref_base == 't')) return true;
    if ((query_base == 'T' || query_base == 't') && (ref_base == 'C' || ref_base == 'c')) return true;
    if ((query_base == 'A' || query_base == 'a') && (ref_base == 'G' || ref_base == 'g')) return true;
    if ((query_base == 'G' || query_base == 'g') && (ref_base == 'A' || ref_base == 'a')) return true;
    return false;
}

// Get the appropriate mismatch penalty based on ancient DNA mode
static inline int get_mismatch_penalty(const AlignmentParameters& params, char query_base, char ref_base) {
    if (params.ancient_dna && is_transition(query_base, ref_base)) {
        return params.transition_penalty;
    }
    return params.ancient_dna ? params.transversion_penalty : params.mismatch;
}

std::optional<AlignmentInfo> Aligner::align(const std::string &query, const std::string &ref) const {
    m_align_calls++;
    AlignmentInfo aln;
    int32_t maskLen = query.length() / 2;
    maskLen = std::max(maskLen, 15);
    if (ref.length() > 2000){
//        std::cerr << "ALIGNMENT TO REF LONGER THAN 2000bp - REPORT TO DEVELOPER. Happened for read: " <<  query << " ref len:" << ref.length() << std::endl;
        return {};
    }

    StripedSmithWaterman::Alignment alignment_ssw;

    // query must be NULL-terminated
    auto flag = ssw_aligner.Align(query.c_str(), ref.c_str(), ref.size(), filter, &alignment_ssw, maskLen);
    if (flag != 0 || alignment_ssw.ref_begin == -1) {
        return {};
    }

    aln.edit_distance = alignment_ssw.mismatches;
    aln.cigar = Cigar(alignment_ssw.cigar);
    aln.sw_score = alignment_ssw.sw_score;
    aln.ref_start = alignment_ssw.ref_begin;
    // end positions are off by 1 in SSW
    aln.ref_end = alignment_ssw.ref_end + 1;
    aln.query_start = alignment_ssw.query_begin;
    aln.query_end = alignment_ssw.query_end + 1;


    // Try to extend to beginning of the query to get an end bonus
    auto qstart = aln.query_start;
    auto rstart = aln.ref_start;
    auto score = aln.sw_score;
    auto edits = aln.edit_distance;
    Cigar front_cigar;
    while (qstart > 0 && rstart > 0) {
        qstart--;
        rstart--;
        if (query[qstart] == ref[rstart]) {
            score += parameters.match;
            front_cigar.push(CIGAR_EQ, 1);
        } else {
            score -= get_mismatch_penalty(parameters, query[qstart], ref[rstart]);
            front_cigar.push(CIGAR_X, 1);
            edits++;
        }
    }
    if (qstart == 0 && score + parameters.end_bonus > aln.sw_score) {
        if (aln.query_start > 0) {
            assert((aln.cigar.m_ops[0] & 0xF) == CIGAR_SOFTCLIP);
            aln.cigar.m_ops.erase(aln.cigar.m_ops.begin());  // remove soft clipping
            front_cigar.reverse();
            front_cigar += aln.cigar;
            aln.cigar = std::move(front_cigar);
        }
        aln.query_start = 0;
        aln.ref_start = rstart;
        aln.sw_score = score + parameters.end_bonus;
        aln.edit_distance = edits;
    }

    // Try to extend to end of query to get an end bonus
    auto qend = aln.query_end;
    auto rend = aln.ref_end;
    score = aln.sw_score;
    edits = aln.edit_distance;
    Cigar back_cigar;
    while (qend < query.length() && rend < ref.length()) {
        if (query[qend] == ref[rend]) {
            score += parameters.match;
            back_cigar.push(CIGAR_EQ, 1);
        } else {
            score -= get_mismatch_penalty(parameters, query[qend], ref[rend]);
            back_cigar.push(CIGAR_X, 1);
            edits++;
        }
        qend++;
        rend++;
    }
    if (qend == query.length() && score + parameters.end_bonus > aln.sw_score) {
        if (aln.query_end < query.length()) {
            assert((aln.cigar.m_ops[aln.cigar.m_ops.size() - 1] & 0xf) == CIGAR_SOFTCLIP);
            aln.cigar.m_ops.pop_back();
            aln.cigar += back_cigar;
        }
        aln.query_end = query.length();
        aln.ref_end = rend;
        aln.sw_score = score + parameters.end_bonus;
        aln.edit_distance = edits;
    }

    // Final validation: ensure CIGAR covers exactly the query length
    // Calculate actual CIGAR length (operations that consume query bases: M, I, S, =, X)
    size_t cigar_query_len = 0;
    for (auto op : aln.cigar.m_ops) {
        int op_type = op & 0xF;
        int op_len = op >> 4;
        if (op_type == 0 || op_type == 1 || op_type == 4 || op_type == 7 || op_type == 8) { // M, I, S, =, X
            cigar_query_len += op_len;
        }
    }

    // Fix CIGAR length mismatch
    if (cigar_query_len < query.length()) {
        // CIGAR is too short - add soft clip at the end
        size_t missing = query.length() - cigar_query_len;
        if (!aln.cigar.m_ops.empty() && (aln.cigar.m_ops.back() & 0xF) == CIGAR_SOFTCLIP) {
            // Extend existing soft clip
            size_t existing_clip = aln.cigar.m_ops.back() >> 4;
            aln.cigar.m_ops.back() = ((existing_clip + missing) << 4) | CIGAR_SOFTCLIP;
        } else {
            // Add new soft clip
            aln.cigar.m_ops.push_back((missing << 4) | CIGAR_SOFTCLIP);
        }
    } else if (cigar_query_len > query.length()) {
        // CIGAR is too long - remove from soft clips
        size_t excess = cigar_query_len - query.length();

        // Try to remove from trailing soft clip first
        if (!aln.cigar.m_ops.empty() && (aln.cigar.m_ops.back() & 0xF) == CIGAR_SOFTCLIP) {
            size_t soft_clip_len = aln.cigar.m_ops.back() >> 4;
            if (soft_clip_len > excess) {
                aln.cigar.m_ops.back() = ((soft_clip_len - excess) << 4) | CIGAR_SOFTCLIP;
                excess = 0;
            } else {
                aln.cigar.m_ops.pop_back();
                excess -= soft_clip_len;
            }
        }

        // If still excess, try leading soft clip
        if (excess > 0 && !aln.cigar.m_ops.empty() && (aln.cigar.m_ops[0] & 0xF) == CIGAR_SOFTCLIP) {
            size_t soft_clip_len = aln.cigar.m_ops[0] >> 4;
            if (soft_clip_len > excess) {
                aln.cigar.m_ops[0] = ((soft_clip_len - excess) << 4) | CIGAR_SOFTCLIP;
            } else {
                aln.cigar.m_ops.erase(aln.cigar.m_ops.begin());
            }
        }
    }

    return aln;
}

/*
 * Find highest-scoring segment between reference and query assuming only matches
 * and mismatches are allowed.
 *
 * The end_bonus is added to the score if the segment extends until the end
 * of the query, once for each end.
 */
std::tuple<size_t, size_t, int> highest_scoring_segment(
    const std::string& query, const std::string& ref, const AlignmentParameters& params
) {
    size_t n = query.length();

    size_t start = 0; // start of the current segment
    int score = params.end_bonus; // accumulated score so far in the current segment

    size_t best_start = 0;
    size_t best_end = 0;
    int best_score = 0;
    for (size_t i = 0; i < n; ++i) {
        if (query[i] == ref[i]) {
            score += params.match;
        } else {
            score -= get_mismatch_penalty(params, query[i], ref[i]);
        }
        if (score < 0) {
            start = i + 1;
            score = 0;
        }
        if (score > best_score) {
            best_start = start;
            best_score = score;
            best_end = i + 1;
        }
    }
    if (score + params.end_bonus > best_score) {
        best_score = score + params.end_bonus;
        best_end = query.length();
        best_start = start;
    }
    return std::make_tuple(best_start, best_end, best_score);
}

AlignmentInfo hamming_align(
    const std::string &query, const std::string &ref, const AlignmentParameters& params
) {
    AlignmentInfo aln;
    if (query.length() != ref.length()) {
        return aln;
    }

    auto [segment_start, segment_end, score] = highest_scoring_segment(query, ref, params);

    Cigar cigar;
    if (segment_start > 0) {
        cigar.push(CIGAR_SOFTCLIP, segment_start);
    }

    // Create CIGAR string and count mismatches
    int counter = 0;
    bool prev_is_match = false;
    int mismatches = 0;
    bool first = true;
    for (size_t i = segment_start; i < segment_end; i++) {
        bool is_match = query[i] == ref[i];
        mismatches += is_match ? 0 : 1;
        if (!first && is_match != prev_is_match) {
            cigar.push(prev_is_match ? CIGAR_EQ : CIGAR_X, counter);
            counter = 0;
        }
        counter++;
        prev_is_match = is_match;
        first = false;
    }
    if (!first) {
        cigar.push(prev_is_match ? CIGAR_EQ : CIGAR_X, counter);
    }

    int soft_right = query.length() - segment_end;
    if (soft_right > 0) {
        cigar.push(CIGAR_SOFTCLIP, soft_right);
    }

    aln.cigar = std::move(cigar);
    aln.sw_score = score;
    aln.edit_distance = mismatches;
    aln.ref_start = segment_start;
    aln.ref_end = segment_end;
    aln.query_start = segment_start;
    aln.query_end = segment_end;
    return aln;
}

std::ostream& operator<<(std::ostream& os, const AlignmentParameters& params) {
    os
        << "AlignmentParameters("
        << "match=" << params.match
        << ", mismatch=" << params.mismatch
        << ", gap_open=" << params.gap_open
        << ", gap_extend=" << params.gap_extend
        << ", end_bonus=" << params.end_bonus;
    if (params.ancient_dna) {
        os << ", ancient_dna=true"
           << ", transition_penalty=" << params.transition_penalty
           << ", transversion_penalty=" << params.transversion_penalty;
    }
    os << ")";
    return os;
}
