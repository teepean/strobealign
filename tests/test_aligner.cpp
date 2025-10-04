#include "doctest.h"
#include "aligner.hpp"

TEST_CASE("hamming_align") {
    // empty sequences
    AlignmentParameters params{7, 5, 0, 0, 0, false, 0, 0};
    auto info = hamming_align(
        "", "",
        params
    );
    CHECK(info.cigar.empty());
    CHECK(info.edit_distance == 0);
    CHECK(info.sw_score == 0);
    CHECK(info.ref_start == 0);
    CHECK(info.ref_span() == 0);
    CHECK(info.query_start == 0);
    CHECK(info.query_end == 0);

    params = {1, 1, 0, 0, 0, false, 0, 0};
    info = hamming_align("AAXGGG", "AAYGGG", params);
    CHECK(info.cigar.to_string() == "2=1X3=");
    CHECK(info.edit_distance == 1);
    CHECK(info.sw_score == 4);
    CHECK(info.ref_start == 0);
    CHECK(info.ref_span() == 6);
    CHECK(info.query_start == 0);
    CHECK(info.query_end == 6);

    params = {1, 4, 0, 0, 0, false, 0, 0};
    info = hamming_align("AXGGG", "AYGGG", params);
    CHECK(info.cigar.to_string() == "2S3=");
    CHECK(info.edit_distance == 0);
    CHECK(info.sw_score == 3);
    CHECK(info.ref_start == 2);
    CHECK(info.ref_span() == 3);
    CHECK(info.query_start == 2);
    CHECK(info.query_end == 5);

    params = {3, 7, 0, 0, 0, false, 0, 0};
    info = hamming_align("NAACCG", "TAACCG", params);
    CHECK(info.cigar.to_string() == "1S5=");
    CHECK(info.edit_distance == 0);
    CHECK(info.sw_score == 5 * 3);
    CHECK(info.ref_start == 1);
    CHECK(info.ref_span() == 5);
    CHECK(info.query_start == 1);
    CHECK(info.query_end == 6);

    info = hamming_align("AACCGN", "AACCGT", params);
    CHECK(info.cigar.to_string() == "5=1S");
    CHECK(info.edit_distance == 0);
    CHECK(info.sw_score == 5 * 3);
    CHECK(info.ref_start == 0);
    CHECK(info.ref_span() == 5);
    CHECK(info.query_start == 0);
    CHECK(info.query_end == 5);

    // negative total score, soft clipping on both ends
    info = hamming_align("NAAAAAAAAAAAAAA", "TAAAATTTTTTTTTT", params);
    CHECK(info.cigar.to_string() == "1S4=10S");
    CHECK(info.edit_distance == 0);
    CHECK(info.sw_score == 4 * 3);
    CHECK(info.ref_start == 1);
    CHECK(info.ref_span() == 4);
    CHECK(info.query_start == 1);
    CHECK(info.query_end == 5);

    info = hamming_align("NAAAAAAAAAAAAAA", "TAAAATTTAAAAAAT", params);
    CHECK(info.cigar.to_string() == "8S6=1S");
    CHECK(info.edit_distance == 0);
    CHECK(info.sw_score == 6 * 3);
    CHECK(info.ref_start == 8);
    CHECK(info.ref_span() == 6);
    CHECK(info.query_start == 8);
    CHECK(info.query_end == 14);

    info = hamming_align("AAAAAAAAAAAAAAA", "TAAAAAATTTAAAAT", params);
    CHECK(info.cigar.to_string() == "1S6=8S");
    CHECK(info.edit_distance == 0);
    CHECK(info.sw_score == 6 * 3);
    CHECK(info.ref_start == 1);
    CHECK(info.ref_span() == 6);
    CHECK(info.query_start == 1);
    CHECK(info.query_end == 7);
}

TEST_CASE("highest_scoring_segment") {
    AlignmentParameters params{5, 7, 0, 0, 0, false, 0, 0};
    auto x = highest_scoring_segment("", "", params);
    CHECK(std::get<0>(x) == 0);
    CHECK(std::get<1>(x) == 0);
    x = highest_scoring_segment("AAAAAAAAAA", "AAAAAATTTT", params);
    CHECK(std::get<0>(x) == 0);
    CHECK(std::get<1>(x) == 6);
    x = highest_scoring_segment("AAAAAAAAAA", "TTTTAAAAAA", params);
    CHECK(std::get<0>(x) == 4);
    CHECK(std::get<1>(x) == 10);
    CHECK(highest_scoring_segment("AAAAAAAAAA", "AAAAAATTTT", params) == std::make_tuple(0ul, 6ul, 30));
    CHECK(highest_scoring_segment("AAAAAAAAAA", "TTTTAAAAAA", params) == std::make_tuple(4ul, 10ul, 30));
    CHECK(highest_scoring_segment("AAAAAAAAAA", "TTAAAAAATT", params) == std::make_tuple(2ul, 8ul, 30));
    CHECK(highest_scoring_segment("AAAAAAAAAAAAAAA", "TAAAAAATTTAAAAT", params) == std::make_tuple(1ul, 7ul, 30));
}

TEST_CASE("highest_scoring_segment with soft clipping") {
    AlignmentParameters params{2, 4, 0, 0, 5, false, 0, 0};
    auto x = highest_scoring_segment("", "", params);
    CHECK(std::get<0>(x) == 0);
    CHECK(std::get<1>(x) == 0);
    CHECK(std::get<2>(x) == 10);

    x = highest_scoring_segment("TAAT", "TAAA", params);
    CHECK(std::get<0>(x) == 0);
    CHECK(std::get<1>(x) == 4);
    CHECK(std::get<2>(x) == 3 * 2 - 4 + 10);

    x = highest_scoring_segment("AAA", "AAA", params);
    CHECK(std::get<0>(x) == 0);
    CHECK(std::get<1>(x) == 3);
    CHECK(std::get<2>(x) == 3 * 2 + 10);

    x = highest_scoring_segment("TAAT", "AAAA", params);
    CHECK(std::get<0>(x) == 0);
    CHECK(std::get<1>(x) == 4);
    CHECK(std::get<2>(x) == 2 * 2 - 2 * 4 + 10);

    x = highest_scoring_segment("ATAATA", "AAAAAA", params);
    CHECK(std::get<0>(x) == 0);
    CHECK(std::get<1>(x) == 6);
    CHECK(std::get<2>(x) == 4 * 2 - 2 * 4 + 10);

    x = highest_scoring_segment("TTAATA", "AAAAAA", params);
    CHECK(std::get<0>(x) == 2);
    CHECK(std::get<1>(x) == 6);
    CHECK(std::get<2>(x) == 3 * 2 - 1 * 4 + 5);

}

TEST_CASE("ssw align no result") {
    AlignmentParameters parameters{2, 8, 12, 1, 10, false, 0, 0};
    Aligner aligner{parameters};
    std::string query = "TCTCTCCCTCTCTCTCTCTCCCTCCCTCTCTCTCCCTCTCTCTCTCTCTCTCCCTCCCTT";
    std::string ref = "GAGGGAGAGAGAGAGAGGGAGAGAGAGAGAGAG";
    auto info = aligner.align(query, ref);
    CHECK(!info.has_value());
}
