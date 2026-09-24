#ifndef _RLBWT_INTERVAL_ENCODING_HPP
#define _RLBWT_INTERVAL_ENCODING_HPP

#include "orbit/common.hpp"
#include "orbit/internal/move/interval_encoding_impl.hpp"
#include "orbit/internal/ds/alphabet.hpp"
#include "orbit/internal/rlbwt/rlbwt_helpers.hpp"

namespace orbit::rlbwt {

template<typename int_vector_t = int_vector_aligned, typename alphabet_t = nucleotide>
class rlbwt_interval_encoding_impl : public interval_encoding_impl<int_vector_t> {
    using base = interval_encoding_impl<int_vector_t>;
public:
    using alphabet_tag = alphabet_t;

    rlbwt_interval_encoding_impl() = default;

    static rlbwt_interval_encoding_impl lf_interval_encoding(const std::vector<uchar>& rlbwt_heads, const std::vector<ulint>& rlbwt_run_lengths, const split_params& sp = split_params()) {
        assert(rlbwt_heads.size() == rlbwt_run_lengths.size());

        rlbwt_interval_encoding_impl enc;

        auto [head_counts, n, max_length] = get_LF_head_counts(rlbwt_heads, rlbwt_run_lengths);
        enc.set_initial_values(n, rlbwt_heads.size(), max_length, sp);
        std::vector<ulint> img_rank_inv = get_LF_img_rank_inv(rlbwt_heads, head_counts);

        enc.init_img_rank_inv(rlbwt_run_lengths, img_rank_inv);

        enc.alphabet_ = alphabet_t(head_counts);
        enc.init_heads(rlbwt_heads, rlbwt_run_lengths);
        return enc;
    }

    // F runs are the BWT runs stably sorted by character, so BWT run i with
    // character c is F run C_head[c] + (number of earlier runs of c), where
    // C_head[c] counts the runs of smaller characters.  One pass fills packed
    // arrays of the F run lengths and of each BWT run's F run.
    static rlbwt_interval_encoding_impl fl_interval_encoding(const std::vector<uchar>& rlbwt_heads, const std::vector<ulint>& rlbwt_run_lengths, const split_params& sp = split_params()) {
        assert(rlbwt_heads.size() == rlbwt_run_lengths.size());

        rlbwt_interval_encoding_impl enc;

        auto [head_counts, n, max_length] = get_LF_head_counts(rlbwt_heads, rlbwt_run_lengths);
        const size_t r = rlbwt_heads.size();
        enc.set_initial_values(n, r, max_length, sp);

        std::vector<size_t> next_f(head_counts.size(), 0);
        std::vector<uchar> F_heads(r);
        for (size_t c = 0, seen = 0; c < head_counts.size(); ++c) {
            next_f[c] = seen;
            std::fill(F_heads.begin() + seen, F_heads.begin() + seen + head_counts[c], static_cast<uchar>(c));
            seen += head_counts[c];
        }
        {
            int_vector_t F_lens(r, bit_width(max_length));
            int_vector_t F_img_rank_inv(r, bit_width(r - 1));
            for (size_t i = 0; i < r; ++i) {
                const size_t f = next_f[rlbwt_heads[i]]++;
                F_lens.set(f, rlbwt_run_lengths[i]);
                F_img_rank_inv.set(i, f);
            }
            enc.init_img_rank_inv(F_lens, F_img_rank_inv);
            enc.alphabet_ = alphabet_t(head_counts);
            enc.init_heads(F_heads, F_lens);
        }
        return enc;
    }

    const int_vector_t& get_heads() const {
        return heads_;
    }

    const alphabet_t& get_alphabet() const {
        return alphabet_;
    }

    const ulint sigma() const {
        return static_cast<ulint>(alphabet_.size());
    }

private:
    alphabet_t alphabet_;
    int_vector_t heads_;

    template<typename collection1_t, typename collection2_t>
    void init_heads(const collection1_t& heads, const collection2_t& original_lengths) {
        heads_ = int_vector_t(base::intervals(), bit_width(alphabet_.size() - 1));
        size_t new_interval_idx = 0;
        for (size_t i = 0; i < original_lengths.size(); ++i) {
            ulint remaining = original_lengths[i];
            while (remaining > 0) {
                ulint chunk = base::get_length(new_interval_idx);
                heads_.set(new_interval_idx, alphabet_.map_char(heads[i]));
                remaining -= chunk;
                ++new_interval_idx;
            }
        }
        assert(new_interval_idx == base::intervals());
    }
};

} // namespace orbit::rlbwt

#endif /* end of include guard: _RLBWT_INTERVAL_ENCODING_HPP */