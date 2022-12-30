#include "copy_number.hpp"
#include "vec_utilities.hpp"

#include <vector>
#include <map>
#include <ostream>

std::ostream& operator<<(std::ostream& os, const genomic_bin& bin) {
    os << bin.chromosome << ":" << bin.allele << ":" << bin.start << "-" << bin.end;
    return os;
}

breakpoint_profile convert_to_breakpoint_profile(const copynumber_profile &p, int diploid_cn) {
    std::map<std::pair<std::string, std::string>, copynumber_profile> chrom_allele_profiles;
    for (size_t i = 0; i < p.bins.size(); i++) {
        auto chrom_allele = std::make_pair(p.bins[i].chromosome, p.bins[i].allele);
        if (!chrom_allele_profiles.count(chrom_allele)) {
            chrom_allele_profiles[chrom_allele] = copynumber_profile();
        }

        chrom_allele_profiles[chrom_allele].bins.push_back(p.bins[i]);
        chrom_allele_profiles[chrom_allele].profile.push_back(p.profile[i]);
    }

    breakpoint_profile bp;
    for (const auto &[chrom_allele, cn_profile] : chrom_allele_profiles) {
        std::vector<size_t> index_vector = argsort(cn_profile.bins);
        std::vector<genomic_bin> bins = select(cn_profile.bins, index_vector);
        std::vector<int> profile = select(cn_profile.profile, index_vector);

        std::vector<int> bp_profile(profile.size());
        for (size_t i = 0; i < profile.size(); i++) {
            if (i == 0) {
                bp_profile[i] = profile[i] - diploid_cn;
            } else {
                bp_profile[i] = profile[i] - profile[i - 1];
            }

            bp.profile.push_back(bp_profile[i]);
            bp.bins.push_back(bins[i]);
        }
    }

    return bp;
}
