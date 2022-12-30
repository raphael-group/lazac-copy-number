#ifndef _COPY_NUMBER_H
#define _COPY_NUMBER_H

#include <vector>
#include <string>

struct genomic_bin {
    std::string chromosome;
    std::string allele;
    int start;
    int end;

    genomic_bin() {};
    genomic_bin(std::string chromosome, std::string allele, int start, int end) :
        chromosome(chromosome), allele(allele), start(start), end(end) {};

    bool operator==(const genomic_bin& other) const {
        return (chromosome == other.chromosome &&
                allele == other.allele &&
                start == other.start &&
                end == other.end);
    }

    /* Lexicographic ordering */
    bool operator<(const genomic_bin& other) const {
        if (chromosome < other.chromosome) return true;
        if (chromosome > other.chromosome) return false;
        if (allele < other.allele) return true;
        if (allele > other.allele) return false;
        if (start < other.start) return true;
        if (start > other.start) return false;
        if (end < other.end) return true;
        if (end > other.end) return false;
        return false;
    }

    friend std::ostream& operator<<(std::ostream& os, const genomic_bin& bin);
};

struct copynumber_profile {
    std::vector<int> profile;
    std::vector<genomic_bin> bins;

    copynumber_profile() {};
    copynumber_profile(std::vector<int> profile, std::vector<genomic_bin> bins) :
        profile(profile), bins(bins) {};

    copynumber_profile operator=(const copynumber_profile &o) const {
        return copynumber_profile(o.profile, o.bins);
    }
};

struct breakpoint_profile {
    std::vector<int> profile;
    std::vector<genomic_bin> bins;

    breakpoint_profile() {};
    breakpoint_profile(std::vector<int> profile, std::vector<genomic_bin> bins) :
        profile(profile), bins(bins) {};

    breakpoint_profile operator=(const breakpoint_profile &o) const {
        return breakpoint_profile(o.profile, o.bins);
    }
};

/*
  Creates a breakpoint profile from a copy number profile
  where each chromosome and allele pair is considered seperately
  and genomic bins are sorted by their starting position.
 */
breakpoint_profile convert_to_breakpoint_profile(const copynumber_profile &p, int diploid_cn);

#endif
