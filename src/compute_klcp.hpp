#ifndef COMPUTE_KLCP_H
#define COMPUTE_KLCP_H

#include "defs.hpp"
#include "AppConfig.hpp"
#include "ReadsDB.hpp"

void print_lcpk(const unsigned& i, const unsigned& j, const ReadsDB& rdb,
                const ivec_t lcpKXY[2][2], const unsigned& k,
                std::ostream& lfs, const char* key = "lcpk");
void compute_klcp(ReadsDB& rdb, AppConfig& cfg);

#endif /* COMPUTE_KLCP_H */
