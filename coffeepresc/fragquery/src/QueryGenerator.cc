#include "QueryGenerator.hpp"

#include <set>

namespace query {
  // index is 1-origin (only positive)
  int QueryGenerator::distToIndex(fltype dist) const {
    int q = dist/width;
    fltype r = dist - q*width;
    if (r >= width/2) q++;
		return (q <= 0) ? 1 : q;
	}

  pair_info QueryGenerator::getPairInfo(const std::vector<pair_info>& scores) const {
    std::vector<pair_info> sorted_scores = scores;
		sort(sorted_scores.begin(), sorted_scores.end());
		pair_info pinfo = pair_info();
    if (scores.size() > 0) pinfo = sorted_scores[0]; // best pair_score
		return pinfo;
	}

  void QueryGenerator::append(fltype dist, int rank1, int rank2, fltype score1, fltype score2, fltype pair_score) {
    fltype range_min = dist - radius*2;
    fltype range_max = dist + radius*2;

    if (range_min > distance_max || range_max < distance_min) return;
    range_min = std::max(range_min, distance_min);
    range_max = std::min(range_max, distance_max);

    pair_info pinfo(rank1, rank2, score1, score2, pair_score);
    eves.push_back(dist_event(range_min        , pinfo, EventType::START));
    eves.push_back(dist_event(range_max + width, pinfo, EventType::END));
  }

  std::vector<output_query> QueryGenerator::getPairInfoVec() {
    std::multiset<pair_info> scores;
    std::vector<output_query> ret;

    std::sort(eves.begin(), eves.end());

    bool event = false;
    int e_ind = 0, d_ind = 0;
    while (e_ind < eves.size()) {
      while (distToIndex(eves[e_ind].dist) <= d_ind && e_ind < eves.size()) {
        if (eves[e_ind].type == EventType::START) {
          scores.insert(eves[e_ind].pinfo);
        } else {
          scores.erase(scores.find(eves[e_ind].pinfo));
        }
        e_ind++;
        event = true;
      }
      
      if (!scores.empty()) {
        pair_info pinfo = getPairInfo(std::vector<pair_info>(scores.begin(), scores.end()));
        if (!ret.empty() && !event) {
          ret.back().dist_max = indexToDist(d_ind);
        } else {
          ret.push_back(output_query(indexToDist(d_ind), indexToDist(d_ind), pinfo));
        }
      }
      d_ind++;
      event = false;
    }
    
    for (int i = 1; i < ret.size(); i++) {
      if (ret[i-1].is_next(ret[i], width)) {
        ret[i-1].dist_max = ret[i].dist_max;
        ret.erase(ret.begin() + i);
        i--;
      }
    }

    return ret;
  }
}