#ifndef QUERY_H_
#define QUERY_H_

#include "common.hpp"
#include "infile_reader.hpp"
#include "utils.hpp"

#include <vector>

namespace query {
  enum EventType { START, END };

  struct pair_info {
    int rank1, rank2;
    fltype score1, score2, pair_score;
    pair_info() : rank1(0), rank2(0), score1(0), score2(0), pair_score(0) {}
    pair_info(int rank1, int rank2, fltype score1, fltype score2, fltype pair_score)
      : rank1(rank1), rank2(rank2), score1(score1), score2(score2), pair_score(pair_score) {}
    bool operator<(const pair_info& o) const { return pair_score < o.pair_score; }
  };

  struct dist_event {
    fltype dist;
    pair_info pinfo;
    EventType type;
    dist_event(fltype dist, pair_info pinfo, EventType type) : dist(dist), pinfo(pinfo), type(type) {}
    bool operator<(const dist_event& o) const { return dist < o.dist; }
  };

  struct output_query {
    fltype dist_min, dist_max;
    int rank1, rank2;
    fltype score1, score2, pair_score;
    output_query(fltype dist_min, fltype dist_max, pair_info pinfo)
      : dist_min(dist_min), dist_max(dist_max), rank1(pinfo.rank1), rank2(pinfo.rank2), score1(pinfo.score1), score2(pinfo.score2), pair_score(pinfo.pair_score) {}
    bool is_next(const output_query& o, fltype width) const { return abs(dist_max - o.dist_min) < width + EPS && abs(pair_score - o.pair_score) < EPS; }
  };

	class QueryGenerator {
    const fltype width, radius, distance_min, distance_max;
		std::vector<dist_event> eves;
		int distToIndex(fltype dist) const; // index is 1-origin (only positive)
    fltype indexToDist(int index) const { return index*width; }
    pair_info getPairInfo(const std::vector<pair_info>& scores) const;
	public:
		QueryGenerator(const format::QueryParams& query_params)
      : width(query_params.distance_width), radius(query_params.cluster_size), distance_min(query_params.distance_min), distance_max(query_params.distance_max) {}
		void append(fltype dist, int rank1, int rank2, fltype score1, fltype score2, fltype pair_score);
		std::vector<output_query> getPairInfoVec();
	};
}
#endif