
#ifndef VC_APX_H
#define VC_APX_H

#include "sr_apx/setmap/setmap.hpp"
#include "sr_apx/graph/graph.hpp"

namespace sr_apx::vc::apx {
	
Set dfs_apx(const Graph&);
Set std_apx(const Graph&);
Set heuristic_apx(const Graph&);

}

#endif
