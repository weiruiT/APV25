#ifndef PADRECONSTRUCTION_H
#define PADRECONSTRUCTION_H 1

#include "TString.h"

// Read selected BM channels from Tracker-step2-ana.root, merge adjacent pads
// separately for chip12 (BM1) and chip13 (BM2), and write BM1Events/BM2Events.
// Each tree contains only events with hits in that BM. A disabled or empty
// BM has a zero-entry tree, with no placeholder event=0 rows.
// Additionally reads the step4 Track tree and writes track-matched event
// arrival times (BM pads in the 3x3 window around the projected track hit
// and tracker strips in clusters near the predicted hit strips) into the
// TrackArrivalEvents tree of the same output file.
// trXY{A,E,S}_{mean,min,n}[layer] pool the matched X/Y strips of ONE
// detector (physical Tracker2 = index 1). ev* retain their legacy meaning
// of ALL three layers pooled. Stored times are in 25 ns sample units.
// BM event trees additionally retain cluster membership and per-cluster t20.
// TrackArrivalEvents.bmClusterS_* use all valid pads of a unique spatially
// matched cluster; bm_matched_cluster_count counts all spatial candidates.
// The existing bmS_* 3x3 pad timing is unchanged.
bool ReconstructPadEvents(TString step2Name, TString step4Name,
                          TString outputName, bool bm1Enabled = true,
                          bool bm2Enabled = true);
#endif
