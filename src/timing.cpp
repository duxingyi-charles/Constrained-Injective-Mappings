//
// Created by Charles Du on 8/15/21.
//

#include "timing.h"

Time_duration global_TLC_time;
Time_duration global_arc_seg_time;
Time_duration global_arc_occupancy_time;

void reset_timings() {
    global_TLC_time = Time_duration::zero();
    global_arc_seg_time = Time_duration::zero();
    global_arc_occupancy_time = Time_duration::zero();
}
