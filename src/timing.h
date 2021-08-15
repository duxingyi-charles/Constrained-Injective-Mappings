//
// Created by Charles Du on 8/15/21.
//

#ifndef TLC_TIMING_H
#define TLC_TIMING_H

#include <chrono>
typedef std::chrono::duration<double> Time_duration;

// computing TLC func/grad
extern Time_duration global_TLC_time;
// computing arc segment area func/grad
extern Time_duration global_arc_seg_time;
// computing arc occupancy func/grad
extern Time_duration global_arc_occupancy_time;


void reset_timings();



#endif //TLC_TIMING_H
