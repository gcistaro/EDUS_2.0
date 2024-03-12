#include <iostream>
#include "core/print_timing.hpp"

void print_timing(const int& val)
{
    auto timing_result = global_rtgraph_timer.process();

    if (val == 1) {
        std::cout << timing_result.print({rt_graph::Stat::Count, rt_graph::Stat::Total, rt_graph::Stat::Percentage,
                                          rt_graph::Stat::SelfPercentage, rt_graph::Stat::Median,
                                          rt_graph::Stat::Min, rt_graph::Stat::Max});
    }
    if (val == 2) {
        timing_result = timing_result.flatten(1).sort_nodes();
        std::cout << timing_result.print({rt_graph::Stat::Count, rt_graph::Stat::Total, rt_graph::Stat::Percentage,
                                          rt_graph::Stat::SelfPercentage, rt_graph::Stat::Median,
                                          rt_graph::Stat::Min, rt_graph::Stat::Max});
    }   
}

