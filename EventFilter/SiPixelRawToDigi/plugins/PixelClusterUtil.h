#ifndef CLUSTER_INFO
#define CLUSTER_INFO


typedef unsigned int uint; 
typedef unsigned long long uint64;

const uint EVENT_bits  = 11; // max 2k events 
const uint MODULE_bits = 11;
const uint XCOR_bits   = 9;
const uint YCOR_bits   = 9;

constexpr uint XCOR_shift   = 0;
constexpr uint YCOR_shift   = XCOR_shift + XCOR_bits; 
constexpr uint MODULE_shift = YCOR_shift + YCOR_bits; 
constexpr uint EVENT_shift  = MODULE_shift + MODULE_bits; 

constexpr uint64 XCOR_mask   = ~(~uint64(0) << XCOR_bits);
constexpr uint64 YCOR_mask   = ~(~uint64(0) << YCOR_bits);
constexpr uint64 MODULE_mask = ~(~uint64(0) << MODULE_bits);
constexpr uint64 EVENT_mask  = ~(~uint64(0) << EVENT_bits);

#endif 