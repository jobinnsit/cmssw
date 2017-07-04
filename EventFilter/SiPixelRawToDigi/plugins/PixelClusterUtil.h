/*Copyright 2017 Sushil Dubey, Shashi Dugad, TIFR

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

    http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
*/

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