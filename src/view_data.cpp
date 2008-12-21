
const int num_pal_entries = 256;

const float palette_data[num_pal_entries+1][3] = 
{
  { 0.0000, 0.0000, 0.4650 }, { 0.0000, 0.0012, 0.4760 }, { 0.0000, 0.0024, 0.4870 }, { 0.0000, 0.0036, 0.4980 }, 
  { 0.0000, 0.0048, 0.5089 }, { 0.0000, 0.0060, 0.5199 }, { 0.0000, 0.0072, 0.5309 }, { 0.0000, 0.0084, 0.5419 }, 
  { 0.0000, 0.0096, 0.5529 }, { 0.0000, 0.0108, 0.5639 }, { 0.0000, 0.0120, 0.5749 }, { 0.0000, 0.0132, 0.5858 }, 
  { 0.0000, 0.0143, 0.5968 }, { 0.0000, 0.0170, 0.6078 }, { 0.0000, 0.0215, 0.6188 }, { 0.0000, 0.0259, 0.6298 }, 
  { 0.0000, 0.0303, 0.6408 }, { 0.0000, 0.0348, 0.6511 }, { 0.0000, 0.0392, 0.6611 }, { 0.0000, 0.0471, 0.6710 }, 
  { 0.0000, 0.0557, 0.6810 }, { 0.0000, 0.0643, 0.6910 }, { 0.0000, 0.0730, 0.7009 }, { 0.0000, 0.0816, 0.7109 }, 
  { 0.0000, 0.0902, 0.7208 }, { 0.0000, 0.0988, 0.7308 }, { 0.0000, 0.1074, 0.7407 }, { 0.0000, 0.1163, 0.7507 }, 
  { 0.0000, 0.1272, 0.7607 }, { 0.0000, 0.1380, 0.7706 }, { 0.0000, 0.1489, 0.7814 }, { 0.0000, 0.1597, 0.7927 }, 
  { 0.0000, 0.1706, 0.8041 }, { 0.0000, 0.1814, 0.8154 }, { 0.0000, 0.1923, 0.8268 }, { 0.0000, 0.2031, 0.8381 }, 
  { 0.0000, 0.2140, 0.8495 }, { 0.0000, 0.2248, 0.8608 }, { 0.0000, 0.2357, 0.8722 }, { 0.0000, 0.2471, 0.8836 }, 
  { 0.0000, 0.2590, 0.8949 }, { 0.0000, 0.2708, 0.9046 }, { 0.0000, 0.2827, 0.9129 }, { 0.0000, 0.2945, 0.9211 }, 
  { 0.0000, 0.3064, 0.9294 }, { 0.0000, 0.3183, 0.9377 }, { 0.0000, 0.3308, 0.9460 }, { 0.0000, 0.3444, 0.9543 }, 
  { 0.0000, 0.3579, 0.9626 }, { 0.0000, 0.3714, 0.9705 }, { 0.0000, 0.3849, 0.9756 }, { 0.0000, 0.3984, 0.9807 }, 
  { 0.0000, 0.4120, 0.9858 }, { 0.0000, 0.4277, 0.9909 }, { 0.0000, 0.4440, 0.9960 }, { 0.0000, 0.4604, 1.0000 }, 
  { 0.0000, 0.4767, 1.0000 }, { 0.0000, 0.4931, 1.0000 }, { 0.0000, 0.5094, 1.0000 }, { 0.0000, 0.5257, 1.0000 }, 
  { 0.0000, 0.5421, 1.0000 }, { 0.0000, 0.5584, 1.0000 }, { 0.0000, 0.5748, 1.0000 }, { 0.0000, 0.5911, 1.0000 }, 
  { 0.0000, 0.6074, 1.0000 }, { 0.0000, 0.6239, 1.0000 }, { 0.0000, 0.6409, 1.0000 }, { 0.0000, 0.6579, 1.0000 }, 
  { 0.0000, 0.6748, 1.0000 }, { 0.0000, 0.6918, 1.0000 }, { 0.0000, 0.7087, 1.0000 }, { 0.0000, 0.7257, 1.0000 }, 
  { 0.0000, 0.7427, 1.0000 }, { 0.0000, 0.7596, 1.0000 }, { 0.0000, 0.7766, 1.0000 }, { 0.0000, 0.7919, 1.0000 }, 
  { 0.0000, 0.8057, 1.0000 }, { 0.0000, 0.8194, 1.0000 }, { 0.0000, 0.8332, 1.0000 }, { 0.0000, 0.8469, 1.0000 }, 
  { 0.0000, 0.8606, 1.0000 }, { 0.0000, 0.8744, 1.0000 }, { 0.0000, 0.8872, 1.0000 }, { 0.0000, 0.8992, 0.9990 }, 
  { 0.0000, 0.9113, 0.9960 }, { 0.0000, 0.9234, 0.9929 }, { 0.0000, 0.9355, 0.9898 }, { 0.0000, 0.9465, 0.9867 }, 
  { 0.0000, 0.9536, 0.9819 }, { 0.0000, 0.9607, 0.9749 }, { 0.0000, 0.9678, 0.9679 }, { 0.0000, 0.9749, 0.9609 }, 
  { 0.0000, 0.9820, 0.9511 }, { 0.0000, 0.9866, 0.9408 }, { 0.0000, 0.9894, 0.9306 }, { 0.0000, 0.9922, 0.9203 }, 
  { 0.0000, 0.9950, 0.9100 }, { 0.0000, 0.9978, 0.8950 }, { 0.0000, 1.0000, 0.8800 }, { 0.0000, 1.0000, 0.8649 }, 
  { 0.0000, 1.0000, 0.8449 }, { 0.0000, 1.0000, 0.8224 }, { 0.0000, 1.0000, 0.8000 }, { 0.0000, 1.0000, 0.7775 }, 
  { 0.0000, 1.0000, 0.7551 }, { 0.0001, 1.0000, 0.7323 }, { 0.0029, 1.0000, 0.7092 }, { 0.0057, 1.0000, 0.6861 }, 
  { 0.0085, 1.0000, 0.6630 }, { 0.0113, 1.0000, 0.6399 }, { 0.0141, 1.0000, 0.6168 }, { 0.0169, 1.0000, 0.5913 }, 
  { 0.0196, 1.0000, 0.5656 }, { 0.0266, 1.0000, 0.5398 }, { 0.0342, 1.0000, 0.5141 }, { 0.0418, 1.0000, 0.4883 }, 
  { 0.0494, 1.0000, 0.4626 }, { 0.0586, 1.0000, 0.4371 }, { 0.0721, 1.0000, 0.4116 }, { 0.0856, 1.0000, 0.3860 }, 
  { 0.0991, 1.0000, 0.3605 }, { 0.1127, 1.0000, 0.3349 }, { 0.1262, 1.0000, 0.3094 }, { 0.1397, 1.0000, 0.2855 }, 
  { 0.1581, 1.0000, 0.2625 }, { 0.1795, 1.0000, 0.2395 }, { 0.2010, 1.0000, 0.2165 }, { 0.2225, 1.0000, 0.1936 }, 
  { 0.2440, 1.0000, 0.1706 }, { 0.2665, 1.0000, 0.1476 }, { 0.2901, 1.0000, 0.1246 }, { 0.3137, 1.0000, 0.1017 }, 
  { 0.3373, 1.0000, 0.0912 }, { 0.3609, 1.0000, 0.0818 }, { 0.3845, 1.0000, 0.0723 }, { 0.4081, 1.0000, 0.0629 }, 
  { 0.4317, 1.0000, 0.0534 }, { 0.4553, 1.0000, 0.0440 }, { 0.4789, 1.0000, 0.0345 }, { 0.5025, 1.0000, 0.0251 }, 
  { 0.5261, 1.0000, 0.0220 }, { 0.5500, 1.0000, 0.0189 }, { 0.5750, 1.0000, 0.0159 }, { 0.6000, 1.0000, 0.0128 }, 
  { 0.6249, 1.0000, 0.0098 }, { 0.6499, 1.0000, 0.0067 }, { 0.6749, 1.0000, 0.0037 }, { 0.6999, 1.0000, 0.0006 }, 
  { 0.7249, 1.0000, 0.0000 }, { 0.7499, 1.0000, 0.0000 }, { 0.7748, 1.0000, 0.0000 }, { 0.7998, 1.0000, 0.0000 }, 
  { 0.8235, 1.0000, 0.0000 }, { 0.8417, 1.0000, 0.0000 }, { 0.8600, 1.0000, 0.0000 }, { 0.8782, 1.0000, 0.0000 }, 
  { 0.8946, 1.0000, 0.0000 }, { 0.9076, 1.0000, 0.0000 }, { 0.9206, 1.0000, 0.0000 }, { 0.9336, 1.0000, 0.0000 }, 
  { 0.9439, 1.0000, 0.0000 }, { 0.9516, 1.0000, 0.0000 }, { 0.9592, 1.0000, 0.0000 }, { 0.9668, 1.0000, 0.0000 }, 
  { 0.9745, 1.0000, 0.0000 }, { 0.9821, 1.0000, 0.0000 }, { 0.9863, 1.0000, 0.0000 }, { 0.9883, 0.9991, 0.0000 }, 
  { 0.9903, 0.9976, 0.0000 }, { 0.9923, 0.9961, 0.0000 }, { 0.9943, 0.9946, 0.0000 }, { 0.9964, 0.9931, 0.0000 }, 
  { 0.9984, 0.9916, 0.0000 }, { 1.0000, 0.9901, 0.0000 }, { 1.0000, 0.9839, 0.0000 }, { 1.0000, 0.9773, 0.0000 }, 
  { 1.0000, 0.9708, 0.0000 }, { 1.0000, 0.9643, 0.0000 }, { 1.0000, 0.9567, 0.0000 }, { 1.0000, 0.9470, 0.0000 }, 
  { 1.0000, 0.9372, 0.0000 }, { 1.0000, 0.9274, 0.0000 }, { 1.0000, 0.9177, 0.0000 }, { 1.0000, 0.9079, 0.0000 }, 
  { 1.0000, 0.8981, 0.0000 }, { 1.0000, 0.8855, 0.0000 }, { 1.0000, 0.8716, 0.0000 }, { 1.0000, 0.8576, 0.0000 }, 
  { 1.0000, 0.8437, 0.0000 }, { 1.0000, 0.8297, 0.0000 }, { 1.0000, 0.8155, 0.0000 }, { 1.0000, 0.8005, 0.0000 }, 
  { 1.0000, 0.7856, 0.0000 }, { 1.0000, 0.7707, 0.0000 }, { 1.0000, 0.7557, 0.0000 }, { 1.0000, 0.7408, 0.0000 }, 
  { 1.0000, 0.7258, 0.0000 }, { 1.0000, 0.7109, 0.0000 }, { 1.0000, 0.6960, 0.0000 }, { 1.0000, 0.6794, 0.0000 }, 
  { 1.0000, 0.6618, 0.0000 }, { 1.0000, 0.6442, 0.0000 }, { 1.0000, 0.6265, 0.0000 }, { 1.0000, 0.6089, 0.0000 }, 
  { 1.0000, 0.5913, 0.0000 }, { 1.0000, 0.5737, 0.0000 }, { 1.0000, 0.5560, 0.0000 }, { 1.0000, 0.5384, 0.0000 }, 
  { 1.0000, 0.5208, 0.0000 }, { 1.0000, 0.5028, 0.0000 }, { 1.0000, 0.4820, 0.0000 }, { 1.0000, 0.4611, 0.0000 }, 
  { 0.9988, 0.4402, 0.0000 }, { 0.9948, 0.4194, 0.0000 }, { 0.9907, 0.3985, 0.0000 }, { 0.9866, 0.3777, 0.0000 }, 
  { 0.9826, 0.3568, 0.0000 }, { 0.9785, 0.3359, 0.0000 }, { 0.9739, 0.3151, 0.0000 }, { 0.9661, 0.2942, 0.0000 }, 
  { 0.9583, 0.2733, 0.0000 }, { 0.9504, 0.2570, 0.0000 }, { 0.9426, 0.2415, 0.0000 }, { 0.9347, 0.2260, 0.0000 }, 
  { 0.9231, 0.2104, 0.0000 }, { 0.9115, 0.1949, 0.0000 }, { 0.8999, 0.1794, 0.0000 }, { 0.8883, 0.1639, 0.0000 }, 
  { 0.8767, 0.1484, 0.0000 }, { 0.8651, 0.1329, 0.0000 }, { 0.8535, 0.1174, 0.0000 }, { 0.8419, 0.1019, 0.0000 }, 
  { 0.8288, 0.0864, 0.0000 }, { 0.8154, 0.0709, 0.0000 }, { 0.8020, 0.0553, 0.0000 }, { 0.7885, 0.0398, 0.0000 }, 
  { 0.7751, 0.0243, 0.0000 }, { 0.7617, 0.0088, 0.0000 }, { 0.7483, 0.0000, 0.0000 }, { 0.7348, 0.0000, 0.0000 }, 
  { 0.7210, 0.0000, 0.0000 }, { 0.7068, 0.0000, 0.0000 }, { 0.6927, 0.0000, 0.0000 }, { 0.6786, 0.0000, 0.0000 }, 
  { 0.6645, 0.0000, 0.0000 }, { 0.6503, 0.0000, 0.0000 }, { 0.6362, 0.0000, 0.0000 }, { 0.6221, 0.0000, 0.0000 }, 
  { 0.6080, 0.0000, 0.0000 }, { 0.5939, 0.0000, 0.0000 }, { 0.5797, 0.0000, 0.0000 }, { 0.5656, 0.0000, 0.0000 }, 
  { 0.5515, 0.0000, 0.0000 }, { 0.5374, 0.0000, 0.0000 }, { 0.5232, 0.0000, 0.0000 }, { 0.5091, 0.0000, 0.0000 }, 
  { 0.4950, 0.0000, 0.0000 }
};

