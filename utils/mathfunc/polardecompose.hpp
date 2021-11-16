#pragma once

#include "utils/mathfunc/mathfunc.hpp"

//fquaternion ExtractRotation(const fmat3& A, const uint32_t& numitr, const fquaternion initq = fquaternion(0.0, 0.0, 0.0, 1.0));
fquaternion ExtractRotation(const fmat3& A, const uint32_t& numitr, const fquaternion initq = fquaternion(0.0, 0.0, 0.0, 1.0));
float ExtractRotation(const fmat2& A, const uint32_t& numitr, const float& initomega = 0.0);
