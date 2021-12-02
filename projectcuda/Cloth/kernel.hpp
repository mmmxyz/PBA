#pragma once

void Init(const fvec2* const RestPositionList, uint32_t* const TriIndList, const uint32_t* const InnerEdgeIndList, const fvec4* const InnerEdgeCList, const fmat2* const AList, const float* const VList, const uint32_t N);

void FemElasticProjectGPU(fvec3* const tempp, float* const ELambdaList, const uint32_t N, const float lambda, const float mu, const float mass, const float dt);
