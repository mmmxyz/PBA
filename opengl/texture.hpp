#pragma once

#include <cstdint>

class texture {
	uint32_t tex;
	uint32_t texw, texh;

    public:
	texture(const uint32_t& width, const uint32_t& height, const uint8_t* data);
	void activate(const uint32_t index) const;
	static void deactivate(const uint32_t index);
};
