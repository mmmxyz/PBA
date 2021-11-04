#pragma once

class shader {

    public:
	uint32_t program;
	shader(const char* vsname, const char* fsname);
	shader(void);
	void setprogram(const char* vsname, const char* fsname);

	void useprogram() const;
};
