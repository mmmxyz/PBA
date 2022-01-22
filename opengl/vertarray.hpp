#pragma once

#include <cstdint>

#include "utils/mathfunc/mathfunc.hpp"

struct vertex {
	float position[3];
	float color[4];
	float uv[2];
	float normal[3];
	uint32_t type;
};

//vaが同じ iff vaoが同じ iff vboが

class vertarray {
    protected:
	uint32_t size, isize;
	vertex* va;

	uint32_t vao = 0, vbo = 0, ibo = 0;

    public:
	vertarray()
	    : size(0)
	    , va(nullptr)
	    , isize(0)
	{
	}

	vertarray(uint32_t size, vertex* data, uint32_t isize = 0, uint32_t* ilist = nullptr);
	~vertarray();

	vertarray(const vertarray& v);
	vertarray(vertarray&& v)
	    : size(v.size)
	    , isize(v.isize)
	    , va(v.va)
	    , vao(v.vao)
	    , vbo(v.vbo)
	    , ibo(v.ibo)
	{
		v.va = nullptr;

		v.size	= 0;
		v.isize = 0;
		v.vao	= 0;
		v.vbo	= 0;
		v.ibo	= 0;
	}

	vertarray& operator=(const vertarray& v);
	vertarray& operator=(vertarray&& v)
	{
		if (vao == v.vao)
			return *this;

		size	= v.size;
		isize	= v.isize;
		v.size	= 0;
		v.isize = 0;

		va   = v.va;
		v.va = nullptr;

		vao = v.vao;
		vbo = v.vbo;
		ibo = v.ibo;

		v.vao = 0;
		v.vbo = 0;
		v.ibo = 0;

		return *this;
	}

	vertex& operator[](const uint32_t& index)
	{
		return va[index];
	}

	const vertex& operator[](const uint32_t& index) const
	{
		return va[index];
	}

	void resetvertarray(uint32_t size, vertex* data, uint32_t isize = 0, uint32_t* ilist = nullptr);
	void vboupdate();
	void setdata(vertex* data);
	void setilist(uint32_t* ilist);

	void setposition(uint32_t index, float x, float y, float z);
	void setposition(uint32_t index, fvec3 V);
	void setposition(uint32_t index, fvec2 V);

	void setallposition(const fvec3* PositionList, const uint32_t start, const uint32_t end);
	void setallposition(const fvec3* PositionList);
	void setallposition(const fvec3 Position, const uint32_t start, const uint32_t end);
	void setallposition(const fvec3 Position);

	void setallposition(const fvec2* PositionList, const uint32_t start, const uint32_t end);
	void setallposition(const fvec2* PositionList);
	void setallposition(const fvec2 Position, const uint32_t start, const uint32_t end);
	void setallposition(const fvec2 Position);

	void setuv(uint32_t index, float u, float v);
	void setnormal(uint32_t index, float nx, float ny, float nz);
	void setcolor(uint32_t index, float r, float g, float b, float alpha);
	void setcolor(float r, float g, float b, float alpha);
	void settype(uint32_t type);
	void bind() const;
	virtual void draw() const = 0;
	void unbind() const;
};
