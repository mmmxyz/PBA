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

	vertarray(uint32_t size, uint32_t isize = 0, uint32_t* ilist = nullptr);
	vertarray(uint32_t size)
	    : vertarray(size, 0, nullptr)
	{
	}

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

	void resetvertarray(uint32_t size, uint32_t isize = 0, uint32_t* ilist = nullptr);
	void vboupdate();
	void setilist(uint32_t* ilist);

	void setposition(uint32_t index, float x, float y, float z)
	{
		va[index].position[0] = x;
		va[index].position[1] = y;
		va[index].position[2] = z;
	}
	void setposition(uint32_t index, fvec3 V)
	{
		this->setposition(index, V.x, V.y, V.z);
	}
	void setposition(uint32_t index, fvec2 V)
	{
		this->setposition(index, V.x, V.y, 0.0);
	}

	void setallposition(const fvec3* PositionList, const uint32_t start, const uint32_t end)
	{
		//end > sizeだとやばいことになる
		for (uint32_t i = start; i < end; i++) {
			va[i].position[0] = PositionList[i].x;
			va[i].position[1] = PositionList[i].y;
			va[i].position[2] = PositionList[i].z;
		}
	}
	void setallposition(const fvec3* PositionList)
	{
		for (uint32_t i = 0; i < size; i++) {
			va[i].position[0] = PositionList[i].x;
			va[i].position[1] = PositionList[i].y;
			va[i].position[2] = PositionList[i].z;
		}
	}
	void setallposition(const fvec3 Position, const uint32_t start, const uint32_t end)
	{
		//end > sizeだとやばいことになる
		for (uint32_t i = start; i < end; i++) {
			va[i].position[0] = Position.x;
			va[i].position[1] = Position.y;
			va[i].position[2] = Position.z;
		}
	}
	void setallposition(const fvec3 Position)
	{
		for (uint32_t i = 0; i < size; i++) {
			va[i].position[0] = Position.x;
			va[i].position[1] = Position.y;
			va[i].position[2] = Position.z;
		}
	}

	void setallposition(const fvec2* PositionList, const uint32_t start, const uint32_t end)
	{
		//end > sizeだとやばいことになる
		for (uint32_t i = start; i < end; i++) {
			va[i].position[0] = PositionList[i].x;
			va[i].position[1] = PositionList[i].y;
			va[i].position[2] = 0.0;
		}
	}
	void setallposition(const fvec2* PositionList)
	{
		for (uint32_t i = 0; i < size; i++) {
			va[i].position[0] = PositionList[i].x;
			va[i].position[1] = PositionList[i].y;
			va[i].position[2] = 0.0;
		}
	}
	void setallposition(const fvec2 Position, const uint32_t start, const uint32_t end)
	{
		//end > sizeだとやばいことになる
		for (uint32_t i = start; i < end; i++) {
			va[i].position[0] = Position.x;
			va[i].position[1] = Position.y;
			va[i].position[2] = 0.0;
		}
	}
	void setallposition(const fvec2 Position)
	{
		for (uint32_t i = 0; i < size; i++) {
			va[i].position[0] = Position.x;
			va[i].position[1] = Position.y;
			va[i].position[2] = 0.0;
		}
	}

	void setuv(uint32_t index, float u, float v)
	{
		va[index].uv[0] = u;
		va[index].uv[1] = v;
	}
	void setnormal(uint32_t index, float nx, float ny, float nz)
	{
		va[index].normal[0] = nx;
		va[index].normal[1] = ny;
		va[index].normal[2] = nz;
	}
	void setcolor(uint32_t index, float r, float g, float b, float alpha)
	{
		va[index].color[0] = r;
		va[index].color[1] = g;
		va[index].color[2] = b;
		va[index].color[3] = alpha;
	}
	void setcolor(float r, float g, float b, float alpha)
	{
		for (uint32_t i = 0; i < size; i++) {
			va[i].color[0] = r;
			va[i].color[1] = g;
			va[i].color[2] = b;
			va[i].color[3] = alpha;
		}
	}
	void settype(uint32_t type)
	{
		for (uint32_t i = 0; i < size; i++) {
			va[i].type = type;
		}
	}

	void bind() const;
	virtual void draw() const = 0;
	void unbind() const;

	uint32_t getsize() const
	{
		return size;
	}
	uint32_t getisize() const
	{
		return isize;
	}

	static void checkRefCounter();
};
