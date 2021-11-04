#include <glad/glad.h>
#include <KHR/khrplatform.h>

#include "opengl/vertarray.hpp"

#include "utils/mathfunc/mathfunc.hpp"

vertarray::vertarray()
    : size(0)
    , va(nullptr)
    , isize(0)
{
}

vertarray::vertarray(uint32_t size, vertex* data, uint32_t isize, uint32_t* ilist)
    : size(size)
    , va(nullptr)
    , isize(isize)
{

	if (data != nullptr) {
		va = data;
	} else {
		va = new vertex[size];
		for (uint32_t i = 0; i < size; i++) {
			va[i].position[0] = 0.0f;
			va[i].position[1] = 0.0f;
			va[i].position[2] = 0.0f;
			va[i].color[0]	  = 1.0f;
			va[i].color[1]	  = 1.0f;
			va[i].color[2]	  = 1.0f;
			va[i].color[3]	  = 1.0f;
			va[i].normal[0]	  = 0.0f;
			va[i].normal[1]	  = 0.0f;
			va[i].normal[2]	  = 0.0f;
			va[i].type	  = 0;
		}
	}

	glGenVertexArrays(1, &vao);
	glBindVertexArray(vao);

	glGenBuffers(1, &vbo);
	glBindBuffer(GL_ARRAY_BUFFER, vbo);
	glBufferData(GL_ARRAY_BUFFER, size * sizeof(vertex), va, GL_DYNAMIC_DRAW);

	glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, sizeof(vertex), reinterpret_cast<vertex*>(0)->position);
	glEnableVertexAttribArray(0);

	glVertexAttribPointer(1, 4, GL_FLOAT, GL_FALSE, sizeof(vertex), reinterpret_cast<vertex*>(0)->color);
	glEnableVertexAttribArray(1);

	glVertexAttribPointer(2, 2, GL_FLOAT, GL_FALSE, sizeof(vertex), reinterpret_cast<vertex*>(0)->uv);
	glEnableVertexAttribArray(2);

	glVertexAttribPointer(3, 3, GL_FLOAT, GL_FALSE, sizeof(vertex), reinterpret_cast<vertex*>(0)->normal);
	glEnableVertexAttribArray(3);

	glVertexAttribIPointer(4, 1, GL_UNSIGNED_INT, sizeof(vertex), (&reinterpret_cast<vertex*>(0)->type));
	glEnableVertexAttribArray(4);

	glBindBuffer(GL_ARRAY_BUFFER, 0);

	if (isize != 0) {
		glGenBuffers(1, &ibo);
		glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, ibo);
		glBufferData(GL_ELEMENT_ARRAY_BUFFER, isize * sizeof(GLuint), ilist, GL_DYNAMIC_DRAW);
		//glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);
	}

	glBindVertexArray(0);
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);
}

void vertarray::resetvertarray(uint32_t size, vertex* data, uint32_t isize, uint32_t* ilist)
{
	//格納しているvao,vboの解放、vaの解放
	//glDeleteBuffersは0を入力に取ると何もしない。
	if (va != nullptr)
		delete[] va;
	glDeleteBuffers(1, &vao);
	glDeleteBuffers(1, &vbo);
	glDeleteBuffers(1, &ibo);

	this->size  = size;
	this->isize = isize;

	if (data != nullptr) {
		va = data;
	} else {
		va = new vertex[size];
		for (uint32_t i = 0; i < size; i++) {
			va[i].position[0] = 0.0f;
			va[i].position[1] = 0.0f;
			va[i].position[2] = 0.0f;
			va[i].color[0]	  = 1.0f;
			va[i].color[1]	  = 1.0f;
			va[i].color[2]	  = 1.0f;
			va[i].color[3]	  = 1.0f;
			va[i].normal[0]	  = 0.0f;
			va[i].normal[1]	  = 0.0f;
			va[i].normal[2]	  = 0.0f;
			va[i].type	  = 0;
		}
	}

	glGenVertexArrays(1, &vao);
	glBindVertexArray(vao);

	glGenBuffers(1, &vbo);
	glBindBuffer(GL_ARRAY_BUFFER, vbo);
	glBufferData(GL_ARRAY_BUFFER, size * sizeof(vertex), va, GL_DYNAMIC_DRAW);

	glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, sizeof(vertex), reinterpret_cast<vertex*>(0)->position);
	glEnableVertexAttribArray(0);

	glVertexAttribPointer(1, 4, GL_FLOAT, GL_FALSE, sizeof(vertex), reinterpret_cast<vertex*>(0)->color);
	glEnableVertexAttribArray(1);

	glVertexAttribPointer(2, 2, GL_FLOAT, GL_FALSE, sizeof(vertex), reinterpret_cast<vertex*>(0)->uv);
	glEnableVertexAttribArray(2);

	glVertexAttribPointer(3, 3, GL_FLOAT, GL_FALSE, sizeof(vertex), reinterpret_cast<vertex*>(0)->normal);
	glEnableVertexAttribArray(3);

	glVertexAttribIPointer(4, 1, GL_UNSIGNED_INT, sizeof(vertex), (&reinterpret_cast<vertex*>(0)->type));
	glEnableVertexAttribArray(4);

	glBindBuffer(GL_ARRAY_BUFFER, 0);

	if (isize != 0) {
		glGenBuffers(1, &ibo);
		glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, ibo);
		glBufferData(GL_ELEMENT_ARRAY_BUFFER, isize * sizeof(GLuint), ilist, GL_DYNAMIC_DRAW);
		//glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);
	}

	glBindVertexArray(0);
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);
}

void vertarray::vboupdate()
{
	glBindVertexArray(vao);

	glBindBuffer(GL_ARRAY_BUFFER, vbo);
	glBufferData(GL_ARRAY_BUFFER, size * sizeof(vertex), va, GL_DYNAMIC_DRAW);
	glBindBuffer(GL_ARRAY_BUFFER, 0);

	glBindVertexArray(0);
}

void vertarray::setdata(vertex* data)
{
	if (va != nullptr)
		delete[] va;

	if (data != nullptr) {
		va = data;
	} else {
		va = new vertex[size];
		for (uint32_t i = 0; i < size; i++) {
			va[i].position[0] = 0.0f;
			va[i].position[1] = 0.0f;
			va[i].position[2] = 0.0f;
			va[i].color[0]	  = 1.0f;
			va[i].color[1]	  = 1.0f;
			va[i].color[2]	  = 1.0f;
			va[i].color[3]	  = 1.0f;
			va[i].normal[0]	  = 0.0f;
			va[i].normal[1]	  = 0.0f;
			va[i].normal[2]	  = 0.0f;
			va[i].type	  = 0;
		}
	}

	glBindVertexArray(vao);

	glBindBuffer(GL_ARRAY_BUFFER, vbo);
	glBufferData(GL_ARRAY_BUFFER, size * sizeof(vertex), va, GL_DYNAMIC_DRAW);
	glBindBuffer(GL_ARRAY_BUFFER, 0);

	glBindVertexArray(0);
}

void vertarray::setilist(uint32_t* ilist)
{

	glBindBuffer(GL_ARRAY_BUFFER, vbo);
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, ibo);
	glBufferData(GL_ELEMENT_ARRAY_BUFFER, isize * sizeof(GLuint), ilist, GL_DYNAMIC_DRAW);
	glBindVertexArray(0);
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);
}

void vertarray::setposition(uint32_t index, float x, float y, float z)
{

	va[index].position[0] = x;
	va[index].position[1] = y;
	va[index].position[2] = z;

	glBindVertexArray(vao);

	glBindBuffer(GL_ARRAY_BUFFER, vbo);
	glBufferSubData(GL_ARRAY_BUFFER, index * sizeof(vertex), sizeof(vertex), va + index);
	glBindBuffer(GL_ARRAY_BUFFER, 0);

	glBindVertexArray(0);
}

void vertarray::setposition(uint32_t index, fvec3 V)
{
	this->setposition(index, V.x, V.y, V.z);
}

void vertarray::setuv(uint32_t index, float u, float v)
{
	va[index].uv[0] = u;
	va[index].uv[1] = v;

	glBindVertexArray(vao);

	glBindBuffer(GL_ARRAY_BUFFER, vbo);
	glBufferSubData(GL_ARRAY_BUFFER, index * sizeof(vertex), sizeof(vertex), va + index);
	glBindBuffer(GL_ARRAY_BUFFER, 0);

	glBindVertexArray(0);
}

void vertarray::setnormal(uint32_t index, float nx, float ny, float nz)
{

	va[index].normal[0] = nx;
	va[index].normal[1] = ny;
	va[index].normal[2] = nz;

	glBindVertexArray(vao);

	glBindBuffer(GL_ARRAY_BUFFER, vbo);
	glBufferSubData(GL_ARRAY_BUFFER, index * sizeof(vertex), sizeof(vertex), va + index);
	glBindBuffer(GL_ARRAY_BUFFER, 0);

	glBindVertexArray(0);
}

void vertarray::setcolor(uint32_t index, float r, float g, float b, float alpha)
{

	va[index].color[0] = r;
	va[index].color[1] = g;
	va[index].color[2] = b;
	va[index].color[3] = alpha;

	glBindVertexArray(vao);

	glBindBuffer(GL_ARRAY_BUFFER, vbo);
	glBufferSubData(GL_ARRAY_BUFFER, index * sizeof(vertex), sizeof(vertex), va + index);
	glBindBuffer(GL_ARRAY_BUFFER, 0);

	glBindVertexArray(0);
}

void vertarray::setcolor(float r, float g, float b, float alpha)
{

	for (uint32_t i = 0; i < size; i++) {
		va[i].color[0] = r;
		va[i].color[1] = g;
		va[i].color[2] = b;
		va[i].color[3] = alpha;
	}

	glBindVertexArray(vao);

	glBindBuffer(GL_ARRAY_BUFFER, vbo);
	glBufferData(GL_ARRAY_BUFFER, size * sizeof(vertex), va, GL_DYNAMIC_DRAW);
	glBindBuffer(GL_ARRAY_BUFFER, 0);

	glBindVertexArray(0);
}

void vertarray::settype(uint32_t type)
{

	for (uint32_t i = 0; i < size; i++) {
		va[i].type = type;
	}

	glBindVertexArray(vao);

	glBindBuffer(GL_ARRAY_BUFFER, vbo);
	glBufferData(GL_ARRAY_BUFFER, size * sizeof(vertex), va, GL_DYNAMIC_DRAW);
	glBindBuffer(GL_ARRAY_BUFFER, 0);

	glBindVertexArray(0);
}

vertarray::~vertarray()
{
	if (va != nullptr)
		delete[] va;
	glDeleteBuffers(1, &vao);
	glDeleteBuffers(1, &vbo);
	glDeleteBuffers(1, &ibo);
}

void vertarray::deletedata()
{
	if (va != nullptr)
		delete[] va;
	glDeleteBuffers(1, &vao);
	glDeleteBuffers(1, &vbo);
	glDeleteBuffers(1, &ibo);
}

void vertarray::bind() const
{
	glBindVertexArray(vao);
}

void vertarray::unbind() const
{
	glBindVertexArray(0);
}
