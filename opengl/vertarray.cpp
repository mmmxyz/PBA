#include <iostream>

#include <glad/glad.h>
#include <KHR/khrplatform.h>

#include "opengl/vertarray.hpp"

#include "utils/mathfunc/mathfunc.hpp"

constexpr uint32_t REFCOUNTERSIZE	  = 256;
static uint8_t RefCounter[REFCOUNTERSIZE] = {};

vertarray::vertarray(uint32_t size, uint32_t isize, uint32_t* ilist)
    : size(size)
    , va(new vertex[size])
    , isize(isize)
{
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

	glGenVertexArrays(1, &vao);
	glBindVertexArray(vao);
	RefCounter[vao] = 1;

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
	}

	glBindVertexArray(0);
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);
}

vertarray::vertarray(const vertarray& v)
    : size(v.size)
    , isize(v.isize)
    , va(v.va)
    , vao(v.vao)
    , vbo(v.vbo)
    , ibo(v.ibo)
{
	RefCounter[vao]++;
}

vertarray& vertarray::operator=(const vertarray& v)
{

	if (vao == v.vao)
		return *this;

	size  = v.size;
	isize = v.isize;

	glBindVertexArray(vao);
	glBindBuffer(GL_ARRAY_BUFFER, 0);
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);
	glBindVertexArray(0);

	if (vao != 0)
		RefCounter[vao]--;
	if (vao != 0 && RefCounter[vao] == 0) {
		if (va != nullptr)
			delete[] va;
		glDeleteVertexArrays(1, &vao);
		glDeleteBuffers(1, &vbo);
		glDeleteBuffers(1, &ibo);
	}

	va = v.va;

	vao = v.vao;
	vbo = v.vbo;
	ibo = v.ibo;

	RefCounter[vao]++;

	return *this;
}

void vertarray::resetvertarray(uint32_t size, uint32_t isize, uint32_t* ilist)
{
	//std::cout << "start" << std::endl;
	//格納しているvao,vboの解放、vaの解放
	//glDeleteBuffersは0を入力に取ると何もしない。

	glBindVertexArray(vao);
	glBindBuffer(GL_ARRAY_BUFFER, 0);
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);
	glBindVertexArray(0);

	if (vao != 0)
		RefCounter[vao]--;

	if (vao != 0 && RefCounter[vao] == 0) {
		if (va != nullptr)
			delete[] va;
		glDeleteVertexArrays(1, &vao);
		glDeleteBuffers(1, &vbo);
		glDeleteBuffers(1, &ibo);
	}

	this->size  = size;
	this->isize = isize;

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

	glGenVertexArrays(1, &vao);
	glBindVertexArray(vao);
	RefCounter[vao]++;

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
	}

	glBindVertexArray(0);
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);

	//std::cout << "end" << std::endl;
}

void vertarray::vboupdate()
{
	//glBindVertexArray(vao);

	glBindBuffer(GL_ARRAY_BUFFER, vbo);
	glBufferData(GL_ARRAY_BUFFER, size * sizeof(vertex), va, GL_DYNAMIC_DRAW);
	glBindBuffer(GL_ARRAY_BUFFER, 0);

	//glBindVertexArray(0);
}

void vertarray::setilist(uint32_t* ilist)
{

	glBindVertexArray(vao);
	//glBindBuffer(GL_ARRAY_BUFFER, vbo);
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, ibo);
	glBufferData(GL_ELEMENT_ARRAY_BUFFER, isize * sizeof(GLuint), ilist, GL_DYNAMIC_DRAW);
	glBindVertexArray(0);
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);
}

vertarray::~vertarray()
{
	if (vao != 0)
		RefCounter[vao]--;

	glBindVertexArray(vao);
	glBindBuffer(GL_ARRAY_BUFFER, 0);
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);
	glBindVertexArray(0);

	if (vao != 0 && RefCounter[vao] == 0) {
		if (va != nullptr)
			delete[] va;
		glDeleteVertexArrays(1, &vao);
		glDeleteBuffers(1, &vbo);
		glDeleteBuffers(1, &ibo);
	}
}

void vertarray::bind() const
{
	glBindVertexArray(vao);
}

void vertarray::unbind() const
{
	glBindVertexArray(0);
}

void vertarray::checkRefCounter()
{
	for (uint32_t i = 0; i < REFCOUNTERSIZE; i++)
		if (RefCounter[i] != 0)
			std::cout << i << " " << uint32_t(RefCounter[i]) << std::endl;
}
