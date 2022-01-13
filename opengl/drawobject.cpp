#include <glad/glad.h>
#include <KHR/khrplatform.h>
#define GLFW_INCLUDE_NONE
#include <GLFW/glfw3.h>

#include <cmath>

#include "opengl/drawobject.hpp"

void linevertarray::draw() const
{
	if (ibo != 0) {
		this->bind();
		glDrawElements(GL_LINES, isize, GL_UNSIGNED_INT, (void*)0);
		this->unbind();
	} else {
		this->bind();
		glDrawArrays(GL_LINES, 0, size);
		this->unbind();
	}
}

void linestripvertarray::draw() const
{
	if (ibo != 0) {
		this->bind();
		glDrawElements(GL_LINES, isize, GL_UNSIGNED_INT, (void*)0);
		this->unbind();
	} else {
		this->bind();
		glDrawArrays(GL_LINES, 0, size);
		this->unbind();
	}
}

void pointvertarray::draw() const
{
	if (ibo != 0) {
		this->bind();
		glDrawElements(GL_POINTS, isize, GL_UNSIGNED_INT, (void*)0);
		this->unbind();
	} else {
		this->bind();
		glDrawArrays(GL_POINTS, 0, size);
		this->unbind();
	}
}

void trianglevertarray::draw() const
{
	if (ibo != 0) {
		this->bind();
		glDrawElements(GL_TRIANGLES, isize, GL_UNSIGNED_INT, (void*)0);
		this->unbind();
	} else {
		this->bind();
		glDrawArrays(GL_TRIANGLES, 0, size);
		this->unbind();
	}
}

spherevertarray::spherevertarray(const float& cx, const float& cy, const float& cz, const float& radi, const uint32_t& N,
    const uint32_t& M)
    : vertarray(N * (M + 1), nullptr, 6 * N * M, nullptr)
    , N(N)
    , M(M)
    , radi(radi)
{
	vertex vdata[N * (M + 1)];
	for (uint32_t i = 0; i < N; i++) {
		for (uint32_t j = 0; j <= M; j++) {
			float theta = 2.0 * 3.1415 * (float(i)) / (float(N));
			float phi   = 3.1415 * ((float(j)) / (float(M)) - 0.5);

			float sint = std::sin(theta);
			float cost = std::cos(theta);
			float sinp = std::sin(phi);
			float cosp = std::cos(phi);

			float x = cosp * cost;
			float z = -cosp * sint;
			float y = sinp;

			vdata[j * N + i].position[0] = cx + radi * x;
			vdata[j * N + i].position[1] = cy + radi * y;
			vdata[j * N + i].position[2] = cz + radi * z;
			vdata[j * N + i].color[0]    = 1.0f;
			vdata[j * N + i].color[1]    = 1.0f;
			vdata[j * N + i].color[2]    = 1.0f;
			vdata[j * N + i].color[3]    = 1.0f;
			vdata[j * N + i].uv[0]	     = ((float)i) / ((float)N);
			vdata[j * N + i].uv[1]	     = ((float)j) / ((float)M);
			vdata[j * N + i].normal[0]   = x;
			vdata[j * N + i].normal[1]   = y;
			vdata[j * N + i].normal[2]   = z;
		}
	}

	this->setdata(vdata);

	uint32_t idata[6 * N * M];
	for (uint32_t i = 0; i < N; i++) {
		for (uint32_t j = 0; j < M; j++) {
			uint32_t k0 = N * j + i;
			uint32_t k1 = k0 + 1;
			if (i == N - 1)
				k1 = N * j;
			uint32_t k2 = k0 + N;
			uint32_t k3 = k2 + 1;
			if (i == N - 1)
				k3 = N * (j + 1);

			idata[6 * (N * j + i) + 0] = k0;
			idata[6 * (N * j + i) + 1] = k1;
			idata[6 * (N * j + i) + 2] = k2;

			idata[6 * (N * j + i) + 3] = k2;
			idata[6 * (N * j + i) + 4] = k1;
			idata[6 * (N * j + i) + 5] = k3;
		}
	}

	this->setilist(idata);
}

void spherevertarray::draw() const
{
	if (ibo != 0) {
		this->bind();
		glDrawElements(GL_TRIANGLES, isize, GL_UNSIGNED_INT, (void*)0);
		this->unbind();
	} else {
		this->bind();
		glDrawArrays(GL_TRIANGLES, 0, size);
		this->unbind();
	}
}

void spherevertarray::setcenter(const float& x, const float& y, const float& z)
{
	cx = x;
	cy = y;
	cz = z;
}

void spherevertarray::setradi(const float& r)
{
	radi = r;
}

void spherevertarray::update()
{
	for (uint32_t i = 0; i < N; i++) {
		for (uint32_t j = 0; j <= M; j++) {
			float theta = 2.0 * 3.1415 * (float(i)) / (float(N));
			float phi   = 3.1415 * ((float(j)) / (float(M)) - 0.5);

			float sint = std::sin(theta);
			float cost = std::cos(theta);
			float sinp = std::sin(phi);
			float cosp = std::cos(phi);

			float x = cosp * cost;
			float z = -cosp * sint;
			float y = sinp;

			va[j * N + i].position[0] = cx + radi * x;
			va[j * N + i].position[1] = cy + radi * y;
			va[j * N + i].position[2] = cz + radi * z;
			va[j * N + i].normal[0]	  = x;
			va[j * N + i].normal[1]	  = y;
			va[j * N + i].normal[2]	  = z;
		}
	}

	this->vboupdate();
}
