#version 460 core

in vec4 position;
in vec3 normal;

uniform mat4 pers;
uniform mat4 euclid;

void main()
{
	vec3 pos    = (position.xyz) / position.w;
	vec4 hpos   = vec4(pos + normal * 0.05, 1.0);
	gl_Position = pers * euclid * hpos;
}
