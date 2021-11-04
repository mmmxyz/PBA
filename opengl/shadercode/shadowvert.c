#version 460 core

in vec4 position;
in vec4 color;
in vec2 uv;
in vec3 normal;
in uint type;

out vec4 vcolor;

uniform mat4 pers;
uniform mat4 euclid;

void main()
{
		gl_Position = pers * euclid * position;
}
