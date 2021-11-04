#version 460 core

in vec4 position;
in vec4 color;
in vec2 uv;
in vec3 normal;
in uint type;

out vec4 vcolor;
out vec3 worldpos;
out vec3 vnormal;
out vec2 vuv;
flat out uint vtype;
out vec4 Lc;

uniform mat4 pers;
uniform mat4 euclid;
uniform mat4 extraeuclid;

uniform mat4 Lpers;
uniform mat4 Leuclid;

void main()
{
	gl_Position = pers * euclid * extraeuclid * position;

	vcolor	 = color;
	worldpos = (extraeuclid * position).xyz;
	vnormal	 = (extraeuclid * vec4(normal, 0.0)).xyz;
	vuv	 = uv;
	vtype	 = type;

	Lc = Lpers * Leuclid * extraeuclid * position;
}
