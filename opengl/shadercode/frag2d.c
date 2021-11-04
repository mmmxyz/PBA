#version 460 core

out vec4 fragment;

in vec4 vcolor;

void main()
{
	fragment = vcolor;
}
