#version 460 core

out vec4 fragment;

in vec4 vcolor;
in vec3 worldpos;
in vec3 vnormal;
in vec2 vuv;
flat in uint vtype;
in vec4 Lc;

uniform vec3 Lpos;
uniform float Lint;
uniform float Lamb;

uniform sampler2DShadow shadowtex;
uniform sampler2D colortexture;

void main()
{
	fragment = vcolor;

	if (vtype == 1) {
		vec3 L	    = normalize(Lpos - worldpos);
		float dist  = length(Lpos - worldpos);
		float dist2 = max(dist * dist, 0.000001);
		vec3 N	    = normalize(vnormal);
		//float Light = Lint * max(abs(dot(L, N)) / (dist2), 0.0001);
		float Light = Lint * max(dot(L, N) / (dist2), 0.0001);

		vec3 LightC  = Lc.xyz / Lc.w;
		float shadow = 1.0;
		if (texture(shadowtex, vec3(0.5 * LightC.xy + vec2(0.5, 0.5), 0.5 * LightC.z + 0.5 - 0.00001)).r == 0.0)
			shadow = 0.3;

		fragment = vec4(vcolor.xyz * (shadow * Light + Lamb), vcolor.w);
	} else if (vtype == 2) {

		vec3 L	    = normalize(Lpos - worldpos);
		float dist  = length(Lpos - worldpos);
		float dist2 = max(dist * dist, 0.000001);
		vec3 N	    = normalize(vnormal);
		//float Light = Lint * max(abs(dot(L, N)) / (dist2), 0.0001);
		float Light = Lint * max(dot(L, N) / (dist2), 0.0001);

		vec3 LightC  = Lc.xyz / Lc.w;
		float shadow = 1.0;
		if (texture(shadowtex, vec3(0.5 * LightC.xy + vec2(0.5, 0.5), 0.5 * LightC.z + 0.5 - 0.00001)).r == 0.0)
			shadow = 0.3;

		vec4 texturecolor = vcolor * texture(colortexture, vuv);

		fragment = vec4(texturecolor.xyz * (shadow * Light + Lamb), texturecolor.w);
	}
}
