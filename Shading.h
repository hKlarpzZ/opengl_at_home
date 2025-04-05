#pragma once
#include <vector>
#include <cmath>
#include "tgaimage.h"
#include "model.h"
#include "geometry.h"
#include "custommath.h"
#include <iostream>

struct IShader {
	virtual ~IShader();
	//Vertex Shader
	virtual Vec3i vertex(int iface, int ivertex) = 0;
	//Fragment/Pixel (DirectX) Shader
	virtual bool fragment(Vec3f bar, TGAColor& color) = 0;
};

struct GouraudShader : public IShader {
	Vec3f varying_intensity; // written by vertex shader, read by fragment shader
	virtual Vec4f vertex(int iface, int ivertex) {
		// get diffuse lighting intensity
		varying_intensity[ivertex] = std::max(0.f, model->normal(iface, ivertex) * light_dir);
		// read the vertex from .obj file
		Vec4f gl_Vertex = embed<4>(model->vert(iface, ivertex));
		// transform it to screen coordinates
		return Viewport * Projection * ModelView * gl_Vertex;
	}
	virtual bool fragment(Vec3f bar, TGAColor& color) {
		// interpolate intensity for the current pixel
		float intensity = varying_intensity * bar;
		color = TGAColor(255 * intensity, 255 * intensity * intensity, 255, 255 * intensity;
		return false; // no, we do not discard this pixel
	}
};