#include <vector>
#include <cmath>
#include "tgaimage.h"
#include "model.h"
#include "geometry.h"
#include "custommath.h"
#include <iostream>

Model* model = NULL;
int* zbuffer = NULL;

const int width = 800;
const int height = 800;
const int depth = 255;

class Camera
{
public:
	Vec3f position;
	Vec3f light_dir;
	Vec3f eye;
	Vec3f center;

	Matrix ModelView;
	Matrix Projection;
	Matrix ViewPort;

	Matrix lookat(Vec3f eye, Vec3f center, Vec3f up) {
		Vec3f z = (eye - center).normalize();
		Vec3f x = (up ^ z).normalize();
		Vec3f y = (z ^ x).normalize();
		Matrix res = Matrix::identity(4);
		for (int i = 0; i < 3; i++) {
			res[0][i] = x.raw[i];
			res[1][i] = y.raw[i];
			res[2][i] = z.raw[i];
			res[i][3] = -center.raw[i];
		}
		return res;
	}

	Matrix viewport(int x, int y, int w, int h) {
		Matrix m = Matrix::identity(4);
		m[0][3] = x + w / 2.f;
		m[1][3] = y + h / 2.f;
		m[2][3] = depth / 2.f;

		m[0][0] = w / 2.f;
		m[1][1] = h / 2.f;
		m[2][2] = depth / 2.f;
		return m;
	}

	Camera() {
		this->position = Vec3f(0, 0, 3);
		this->light_dir = Vec3f(0, 0, 1);
		this->eye = Vec3f(1, 1, 3);
		this->center = Vec3f(0, 0, 0);

		this->ModelView = lookat(eye, center, Vec3f(0, 1, 0));
		this->Projection = Matrix::identity(4);
		this->ViewPort = viewport(width / 8, height / 8, width * 3 / 4, height * 3 / 4);
		this->Projection[3][2] = -1.f / this->position.z;
	}
};

Camera* camera = NULL;

struct IShader {
	//Vertex Shader
	virtual Vec3f vertex(int iface, int ivertex) = 0;
	//Fragment/Pixel (DirectX) Shader
	virtual bool fragment(Vec3f bar, TGAColor& color) = 0;
};

struct GouraudShader : public IShader {
	Vec3f varying_intensity; // written by vertex shader, read by fragment shader
	
	Vec2i varying_uv[3]; // triangle uv coordinates, written by the vertex shader, read by the fragment shader
	Vec3f varying_tri[3]; // triangle coordinates (clip space), written by VS, read by FS

	Vec3f vertex(int iface, int ivertex) override {
		// get diffuse lighting intensity
		varying_intensity.raw[ivertex] = std::max(0.f, model->normal(iface, ivertex) * camera->light_dir); //ambient occlusion?
		varying_uv[ivertex] = model->uv(iface, ivertex);
		// read the vertex from .obj file and transform it to screen coordinates
		Vec3f gl_Vertex = Matrix::mat2vec(camera->ViewPort * camera->Projection * camera->ModelView * Matrix::vec2mat(model->vert(model->face(iface)[ivertex])));
		
		varying_tri[ivertex] = gl_Vertex;

		return gl_Vertex;
	}
	bool fragment(Vec3f bar, TGAColor& color) override {
		// interpolate intensity for the current pixel
		Vec2f uv;
		uv.x = varying_uv[0].x * bar.x + varying_uv[1].x * bar.y + varying_uv[2].x * bar.z;
		uv.y = varying_uv[0].y * bar.x + varying_uv[1].y * bar.y + varying_uv[2].y * bar.z;
		Vec2i uv_int(uv.x, uv.y);
		color = model->diffuse(uv_int);

		float intensity = varying_intensity * bar;
		color.r *= intensity;
		color.g *= intensity;
		color.b *= intensity;

		return false; // no, we do not discard this pixel
	}
};

struct PhongShader : public IShader {
	Vec3f varying_intensity; // written by vertex shader, read by fragment shader

	Vec2i varying_uv[3]; // triangle uv coordinates, written by the vertex shader, read by the fragment shader
	Vec3f varying_tri[3]; // triangle coordinates (clip space), written by VS, read by FS

	Vec3f varying_spec;
	Vec3f varying_diffuse;

	Vec3f vertex(int iface, int ivertex) override {
		// get diffuse lighting intensity
		varying_intensity.raw[ivertex] = std::max(0.f, model->normal(iface, ivertex) * camera->light_dir); //ambient occlusion?
		Vec3f normal_vec = model->normal(iface, ivertex);
		Vec3f reflection_vec = Vec3f(2 * (normal_vec * camera->light_dir) * normal_vec.x - camera->light_dir.x,\
			2 * (normal_vec * camera->light_dir) * normal_vec.y - camera->light_dir.y,\
			2 * (normal_vec * camera->light_dir) * normal_vec.z - camera->light_dir.z);
		varying_spec.raw[ivertex] = std::max(0.f, reflection_vec * Vec3f(camera->position - model->vert(ivertex)).normalize());
		varying_diffuse.raw[ivertex] = std::max(0.f, normal_vec * camera->light_dir);
		varying_uv[ivertex] = model->uv(iface, ivertex);
		// read the vertex from .obj file and transform it to screen coordinates
		Vec3f gl_Vertex = Matrix::mat2vec(camera->ViewPort * camera->Projection * camera->ModelView * Matrix::vec2mat(model->vert(model->face(iface)[ivertex])));

		varying_tri[ivertex] = gl_Vertex;

		return gl_Vertex;
	}
	bool fragment(Vec3f bar, TGAColor& color) override {
		// interpolate intensity for the current pixel
		Vec2f uv;
		uv.x = varying_uv[0].x * bar.x + varying_uv[1].x * bar.y + varying_uv[2].x * bar.z;
		uv.y = varying_uv[0].y * bar.x + varying_uv[1].y * bar.y + varying_uv[2].y * bar.z;
		Vec2i uv_int(uv.x, uv.y);
		color = model->diffuse(uv_int);

		//float intensity = varying_intensity * bar;
		float specular = varying_spec * bar;
		float diffuse = varying_diffuse * bar;
		color.r *= diffuse + specular * .5 + diffuse * .01;
		color.g *= diffuse + specular * .5 + diffuse * .01;
		color.b *= diffuse + specular * .5 + diffuse * .01;

		return false; // no, we do not discard this pixel
	}
};

inline Vec3f cross(const Vec3f& a, const Vec3f& b) {
	return Vec3f(
		a.y * b.z - a.z * b.y,
		a.z * b.x - a.x * b.z,
		a.x * b.y - a.y * b.x
	);
}

Vec3f barycentric(Vec2i A, Vec2i B, Vec2i C, Vec2i P) {
	Vec3f s[2];
	for (int i = 2; i--; ) {
		s[i].raw[0] = C.raw[i] - A.raw[i];
		s[i].raw[1] = B.raw[i] - A.raw[i];
		s[i].raw[2] = A.raw[i] - P.raw[i];
	}
	Vec3f u = cross(s[0], s[1]);
	if (std::abs(u.raw[2]) < 1e-2) return Vec3f(-1, 1, 1);
	return Vec3f(1.f - (u.x + u.y) / u.z, u.y / u.z, u.x / u.z);
}

void triangle(Vec3f ScreenCords[3], IShader& shader, TGAImage& image, int* zbuffer) {
	// finding bounding box
	Vec2i bboxmin(image.get_width() - 1, image.get_height() - 1);
	Vec2i bboxmax(0, 0);
	Vec2i clamp(image.get_width() - 1, image.get_height() - 1);

	Vec2i pts[3] = { Vec2i(ScreenCords[0].x, ScreenCords[0].y), Vec2i(ScreenCords[1].x, ScreenCords[1].y), Vec2i(ScreenCords[2].x, ScreenCords[2].y) };

	for (int i = 0; i < 3; i++) {
		for (int j = 0; j < 2; j++) {
			bboxmin.raw[j] = std::max(0, std::min(bboxmin.raw[j], pts[i].raw[j]));
			bboxmax.raw[j] = std::min(clamp.raw[j], std::max(bboxmax.raw[j], pts[i].raw[j]));
		}
	}

	// drawing pixels in bounding box
	Vec2i P;
	for (P.x = bboxmin.x; P.x <= bboxmax.x; P.x++) {
		for (P.y = bboxmin.y; P.y <= bboxmax.y; P.y++) {
			Vec3f bc = barycentric(pts[0], pts[1], pts[2], P);
			if (bc.x < 0 || bc.y < 0 || bc.z < 0) continue;

			// »нтерпол€ци€ глубины
			float z = ScreenCords[0].z * bc.x + ScreenCords[1].z * bc.y + ScreenCords[2].z * bc.z;

			int idx = P.x + P.y * image.get_width();
			if (zbuffer[idx] < z) {
				zbuffer[idx] = z;

				TGAColor color(255, 255, 255, 255);
				shader.fragment(bc, color);
				image.set(P.x, P.y, color);
			}
		}
	}
}


int main(int argc, char** argv) {
	if (2 == argc) {
		model = new Model(argv[1]);
	}
	else {
		model = new Model("obj/african_head.obj");
	}

	zbuffer = new int[width * height];
	for (int i = 0; i < width * height; i++) {
		zbuffer[i] = std::numeric_limits<int>::min();
	}

	camera = new Camera();

	TGAImage image(width, height, TGAImage::RGB);

	//std::cout << "Total faces " << model->nfaces() << std::endl;

	PhongShader shader;
	//GouraudShader shader;
	for (int iface = 0; iface < model->nfaces(); iface++) {
		//std::cout << "[F] Found face with id [" << iface << "] and data " << model->face(iface)[0] << " " << model->face(iface)[1] << " " << model->face(iface)[2] << ":" << std::endl;
		//Vec3f screen_coords[3];
		for (int ivertex = 0; ivertex < 3; ivertex++) {
			//screen_coords[ivertex] = shader.vertex(iface, ivertex);
			shader.vertex(iface, ivertex);
			/*std::cout << " - Got vertex (" << model->vert(ivertex).x << ", " << model->vert(ivertex).y << ", "\
				<< model->vert(ivertex).z << ") to screen coords:" << std::endl;
			std::cout << "   [V] (" << screen_coords[ivertex].x << ", " << screen_coords[ivertex].y << ", " << screen_coords[ivertex].z << ")" << std::endl;*/
		}
		//std::cout << std::endl;
		triangle(shader.varying_tri, shader, image, zbuffer);
	}

	image.flip_vertically(); // i want to have the origin at the left bottom corner of the image
	image.write_tga_file("output.tga");


	TGAImage zbuffer_image(width, height, TGAImage::RGB);
	for (int i = 0; i < width * height; i++)
	{
		int a = zbuffer[i] == std::numeric_limits<int>::min() ? 0 : zbuffer[i];
		zbuffer_image.set(i % width, i / width, TGAColor(255 * a, 255 * a, 255 * a, 255));
	}

	zbuffer_image.flip_vertically(); // i want to have the origin at the left bottom corner of the image
	zbuffer_image.write_tga_file("zbuffer_image_output.tga");


	delete model;
	delete camera;

	return 0;
}