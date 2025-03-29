#include <vector>
#include <cmath>
#include "tgaimage.h"
#include "model.h"
#include "geometry.h"
#include "custommath.h"

class Camera
{
public:
    Camera(Vec3f coordinates);
    ~Camera();
    Vec3f GetPosition();

private:
    Vec3f coords;
};

Camera::Camera(Vec3f coordinates)
{
    this->coords = coordinates;
}

Camera::~Camera()
{
}

Vec3f Camera::GetPosition()
{
    return coords;
}

const TGAColor white = TGAColor(255, 255, 255, 255);
const TGAColor red   = TGAColor(255, 0,   0,   255);
Model *model = NULL;
const int width  = 800;
const int height = 800;
const int depth = 255;
int* zbuffer = NULL;

Vec3f light_dir(0, 0, -1);
Camera camera(Vec3f(0, 0, 300));


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


//void line(int x0, int y0, int x1, int y1, TGAImage &image, TGAColor color) {
//    bool steep = false;
//    if (std::abs(x0-x1)<std::abs(y0-y1)) {
//        std::swap(x0, y0);
//        std::swap(x1, y1);
//        steep = true;
//    }
//    if (x0>x1) {
//        std::swap(x0, x1);
//        std::swap(y0, y1);
//    }
//
//    for (int x=x0; x<=x1; x++) {
//        float t = (x-x0)/(float)(x1-x0);
//        int y = y0*(1.-t) + y1*t;
//        if (steep) {
//            image.set(y, x, color);
//        } else {
//            image.set(x, y, color);
//        }
//    }
//}

void triangle(Vec3i t0, Vec3i t1, Vec3i t2, Vec2i uv0, Vec2i uv1, Vec2i uv2, TGAImage& image, float intensity, int* zbuffer) {
    if (t0.y == t1.y && t0.y == t2.y) return; // i dont care about degenerate triangles
    if (t0.y > t1.y) { std::swap(t0, t1); std::swap(uv0, uv1); }
    if (t0.y > t2.y) { std::swap(t0, t2); std::swap(uv0, uv2); }
    if (t1.y > t2.y) { std::swap(t1, t2); std::swap(uv1, uv2); }
    int total_height = t2.y - t0.y;
    for (int i = 0; i < total_height; i++) {
        bool second_half = i > t1.y - t0.y || t1.y == t0.y;
        int segment_height = second_half ? t2.y - t1.y : t1.y - t0.y;
        float alpha = (float)i / total_height;
        float beta = (float)(i - (second_half ? t1.y - t0.y : 0)) / segment_height; // be careful: with above conditions no division by zero here
        
        /*Vec3f DA((float)(t2.x - t0.x) * alpha, (float)(t2.y - t0.y) * alpha, (float)(t2.z - t0.z) * alpha);
        Vec3i A = t0 + Vec3i(DA.x, DA.y, DA.z);

        Vec3f DB;
        if (second_half)
        {
            DB = Vec3f((float)(t2.x - t1.x) * beta, (float)(t2.y - t1.y) * beta, (float)(t2.z - t1.z) * beta);
        }
        else
        {
            DB = Vec3f((float)(t1.x - t0.x) * beta, (float)(t1.y - t0.y) * beta, (float)(t1.z - t0.z) * beta);
        }
        Vec3i DBi(DB.x, DB.y, DB.z);
        Vec3i B = second_half ? t1 + DBi : t0 + DBi;*/

        Vec3f A = Vec3f(t0.x, t0.y, t0.z) + Vec3f((t2 - t0).x, (t2 - t0).y, (t2 - t0).z) * alpha;
        Vec3f B = second_half ? Vec3f(t1.x, t1.y, t1.z) + Vec3f((t2 - t1).x, (t2 - t1).y, (t2 - t1).z) * beta : Vec3f(t0.x, t0.y, t0.z) + Vec3f((t1 - t0).x, (t1 - t0).y, (t1 - t0).z) * beta;

        Vec2f uvA = Vec2f(uv0.x, uv0.y) + Vec2f((uv2.x - uv0.x), (uv2.y - uv0.y)) * alpha;
        Vec2f uvB = second_half ? Vec2f(uv1.x, uv1.y) + Vec2f((uv2.x - uv1.x), (uv2.y - uv1.y)) * beta : Vec2f(uv0.x, uv0.y) + Vec2f((uv1.x - uv0.x), (uv1.y - uv0.y)) * beta;

        if (A.x > B.x) std::swap(A, B);
        for (int j = A.x; j <= B.x; j++) {
            float phi = B.x == A.x ? 1. : (float)(j - A.x) / (float)(B.x - A.x);

            Vec3f Pshit = Vec3f(A.x, A.y, A.z) + Vec3f((B.x - A.x), (B.y - A.y), (B.z - A.z)) * phi; //Pshit = A + B * phi; 
            Vec3i P(Pshit.x, Pshit.y, Pshit.z);

            Vec2f uvPshit = Vec2f(uvA.x, uvA.y) + Vec2f((uvB.x - uvA.x), (uvB.y - uvA.y)) * phi;
            Vec2i uvP(uvPshit.x, uvPshit.y);

            int idx = P.x + P.y * width;
            
            if (idx >= height*width or idx < 0)
            {
                continue;
            }

            if (zbuffer[idx] < P.z) {
                zbuffer[idx] = P.z;
                TGAColor color = 0 ? TGAColor(255, 0, 0, 1) : model->diffuse(uvP);
                image.set(P.x, P.y, TGAColor(color.r * intensity, color.g * intensity, color.b * intensity, 1));
            }
        }
    }
}



int main(int argc, char** argv) {
    if (2==argc) {
        model = new Model(argv[1]);
    } else {
        model = new Model("obj/african_head.obj");
    }

    zbuffer = new int[width * height];
    for (int i = 0; i < width * height; i++) {
        zbuffer[i] = std::numeric_limits<int>::min();
    }


    Matrix Projection = Matrix::identity(4);
    Matrix ViewPort = viewport(width / 8, height / 8, width * 3 / 4, height * 3 / 4);
    Projection[3][2] = -1.f / camera.GetPosition().z;


    TGAImage image(width, height, TGAImage::RGB);
    for (int i = 0; i < model->nfaces(); i++) {
        std::vector<int> face = model->face(i);
        Vec3i screen_coords[3];
        Vec3f world_coords[3];

        for (int j = 0; j < 3; j++) {
            Vec3f v = model->vert(face[j]);
            Vec3f sup = Matrix::mat2vec(ViewPort * Projection * Matrix::vec2mat(v));
            screen_coords[j] = Vec3i(sup.x, sup.y, sup.z); //Vec3i((v.x + 1.) * width / 2., (v.y + 1.) * height / 2., (v.z + 1.) * depth / 2.);
            world_coords[j] = v;
        }

        Vec3f n = (world_coords[2] - world_coords[0]) ^ (world_coords[1] - world_coords[0]);
        n.normalize();
        float intensity = n * light_dir;

        if (intensity > 0) {
            Vec2i uv[3];
            for (int k = 0; k < 3; k++) {
                uv[k] = model->uv(i, k);
            }
            triangle(screen_coords[0], screen_coords[1], screen_coords[2], uv[0], uv[1], uv[2], image, intensity, zbuffer);
        }
    }

    TGAImage zbuffer_image(width, height, TGAImage::RGB);
    for (int i = 0; i < width * height; i++)
    {
        int a = zbuffer[i] == std::numeric_limits<int>::min() ? 0 : zbuffer[i];
        zbuffer_image.set(i % width, i / width, TGAColor(255 * a, 255 * a, 255 * a, 255));
    }

    image.flip_vertically(); // i want to have the origin at the left bottom corner of the image
    image.write_tga_file("output.tga");

    zbuffer_image.flip_vertically(); // i want to have the origin at the left bottom corner of the image
    zbuffer_image.write_tga_file("zbuffer_image_output.tga");
    delete model;
    return 0;
}

