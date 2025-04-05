//#include <vector>
//#include <cmath>
//#include "tgaimage.h"
//#include "model.h"
//#include "geometry.h"
//#include "custommath.h"
//
//class Camera
//{
//public:
//    Camera(Vec3f coordinates);
//    ~Camera();
//    Vec3f GetPosition();
//
//private:
//    Vec3f coords;
//};
//
//Camera::Camera(Vec3f coordinates)
//{
//    this->coords = coordinates;
//}
//
//Camera::~Camera()
//{
//}
//
//Vec3f Camera::GetPosition()
//{
//    return coords;
//}
//
//const TGAColor white = TGAColor(255, 255, 255, 255);
//const TGAColor red   = TGAColor(255, 0,   0,   255);
//Model *model = NULL;
//const int width  = 800;
//const int height = 800;
//const int depth = 255;
//int* zbuffer = NULL;
//
//Vec3f light_dir(0, 0, -1);
//Camera camera(Vec3f(0, 0, 3));
//Vec3f eye(1, 1, 3);
//Vec3f center(0, 0, 0);
//
//
//Matrix viewport(int x, int y, int w, int h) {
//    Matrix m = Matrix::identity(4);
//    m[0][3] = x + w / 2.f;
//    m[1][3] = y + h / 2.f;
//    m[2][3] = depth / 2.f;
//
//    m[0][0] = w / 2.f;
//    m[1][1] = h / 2.f;
//    m[2][2] = depth / 2.f;
//    return m;
//}
//
//Matrix lookat(Vec3f eye, Vec3f center, Vec3f up) {
//    Vec3f z = (eye - center).normalize();
//    Vec3f x = (up ^ z).normalize();
//    Vec3f y = (z ^ x).normalize();
//    Matrix res = Matrix::identity(4);
//    for (int i = 0; i < 3; i++) {
//        res[0][i] = x.raw[i];
//        res[1][i] = y.raw[i];
//        res[2][i] = z.raw[i];
//        res[i][3] = -center.raw[i];
//    }
//    return res;
//}
//
////void line(int x0, int y0, int x1, int y1, TGAImage &image, TGAColor color) {
////    bool steep = false;
////    if (std::abs(x0-x1)<std::abs(y0-y1)) {
////        std::swap(x0, y0);
////        std::swap(x1, y1);
////        steep = true;
////    }
////    if (x0>x1) {
////        std::swap(x0, x1);
////        std::swap(y0, y1);
////    }
////
////    for (int x=x0; x<=x1; x++) {
////        float t = (x-x0)/(float)(x1-x0);
////        int y = y0*(1.-t) + y1*t;
////        if (steep) {
////            image.set(y, x, color);
////        } else {
////            image.set(x, y, color);
////        }
////    }
////}
//
//inline Vec3f cross(const Vec3f& a, const Vec3f& b) {
//    return Vec3f(
//        a.y * b.z - a.z * b.y,
//        a.z * b.x - a.x * b.z,
//        a.x * b.y - a.y * b.x
//    );
//}
//
//Vec3f barycentric(Vec2i A, Vec2i B, Vec2i C, Vec2i P) {
//    Vec3f s[2];
//    for (int i = 2; i--; ) {
//        s[i].raw[0] = C.raw[i] - A.raw[i];
//        s[i].raw[1] = B.raw[i] - A.raw[i];
//        s[i].raw[2] = A.raw[i] - P.raw[i];
//    }
//    Vec3f u = cross(s[0], s[1]);
//    if (std::abs(u.raw[2]) < 1e-2) return Vec3f(-1, 1, 1);
//    return Vec3f(1.f - (u.x + u.y) / u.z, u.y / u.z, u.x / u.z);
//}
//
//void triangle(Vec3i t0, Vec3i t1, Vec3i t2, Vec2i uv0, Vec2i uv1, Vec2i uv2, TGAImage& image, float intensity, int* zbuffer) {
//    // Находим ограничивающий прямоугольник
//    Vec2i bboxmin(image.get_width() - 1, image.get_height() - 1);
//    Vec2i bboxmax(0, 0);
//    Vec2i clamp(image.get_width() - 1, image.get_height() - 1);
//
//    Vec2i pts[3] = { Vec2i(t0.x, t0.y), Vec2i(t1.x, t1.y), Vec2i(t2.x, t2.y) };
//
//    for (int i = 0; i < 3; i++) {
//        for (int j = 0; j < 2; j++) {
//            bboxmin.raw[j] = std::max(0, std::min(bboxmin.raw[j], pts[i].raw[j]));
//            bboxmax.raw[j] = std::min(clamp.raw[j], std::max(bboxmax.raw[j], pts[i].raw[j]));
//        }
//    }
//
//    // Перебираем все пиксели в ограничивающем прямоугольнике
//    Vec2i P;
//    for (P.x = bboxmin.x; P.x <= bboxmax.x; P.x++) {
//        for (P.y = bboxmin.y; P.y <= bboxmax.y; P.y++) {
//            Vec3f bc = barycentric(pts[0], pts[1], pts[2], P);
//            if (bc.x < 0 || bc.y < 0 || bc.z < 0) continue;
//
//            // Интерполяция глубины
//            float z = t0.z * bc.x + t1.z * bc.y + t2.z * bc.z;
//
//            int idx = P.x + P.y * image.get_width();
//            if (zbuffer[idx] < z) {
//                zbuffer[idx] = z;
//
//                // Интерполяция UV-координат
//                float u = uv0.x * bc.x + uv1.x * bc.y + uv2.x * bc.z;
//                float v = uv0.y * bc.x + uv1.y * bc.y + uv2.y * bc.z;
//
//                TGAColor color = model->diffuse(Vec2i(u, v));
//                color.r *= intensity;
//                color.g *= intensity;
//                color.b *= intensity;
//                image.set(P.x, P.y, color);
//            }
//        }
//    }
//}
//
//int main(int argc, char** argv) {
//    if (2==argc) {
//        model = new Model(argv[1]);
//    } else {
//        model = new Model("obj/cat.obj");
//    }
//
//    zbuffer = new int[width * height];
//    for (int i = 0; i < width * height; i++) {
//        zbuffer[i] = std::numeric_limits<int>::min();
//    }
//
//    Matrix ModelView = lookat(eye, center, Vec3f(0, 1, 0));
//    Matrix Projection = Matrix::identity(4);
//    Matrix ViewPort = viewport(width / 8, height / 8, width * 3 / 4, height * 3 / 4);
//    Projection[3][2] = -1.f / camera.GetPosition().z;
//
//    Matrix z = (ViewPort * Projection * ModelView);
//
//    TGAImage image(width, height, TGAImage::RGB);
//    for (int i = 0; i < model->nfaces(); i++) {
//        std::vector<int> face = model->face(i);
//        Vec3i screen_coords[3];
//        Vec3f world_coords[3];
//
//        for (int j = 0; j < 3; j++) {
//            Vec3f v = model->vert(face[j]);
//            Matrix gl_Vertex = ViewPort * Projection * ModelView * Matrix::vec2mat(v);
//            Vec3f vert_screen_coords = Matrix::mat2vec(gl_Vertex);
//            screen_coords[j] = Vec3i(vert_screen_coords.x, vert_screen_coords.y, vert_screen_coords.z);
//            world_coords[j] = v;
//        }
//
//        Vec3f n = (world_coords[2] - world_coords[0]) ^ (world_coords[1] - world_coords[0]);
//        n.normalize();
//        float intensity = n * light_dir;
//
//        if (intensity > 0) {
//            Vec2i uv[3];
//            for (int k = 0; k < 3; k++) {
//                uv[k] = model->uv(i, k);
//            }
//            triangle(screen_coords[0], screen_coords[1], screen_coords[2], uv[0], uv[1], uv[2], image, intensity, zbuffer);
//        }
//    }
//
//    TGAImage zbuffer_image(width, height, TGAImage::RGB);
//    for (int i = 0; i < width * height; i++)
//    {
//        int a = zbuffer[i] == std::numeric_limits<int>::min() ? 0 : zbuffer[i];
//        zbuffer_image.set(i % width, i / width, TGAColor(255 * a, 255 * a, 255 * a, 255));
//    }
//
//    image.flip_vertically(); // i want to have the origin at the left bottom corner of the image
//    image.write_tga_file("output.tga");
//
//    zbuffer_image.flip_vertically(); // i want to have the origin at the left bottom corner of the image
//    zbuffer_image.write_tga_file("zbuffer_image_output.tga");
//    delete model;
//    return 0;
//}
//
