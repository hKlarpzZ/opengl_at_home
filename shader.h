#include "tgaimage.h"
#include "model.h"
#include "geometry.h"
#include "custommath.h"

class IShader {
public:
    virtual ~IShader() = default;
    virtual Vec3f vertex(int iface, int nthvert, Matrix& gl_Vertex) = 0;
    virtual TGAColor fragment(Vec3f bar, float intensity) = 0;
};

class PhongShader : public IShader {
public:
    Matrix ModelView;
    Matrix Projection;
    Matrix ViewPort;
    Vec3f light_dir;
    Vec3f eye;
    Model* model;

    Matrix varying_tri;  // Координаты вершин
    Matrix varying_nrm;  // Нормали вершин
    Matrix varying_uv;  // Добавляем матрицу для UV-координат

    Vec3f vertex(int iface, int nthvert, Matrix& gl_Vertex) override;
    TGAColor fragment(Vec3f bar, float intensity) override;

    //Vec3f vertex(int iface, int nthvert, Matrix& gl_Vertex) override {
    //    Vec3f v = model->vert(iface);
    //    gl_Vertex = ViewPort * Projection * ModelView * Matrix::vec2mat(v);
    //    Matrix norm = Matrix::vec2mat(model->normal(iface, nthvert));
    //    varying_nrm.set_col(nthvert, norm);
    //    Vec2i uv_int = model->uv(iface, nthvert);
    //    Vec2f uv(uv_int.x, uv_int.y);
    //    Matrix uvmatrix = Matrix::vec2mat(uv);
    //    varying_uv.set_col(nthvert, uvmatrix); // Сохраняем UV
    //    return Matrix::mat2vec(gl_Vertex);
    //}

    //TGAColor fragment(Vec3f bar, float intensity) override {
    //    Vec3f n = Matrix::mat2vec(varying_nrm * Matrix::vec2mat(bar)).normalize();
    //    Vec3f l = light_dir.normalize();
    //    Vec3f r = (n * 2.0f * (n * l) - l).normalize();
    //    Vec3f v = (eye - Matrix::mat2vec(varying_tri * Matrix::vec2mat(bar))).normalize();

    //    float ambient = 0.1f;
    //    float diffuse = std::max(0.0f, n * l);
    //    float specular = pow(std::max(0.0f, r * v), 50);

    //    // Интерполируем UV
    //    Vec3f uv_dump = Matrix::mat2vec(varying_uv * Matrix::vec2mat(bar));
    //    Vec2i uv(uv_dump.x, uv_dump.y);
    //    TGAColor color = model->diffuse(uv);
    //    TGAColor color = model->diffuse(uv);
    //    color.r *= (ambient + diffuse + 0.5f * specular);
    //    color.g *= (ambient + diffuse + 0.5f * specular);
    //    color.b *= (ambient + diffuse + 0.5f * specular);
    //    return color;
    //}
};