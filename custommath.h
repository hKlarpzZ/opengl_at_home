#include <cmath>

const int DEFAULT_ALLOC = 4;

class Matrix {
    std::vector<std::vector<float> > m;
    int rows, cols;
public:
    Matrix(int r = DEFAULT_ALLOC, int c = DEFAULT_ALLOC);
    inline int nrows();
    inline int ncols();

    static Matrix identity(int dimensions);
    std::vector<float>& operator[](const int i);
    Matrix operator*(const Matrix& a);
    Matrix transpose();
    Matrix inverse();
    
    static Vec3f mat2vec(Matrix m);
    static Matrix vec2mat(Vec3f v);

    friend std::ostream& operator<<(std::ostream& s, Matrix& m);
};
