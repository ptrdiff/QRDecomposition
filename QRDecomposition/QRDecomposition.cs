using System;
using MathNet.Numerics.LinearAlgebra;
using MathNet.Numerics.LinearAlgebra.Double;

namespace QRDecompositionClass
{
    class QRDecomposition
    {
        private Matrix<double> _Q;
        private Matrix<double> _R;

        public Matrix<double> Q { get => _Q; }
        public Matrix<double> R { get => _R; }

        public QRDecomposition(Matrix<double> matrix)
        {
            Vector<double> vectorU = new DenseVector(matrix.RowCount);
            Vector<double> vectorV = new DenseVector(matrix.RowCount);

            Matrix<double> P = DenseMatrix.CreateIdentity(matrix.RowCount);
            Matrix<double> Q = DenseMatrix.CreateIdentity(matrix.RowCount);
            Matrix<double> R = matrix;

            double mag, alpha;

            for (int i = 0; i < matrix.ColumnCount; ++i)
            {
                vectorU.Clear();
                vectorV.Clear();
                mag = 0.0;
                for (int j = i; j < matrix.RowCount; ++j)
                {
                    vectorU[j] = R[j, i];
                    mag += vectorU[j] * vectorU[j];
                }
                mag = Math.Sqrt(mag);

                alpha = vectorU[i] < 0 ? mag : -mag;

                mag = 0.0;
                for (int j = i; j < matrix.RowCount; ++j)
                {
                    vectorV[j] = j == i ? vectorU[j] + alpha : vectorU[j];
                    mag += vectorV[j] * vectorV[j];
                }
                mag = Math.Sqrt(mag);

                if (mag < 1e-16) continue;

                for (int j = i; j < matrix.RowCount; ++j) vectorV[j] /= mag;

                P = DenseMatrix.CreateIdentity(matrix.RowCount) - (2 * (vectorV.ToColumnMatrix() * vectorV.ToRowMatrix()));
                R = P * R;
                Q = Q * P;
            }

            _Q = Q;
            _R = R;
        }

        public Matrix<double> Solve(Matrix<double> B)
        {
            Matrix<double> Y = _Q.Transpose() * B;
            Matrix<double> X = new DenseMatrix(B.RowCount,1);
            for (int i = _R.RowCount - 1; i >= 0; --i)
            {
                X[i,0] = (Y[i, 0] - _R.Row(i) * X.Column(0)) / _R[i, i];
            }
            return X;
        }
    }
}
