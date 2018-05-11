using System;
using MathNet.Numerics.LinearAlgebra;
using MathNet.Numerics.LinearAlgebra.Double;
using MathNet.Numerics.Distributions;

namespace QRDecompositionClass
{

    class Program
    {
        static void Main(string[] args)
        {
            ///Matrix<double> A = DenseMatrix.OfArray(new double[,] { { 1, -2, -1, 1 }, { 1, -8, -2, -3 }, { 2, 2, -1, 7 }, { 1, 1, 2, 1 } });
            Matrix<double> B = DenseMatrix.OfArray(new double[,] { { 1 }, { -2 }, { 7 }, { 1 } });
            Matrix<double> A = DenseMatrix.CreateRandom(4, 4, Normal.WithMeanVariance(5.0, 5.0));
            Console.WriteLine(A);
            QRDecomposition qRDecomposition = new QRDecomposition(A);
            Console.WriteLine(qRDecomposition.Q * qRDecomposition.R);
            Console.WriteLine(A*qRDecomposition.Solve(B));
        }
    }
}