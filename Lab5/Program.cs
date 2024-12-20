using System;
using System.Collections.Generic;
using System.ComponentModel;
using System.Numerics;
using System.Reflection;
using System.Reflection.Metadata;
using System.Security.Cryptography.X509Certificates;
using System.Text;
using System.Xml;
using static System.Runtime.InteropServices.JavaScript.JSType;

public class Program
{
    public static void Main()
    {
        Program program = new Program();
    }
    
    #region Level 1
    public long Task_1_1(int n, int k)
    {
        long answer = 0;

        // code here
    static int Factorial(int n) {
        int fc = 1;
        for (int i = 2; i <= n; i++)
            fc *= i;
        return fc;
        
    }
    static int Combinations(int n, int k) {
        int fc = Factorial(n) / (Factorial(k) * Factorial(n-k));
        return fc;

    }
        // create and use Combinations(n, k);
        // create and use Factorial(n);
        if (n < 1 && Math.Abs(k) == 1) return 0;
        answer = Combinations(n, k);

        // end

        return answer;
    }

    public int Task_1_2(double[] first, double[] second)
    {
        int answer = 0;

        // code here
        double a1 = first[0], b1 = first[1], c1 = first[2];
        double a2 = second[0], b2 = second[1], c2 = second[2];
        if (a1 >= b1 + c1 || b1 >= a1 + c1 || c1 >= a1 + b1) return -1;
        if (a2 >= b2 + c2 || b2 >= a2 + c2 || c2 >= a2 + b2) return -1;
        static double GeronArea(double a, double b, double c) {
            double p = (a + b + c) / 2;
            return Math.Sqrt(p*(p-a)*(p-b)*(p-c));
        }
        // create and use GeronArea(a, b, c);
        double area = GeronArea(a1, b1, b1) - GeronArea(a2, b2, c2);
        switch (area) {
            case > 0: return 1;
            case < 0: return 2;
            default: return 0;
        }

        // end

        // first = 1, second = 2, equal = 0, error = -1
        return answer;
    }
    public static double GetDistance(double v, double a, double t)
        {
            return Math.Round(v * t + a * t * t / 2, 3);
        }
    public int Task_1_3a(double v1, double a1, double v2, double a2, int time)
    {
        int answer = 0;

        // code here
        
        // create and use GetDistance(v, a, t); t - hours
        switch (GetDistance(v1, a1, time) - GetDistance(v2, a2, time)){
            case > 0: return 1;
            case < 0: return 2;
            default: return 0;
        }

        // end

        // first = 1, second = 2, equal = 0
    }

    public int Task_1_3b(double v1, double a1, double v2, double a2)
    {
        int answer = 0;

        // code here

        // use GetDistance(v, a, t); t - hours
        int t = 1;
        while (GetDistance(a1, v1, t) > GetDistance(v2, a2, t)) 
            t++;
        answer = t;
        // end

        return answer;
    }
    #endregion

    #region Level 2
    void FindMaxIndex(int[,] matrix, out int row, out int column) {
            row = 0; column = 0;
            for (int i = 0; i < matrix.GetLength(0); i++) 
                for (int j = 0; j < matrix.GetLength(1); j++) 
                    if (matrix[i, j] > matrix[row, column]) { row = i; column = j; }
        }      
    void FindMinIndex(int[,] matrix, out int row, out int column) {
            row = 0; column = 0;
            for (int i = 0; i < matrix.GetLength(0); i++) 
                for (int j = 0; j < matrix.GetLength(1); j++) 
                    if (matrix[i, j] < matrix[row, column]) { row = i; column = j; }
        }      
    public void Task_2_1(int[,] A, int[,] B)
    {
        // code here
          
        
        // create and use FindMaxIndex(matrix, out row, out column);
        FindMaxIndex(A, out int rowA, out int columnA);
        FindMaxIndex(B, out int rowB, out int columnB);
        (A[rowA, columnA], B[rowB, columnB]) = (B[rowB, columnB], A[rowA, columnA]);
        // end
    }

    public void Task_2_2(double[] A, double[] B)
    {
        // code here

        // create and use FindMaxIndex(array);
        // only 1 array has to be changed!

        // end
    }
        void RemoveRow (ref int[,] matrix, int row) {
            int n = matrix.GetLength(0), m = matrix.GetLength(1);
            var newMatrix = new int[n-1, m];
            for (int i = 0; i < n-1; i++) 
                if (i < row) 
                    for (int j = 0; j < m; j++)
                        newMatrix[i, j] = matrix[i, j];
                else
                    for (int j = 0; j < m; j++) 
                        newMatrix[i, j] = matrix[i+1, j];
            matrix = newMatrix;
        }
    public void Task_2_3(ref int[,] B, ref int[,] C)
    {
        // code here
        int FindDiagonalMaxIndex(int[,] matrix) {
            int index = 0;
            for (int i = 1; i < matrix.GetLength(0); i++)
                if (matrix[i, i] > matrix[index, index]) index = i;
            return index;
        }
        
        //  create and use method FindDiagonalMaxIndex(matrix);
        int rowB = FindDiagonalMaxIndex(B), rowC = FindDiagonalMaxIndex(C);
        
        RemoveRow(ref B, rowB);
        RemoveRow(ref C, rowC);

        // end
    }

    public void Task_2_4(int[,] A, int[,] B)
    {
        // code here

        //  create and use method FindDiagonalMaxIndex(matrix); like in Task_2_3

        // end
    }

    public void Task_2_5(int[,] A, int[,] B)
    {
        // code here
        void FindMaxInColumn(int[,] matrix, int columnIndex, out int rowIndex) {
            rowIndex = 0;
            for (int i = 1; i < matrix.GetLength(0); i++)
                if (matrix[i, columnIndex] > matrix[rowIndex, columnIndex]) rowIndex = i;
        }
        // create and use FindMaxInColumn(matrix, columnIndex, out rowIndex);
        int columnIndex = 0;
        FindMaxInColumn(A, columnIndex, out int rowIndexA);
        FindMaxInColumn(B, columnIndex, out int rowIndexB);
        for (int j = 0; j < A.GetLength(1); j++)
            (A[rowIndexA, j], B[rowIndexB, j]) = (B[rowIndexB, j], A[rowIndexA, j]);
        // end
    }


    public void Task_2_6(ref int[] A, int[] B)
    {
        // code here

        // create and use FindMax(matrix, out row, out column); like in Task_2_1
        // create and use DeleteElement(array, index);

        // end
    }

    public void Task_2_7(ref int[,] B, int[,] C)
    {
        // code here

        // create and use CountRowPositive(matrix, rowIndex);
        int CountRowPositive(int[,] matrix, int rowIndex) {
            int count = 0;
            for (int j = 0; j < matrix.GetLength(1); j++)
                if (matrix[rowIndex, j] > 0) count++;
            return count;
        }
        int CountcolumnPositive(int[,] matrix, int colIndex) {
            int count = 0;
            for (int i = 0; i < matrix.GetLength(0); i++)
                if (matrix[i, colIndex] > 0) count++;
            return count;
       
        }
        // create and use CountColumnPositive(matrix, colIndex);
        int n1 = B.GetLength(0), m1 = B.GetLength(1), n2 = C.GetLength(0), m2 = C.GetLength(1);
        int maxB = 0, maxC = 0;
        int indexB = 0, indexC = 0;
        for (int i = 0; i < n1; i++) {
            int k = CountRowPositive(B, i);
            if (k > maxB) { maxB = k; indexB = i; }
        }
        for (int j = 0; j < m2; j++) {
            int k = CountcolumnPositive(C, j);
            if (k > maxC) { maxC = k; indexC = j; }
        }

        var newB = new int[n1+1, m1];

        for (int i = 0; i < n1; i++) 
            if (i <= indexB)
                for (int j = 0; j < m1; j++)
                    newB[i, j] = B[i, j];
            else
                for (int j = 0; j < m1; j++)
                    newB[i+1, j] = B[i, j];

        for (int j = 0; j < m1; j++)    
            newB[indexB+1, j] = C[j, indexC];
        B = newB;
        // end
    }

    public void Task_2_8(int[] A, int[] B)
    {
        // code here

        // create and use SortArrayPart(array, startIndex);

        // end
    }

    public int[] Task_2_9(int[,] A, int[,] C)
    {
        int[] answer = default(int[]);

        // code here
        int[] SumPositiveElementsInColumns(int[,] matrix) {
            var sumColumns = new int[matrix.GetLength(1)];
            for (int j = 0; j < matrix.GetLength(1); j++) 
                for (int i = 0; i < matrix.GetLength(0); i++)
                    if (matrix[i, j] > 0) sumColumns[j] += matrix[i, j];  
            return sumColumns;
            
        }
        // create and use SumPositiveElementsInColumns(matrix);
        int n1 = A.GetLength(0), m1 = A.GetLength(1), n2 = C.GetLength(0), m2 = C.GetLength(1);

        int[] sumA = new int[m1], sumC = new int[m2];
        answer = new int[m1+m2];
        sumA = SumPositiveElementsInColumns(A);
        sumC = SumPositiveElementsInColumns(C);
        for (int i = 0; i < m1; i++)
            answer[i] = sumA[i];

        for (int i = m1; i < m1+m2; i++)
            answer[i] = sumC[i-m1];
        // end

        return answer;
    }

    public void Task_2_10(ref int[,] matrix)
    {
        // code here

        // create and use RemoveColumn(matrix, columnIndex);

        // end
    }

    public void Task_2_11(int[,] A, int[,] B)
    {
        // code here

        // use FindMaxIndex(matrix, out row, out column); from Task_2_1
        int i1, j1, i2, j2;
        FindMaxIndex(A, out i1, out j1);
        FindMaxIndex(B, out i2, out j2);
        (A[i1, j1], B[i2, j2]) = (B[i2, j2], A[i1, j1]);
        // end
    }
    public void Task_2_12(int[,] A, int[,] B)
    {
        // code here

        // create and use FindMaxColumnIndex(matrix);

        // end
    }

    public void Task_2_13(ref int[,] matrix)
    {
        // code here

        // create and use RemoveRow(matrix, rowIndex);
        FindMaxIndex(matrix, out int rowMax, out int columnMax);
        FindMinIndex(matrix, out int rowMin, out int columnMin);
        switch (rowMax-rowMin){
            case 0: RemoveRow(ref matrix, rowMax); break;
            case > 0: RemoveRow(ref matrix, rowMax); RemoveRow(ref matrix, rowMin); break;
            case < 0: RemoveRow(ref matrix, rowMin); RemoveRow(ref matrix, rowMax); break;
        }

        // end
    }

    public void Task_2_14(int[,] matrix)
    {
        // code here

        // create and use SortRow(matrix, rowIndex);

        // end
    }

    public int Task_2_15(int[,] A, int[,] B, int[,] C)
    {
        int answer = 0;

        // code here
        double GetAverageWithoutMinMax(int[,] matrix) {
            double average = 0;
            int n = matrix.GetLength(0), m = matrix.GetLength(1);
            int max = matrix[0, 0], min = matrix[0, 0];
            for (int i = 0; i < n; i++)
                for (int j = 0; j < m; j++) {
                    if (matrix[i, j] > max) max = matrix[i, j];
                    if (matrix[i, j] < min) min = matrix[i, j];
                    average += matrix[i, j];
                }
            average -= max;
            average -= min;
            average /= n*m - 2;
            return average;
        }
        // create and use GetAverageWithoutMinMax(matrix);
        var array = new double[] {GetAverageWithoutMinMax(A), GetAverageWithoutMinMax(B), GetAverageWithoutMinMax(C)};
        if (array[0] <= array[1] && array[1] <= array[2]) answer = 1;
        else if (array[0] >= array[1] && array[1] >= array[2]) answer = -1;
        else answer = 0;
        // end

        // 1 - increasing   0 - no sequence   -1 - decreasing
        return answer;
    }

    public void Task_2_16(int[] A, int[] B)
    {
        // code here

        // create and use SortNegative(array);

        // end
    }

    public void Task_2_17(int[,] A, int[,] B)
    {
        // code here
        void SortRowsByMaxElement(int[,] matrix) {
            int n = matrix.GetLength(0), m = matrix.GetLength(1);
            int[] arrayOfMax = new int[n]; 
            for (int i = 0; i < n; i++) {
            int indexMax = 0;
            for (int j = 1; j < m; j++)
                if (matrix[i, j] > matrix[i, indexMax]) indexMax = j;
            arrayOfMax[i] = matrix[i, indexMax];
            }
            for (int pos = 1; pos < n; pos++) {
                int i = pos;
                while (i > 0 && arrayOfMax[i-1] < arrayOfMax[i]) {
                    int tmp = arrayOfMax[i-1];
                    arrayOfMax[i-1] = arrayOfMax[i];
                    arrayOfMax[i] = tmp;
                    for (int j = 0; j < m; j++) {
                        tmp = matrix[i-1, j];
                        matrix[i-1, j] = matrix[i, j];
                        matrix[i, j] = tmp;
                    }
                    i--;
                }
            }
                
        }
        // create and use SortRowsByMaxElement(matrix);
        SortRowsByMaxElement(A);
        SortRowsByMaxElement(B);
        // end
    }

    public void Task_2_18(int[,] A, int[,] B)
    {
        // code here

        // create and use SortDiagonal(matrix);

        // end
    }

    public void Task_2_19(ref int[,] matrix)
    {
        // code here
        int n = matrix.GetLength(0), m = matrix.GetLength(1);
        for (int i = n-1; i >= 0; i--) 
            for (int j = m-1; j >= 0; j--) {
                if (matrix[i, j] == 0) { RemoveRow(ref matrix, i); break; }
            }
        
        // use RemoveRow(matrix, rowIndex); from 2_13

        // end
    }
    public void Task_2_20(ref int[,] A, ref int[,] B)
    {
        // code here

        // use RemoveColumn(matrix, columnIndex); from 2_10

        // end
    }

    public void Task_2_21(int[,] A, int[,] B, out int[] answerA, out int[] answerB)
    {
        answerA = null;
        answerB = null;

        // code here
        int[] CreateArrayFromMins(int[,] matrix) {
            int n = matrix.GetLength(0);
            var arrayOfMin = new int[n];
            for (int i = 0; i < n; i++) {
                int min = matrix[i, i];
                for (int j = i+1; j < n; j++)
                    if (matrix[i, j] < min) min = matrix[i, j];
                arrayOfMin[i] = min;
            }
            return arrayOfMin;
        }
        // create and use CreateArrayFromMins(matrix);
        answerA = CreateArrayFromMins(A);
        answerB = CreateArrayFromMins(B);
        // end
    }

    public void Task_2_22(int[,] matrix, out int[] rows, out int[] cols)
    {
        rows = null;
        cols = null;

        // code here

        // create and use CountNegativeInRow(matrix, rowIndex);
        // create and use FindMaxNegativePerColumn(matrix);

        // end
    }

    public void Task_2_23(double[,] A, double[,] B)
    {
        // code here

        // create and use MatrixValuesChange(matrix);

        // end
    }

    public void Task_2_24(int[,] A, int[,] B)
    {
        // code here

        // use FindMaxIndex(matrix, out row, out column); like in 2_1
        // create and use SwapColumnDiagonal(matrix, columnIndex);

        // end
    }

    public void Task_2_25(int[,] A, int[,] B, out int maxA, out int maxB)
    {
        maxA = 0;
        maxB = 0;

        // code here

        // create and use FindRowWithMaxNegativeCount(matrix);
        // in FindRowWithMaxNegativeCount create and use CountNegativeInRow(matrix, rowIndex); like in 2_22

        // end
    }

    public void Task_2_26(int[,] A, int[,] B)
    {
        // code here

        // create and use FindRowWithMaxNegativeCount(matrix); like in 2_25
        // in FindRowWithMaxNegativeCount use CountNegativeInRow(matrix, rowIndex); from 2_22

        // end
    }

    public void Task_2_27(int[,] A, int[,] B)
    {
        // code here

        // create and use FindRowMaxIndex(matrix, rowIndex, out columnIndex);
        // create and use ReplaceMaxElementOdd(matrix, row, column);
        // create and use ReplaceMaxElementEven(matrix, row, column);

        // end
    }

    public void Task_2_28a(int[] first, int[] second, ref int answerFirst, ref int answerSecond)
    {
        // code here

        // create and use FindSequence(array, A, B); // 1 - increasing, 0 - no sequence,  -1 - decreasing
        // A and B - start and end indexes of elements from array for search

        // end
    }

    public void Task_2_28b(int[] first, int[] second, ref int[,] answerFirst, ref int[,] answerSecond)
    {
        // code here

        // use FindSequence(array, A, B); from Task_2_28a or entirely Task_2_28a
        // A and B - start and end indexes of elements from array for search

        // end
    }

    public void Task_2_28c(int[] first, int[] second, ref int[] answerFirst, ref int[] answerSecond)
    {
        // code here

        // use FindSequence(array, A, B); from Task_2_28a or entirely Task_2_28a or Task_2_28b
        // A and B - start and end indexes of elements from array for search

        // end
    }
    #endregion

    #region Level 3
    public void Task_3_1(ref double[,] firstSumAndY, ref double[,] secondSumAndY)
    {
        // code here

        // create and use public delegate SumFunction(x) and public delegate YFunction(x)
        // create and use method GetSumAndY(sFunction, yFunction, a, b, h);
        // create and use 2 methods for both functions calculating at specific x

        // end
    }

    public void Task_3_2(int[,] matrix)
    {
        // SortRowStyle sortStyle = default(SortRowStyle); - uncomment me

        // code here

        // create and use public delegate SortRowStyle(matrix, rowIndex);
        // create and use methods SortAscending(matrix, rowIndex) and SortDescending(matrix, rowIndex)
        // change method in variable sortStyle in the loop here and use it for row sorting

        // end
    }

    public double Task_3_3(double[] array)
    {
        double answer = 0;
        // SwapDirection swapper = default(SwapDirection); - uncomment me

        // code here

        // create and use public delegate SwapDirection(array);
        // create and use methods SwapRight(array) and SwapLeft(array)
        // create and use method GetSum(array, start, h) that sum elements with a negative indexes
        // change method in variable swapper in the if/else and than use swapper(matrix)

        // end

        return answer;
    }

    public int Task_3_4(int[,] matrix, bool isUpperTriangle)
    {
        int answer = 0;

        // code here

        // create and use public delegate GetTriangle(matrix);
        // create and use methods GetUpperTriange(array) and GetLowerTriange(array)
        // create and use GetSum(GetTriangle, matrix)

        // end

        return answer;
    }

    public void Task_3_5(out int func1, out int func2)
    {
        func1 = 0;
        func2 = 0;

        // code here

        // use public delegate YFunction(x, a, b, h) from Task_3_1
        // create and use method CountSignFlips(YFunction, a, b, h);
        // create and use 2 methods for both functions

        // end
    }

    public void Task_3_6(int[,] matrix)
    {
        // code here

        // create and use public delegate FindElementDelegate(matrix);
        // use method FindDiagonalMaxIndex(matrix) like in Task_2_3;
        // create and use method FindFirstRowMaxIndex(matrix);
        // create and use method SwapColumns(matrix, FindDiagonalMaxIndex, FindFirstRowMaxIndex);

        // end
    }

    public void Task_3_7(ref int[,] B, int[,] C)
    {
        // code here

        // create and use public delegate CountPositive(matrix, index);
        // use CountRowPositive(matrix, rowIndex) from Task_2_7
        // use CountColumnPositive(matrix, colIndex) from Task_2_7
        // create and use method InsertColumn(matrixB, CountRow, matrixC, CountColumn);

        // end
    }

    public void Task_3_10(ref int[,] matrix)
    {
        // FindIndex searchArea = default(FindIndex); - uncomment me

        // code here

        // create and use public delegate FindIndex(matrix);
        // create and use method FindMaxBelowDiagonalIndex(matrix);
        // create and use method FindMinAboveDiagonalIndex(matrix);
        // use RemoveColumn(matrix, columnIndex) from Task_2_10
        // create and use method RemoveColumns(matrix, findMaxBelowDiagonalIndex, findMinAboveDiagonalIndex)

        // end
    }

    public void Task_3_13(ref int[,] matrix)
    {
        // code here

        // use public delegate FindElementDelegate(matrix) from Task_3_6
        // create and use method RemoveRows(matrix, findMax, findMin)

        // end
    }

    public void Task_3_22(int[,] matrix, out int[] rows, out int[] cols)
    {

        rows = null;
        cols = null;

        // code here

        // create and use public delegate GetNegativeArray(matrix);
        // use GetNegativeCountPerRow(matrix) from Task_2_22
        // use GetMaxNegativePerColumn(matrix) from Task_2_22
        // create and use method FindNegatives(matrix, searcherRows, searcherCols, out rows, out cols);

        // end
    }

    public void Task_3_27(int[,] A, int[,] B)
    {
        // code here

        // create and use public delegate ReplaceMaxElement(matrix, rowIndex, max);
        // use ReplaceMaxElementOdd(matrix) from Task_2_27
        // use ReplaceMaxElementEven(matrix) from Task_2_27
        // create and use method EvenOddRowsTransform(matrix, replaceMaxElementOdd, replaceMaxElementEven);

        // end
    }

    public void Task_3_28a(int[] first, int[] second, ref int answerFirst, ref int answerSecond)
    {
        // code here

        // create public delegate IsSequence(array, left, right);
        // create and use method FindIncreasingSequence(array, A, B); similar to FindSequence(array, A, B) in Task_2_28a
        // create and use method FindDecreasingSequence(array, A, B); similar to FindSequence(array, A, B) in Task_2_28a
        // create and use method DefineSequence(array, findIncreasingSequence, findDecreasingSequence);

        // end
    }

    public void Task_3_28c(int[] first, int[] second, ref int[] answerFirstIncrease, ref int[] answerFirstDecrease, ref int[] answerSecondIncrease, ref int[] answerSecondDecrease)
    {
        // code here

        // create public delegate IsSequence(array, left, right);
        // use method FindIncreasingSequence(array, A, B); from Task_3_28a
        // use method FindDecreasingSequence(array, A, B); from Task_3_28a
        // create and use method FindLongestSequence(array, sequence);

        // end
    }
    #endregion
    #region bonus part
    public double[,] Task_4(double[,] matrix, int index)
    {
        // MatrixConverter[] mc = new MatrixConverter[4]; - uncomment me

        // code here

        // create public delegate MatrixConverter(matrix);
        // create and use method ToUpperTriangular(matrix);
        // create and use method ToLowerTriangular(matrix);
        // create and use method ToLeftDiagonal(matrix); - start from the left top angle
        // create and use method ToRightDiagonal(matrix); - start from the right bottom angle

        // end

        return matrix;
    }
    #endregion
}
