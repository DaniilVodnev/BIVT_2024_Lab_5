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
    long Factorial(int n) {
        int fc = 1;
        for (int i = 2; i <= n; i++)
            fc *= i;
        return fc;
    }
    long Combinations(int n, int k) {
        long fc = Factorial(n) / (Factorial(k) * Factorial(n-k));
        return fc;

    }
    public long Task_1_1(int n, int k)
    {
        long answer = 0;

        // code here
    
    
        // create and use Combinations(n, k);
        // create and use Factorial(n);
        if (n < 1 && Math.Abs(k) == 1) return 0;
        answer = Combinations(n, k);

        // end

        return answer;
    }
    double GeronArea(double a, double b, double c) {
            double p = (a + b + c) / 2;
        return Math.Sqrt(p*(p-a)*(p-b)*(p-c));
    }
    public int Task_1_2(double[] first, double[] second)
    {
        int answer = 0;

        // code here
        double a1 = first[0], b1 = first[1], c1 = first[2];
        double a2 = second[0], b2 = second[1], c2 = second[2];
        if (a1 >= b1 + c1 || b1 >= a1 + c1 || c1 >= a1 + b1) return -1;
        if (a2 >= b2 + c2 || b2 >= a2 + c2 || c2 >= a2 + b2) return -1;
        
        // create and use GeronArea(a, b, c);
        double area1 = GeronArea(a1, b1, c1), area2 = GeronArea(a2, b2, c2);
        if (area1 == area2) return 0;
        if (area1 > area2) return 1;
        if (area2 > area1) return 2;

        // end

        // first = 1, second = 2, equal = 0, error = -1
        return answer;
    }
    double GetDistance(double v, double a, double t)
        {
            return v * t + a * t * t / 2;
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
        answer++;
        while (GetDistance(v1, a1, answer) > GetDistance(v2, a2, answer)) 
            answer++;
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
    int FindDiagonalMaxIndex(int[,] matrix) {
        int index = 0;
        for (int i = 1; i < matrix.GetLength(0); i++)
            if (matrix[i, i] > matrix[index, index]) index = i;
        return index;
    }
    public void Task_2_3(ref int[,] B, ref int[,] C)
    {
        // code here
        
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
    void FindMaxInColumn(int[,] matrix, int columnIndex, out int rowIndex) {
        rowIndex = 0;
        for (int i = 1; i < matrix.GetLength(0); i++)
            if (matrix[i, columnIndex] > matrix[rowIndex, columnIndex]) rowIndex = i;
    }
    public void Task_2_5(int[,] A, int[,] B)
    {
        // code here
        
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
    public void Task_2_7(ref int[,] B, int[,] C)
    {
        // code here

        // create and use CountRowPositive(matrix, rowIndex);
        
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
    int[] SumPositiveElementsInColumns(int[,] matrix) {
        var sumColumns = new int[matrix.GetLength(1)];
        for (int j = 0; j < matrix.GetLength(1); j++) 
            for (int i = 0; i < matrix.GetLength(0); i++)
                if (matrix[i, j] > 0) sumColumns[j] += matrix[i, j];  
        return sumColumns;
        
    }
    public int[] Task_2_9(int[,] A, int[,] C)
    {
        int[] answer = default(int[]);

        // code here
        
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
    public int Task_2_15(int[,] A, int[,] B, int[,] C)
    {
        int answer = 0;

        // code here
        
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
    public void Task_2_17(int[,] A, int[,] B)
    {
        // code here
        
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
            for (int j = m-1; j >= 0; j--) 
                if (matrix[i, j] == 0) { RemoveRow(ref matrix, i); break; }
            
        
        // use RemoveRow(matrix, rowIndex); from 2_13

        // end
    }
    public void Task_2_20(ref int[,] A, ref int[,] B)
    {
        // code here

        // use RemoveColumn(matrix, columnIndex); from 2_10

        // end
    }
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
    public void Task_2_21(int[,] A, int[,] B, out int[] answerA, out int[] answerB)
    {
        answerA = null;
        answerB = null;
        // code here
        
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
    double[] MatrixToArray(double[,] matrix) {
        var array = new double[matrix.GetLength(0)*matrix.GetLength(1)];
        int i = 0;
        foreach (var el in matrix) {
            array[i] = el;
            i++;
        }
        return array;
    }
    void InsertSort(double[] array) {
        int n = array.Length;
        for (int pos = 1; pos < n; pos++) {
            int i = pos;
            while (i > 0 && array[i-1] < array[i]) {
                (array[i-1], array[i]) = (array[i], array[i-1]);
                i--;
            }
        }
    
    }
    double[] FiveMax(double[] array, int size) {
        int five = size >= 5 ? 5 : size;
        var fiveMax = new double[five];
        for (int i = 0; i < five; i++)
            fiveMax[i] = array[i];
        return fiveMax;   
    }
    void MatrixValuesChange(double[,] matrix, double[] fiveMax) {
        for (int i = 0; i < matrix.GetLength(0); i++) 
            for (int j = 0; j < matrix.GetLength(1); j++) 
                if (matrix[i, j] == fiveMax[0] || matrix[i, j] == fiveMax[1] || matrix[i, j] == fiveMax[2] || matrix[i, j] == fiveMax[3] || matrix[i, j] == fiveMax[4]) 
                    matrix[i, j] *= matrix[i, j] > 0 ? 2 : 0.5;
                else    
                    matrix[i, j] *= matrix[i, j] > 0 ? 0.5 : 2;
    }
    public void Task_2_23(double[,] A, double[,] B)
    {
        // code here
        double[] arrayA = MatrixToArray(A), arrayB = MatrixToArray(B);
        InsertSort(arrayA); InsertSort(arrayB);
        int sizeA = arrayA.Length, sizeB = arrayB.Length;
        double[] fiveMaxA = FiveMax(arrayA, sizeA), fiveMaxB = FiveMax(arrayB, sizeB);
        MatrixValuesChange(A, fiveMaxA); MatrixValuesChange(B, fiveMaxB);
        
        
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
    int[] CountNegativeInRow(int[,] matrix) {
        int n = matrix.GetLength(0), m = matrix.GetLength(1);
        var arrayOfCountNegative = new int[n];
        for (int i = 0; i < n; i++)
            for (int j = 0; j < m; j++)
                if (matrix[i, j] < 0) arrayOfCountNegative[i]++;
        return arrayOfCountNegative;
    }
    int FindRowWithMaxNegativeCount(int[] array) {
        int indexMax = 0;
        for (int i = 1; i < array.Length; i++)
            if (array[i] > array[indexMax]) indexMax = i;
        return indexMax;
    }
    public void Task_2_25(int[,] A, int[,] B, out int maxA, out int maxB)
    {
        maxA = 0;
        maxB = 0;

        // code here

        // create and use FindRowWithMaxNegativeCount(matrix);
        // in FindRowWithMaxNegativeCount create and use CountNegativeInRow(matrix, rowIndex); like in 2_22
        int[] arrayA = CountNegativeInRow(A), arrayB = CountNegativeInRow(B);
        maxA = FindRowWithMaxNegativeCount(arrayA); maxB = FindRowWithMaxNegativeCount(arrayB);

        // end
    }

    public void Task_2_26(int[,] A, int[,] B)
    {
        // code here

        // create and use FindRowWithMaxNegativeCount(matrix); like in 2_25
        // in FindRowWithMaxNegativeCount use CountNegativeInRow(matrix, rowIndex); from 2_22

        // end
    }
    void FindRowMaxIndex(int[,] matrix, int rowIndex, out int columnIndex) {
        columnIndex = 0;
        for (int j = 1; j < matrix.GetLength(1); j++) 
            if (matrix[rowIndex, j] > matrix[rowIndex, columnIndex]) columnIndex = j;
    }
    void ReplaceMaxElementOdd(int[,] matrix, int row, int column) {
        matrix[row, column] *= column + 1;
    }
    void ReplaceMaxElementEven(int[,] matrix, int row, int column) {
        matrix[row, column] = 0;
    }
    void MatrixChange(int[,] matrix) {
        int n = matrix.GetLength(0), m = matrix.GetLength(1);
        for (int i = 0; i < n; i++) {
            int columnIndex;
            FindRowMaxIndex(matrix, i, out columnIndex);
            if (i % 2 == 0) ReplaceMaxElementOdd(matrix, i, columnIndex);
            else ReplaceMaxElementEven(matrix, i, columnIndex);
        }
    }
    public void Task_2_27(int[,] A, int[,] B)
    {
        // code here

        // create and use FindRowMaxIndex(matrix, rowIndex, out columnIndex);
        // create and use ReplaceMaxElementOdd(matrix, row, column);
        // create and use ReplaceMaxElementEven(matrix, row, column);
        MatrixChange(A);
        MatrixChange(B);
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
    //public delegate double YFunction1(double x);
    //public delegate double SumFunction(int i, double x, ref int change);
    public delegate double YFunction(double x, double a, double b, double h);
    public delegate double SumFunction(int i, double x, ref int change);
    public delegate double YFunction1(double x);
    public double Factor(int i, double x, ref int f) {
        if (i > 0) f *= i;
        return Math.Cos(i*x)/f;
    }
    public double YFunc1(double x) {
        return Math.Exp(Math.Cos(x)) * Math.Cos(Math.Sin(x));
    }
    public double YFunc2(double x) {
        return (((x*x) - Math.PI * Math.PI/3))/4;
    }
    public double Siign(int i, double x, ref int sign) {
        sign *= -1;
        return sign * Math.Cos(i*x)/(i*i);
    }
    public double Calculate(SumFunction SFunc, double x, int i) {
        double s = 0;
        const double e = 0.0001;
        int g = 1;
        double member = SFunc(i, x, ref g);
        while (Math.Abs(member) > e) {
            s += member;
            member = SFunc(++i, x, ref g);
        }
        return s;
    }
    public void GetSumAndY(SumFunction SFunc, YFunction1 YFunc, double a, double b, double h, double[,] SumAndY, int first) {
        double x;
        double sum;
        double y;
        int lim = (int)((b - a) / h) + 1;
        for (int i = 0; i < lim; i++) {
            x = a + i * h;
            sum = Calculate(SFunc, x, first);
            y = YFunc(x);
            SumAndY[i, 0] = sum;
            SumAndY[i, 1] = y;
        }
    }
    public void Task_3_1(ref double[,] firstSumAndY, ref double[,] secondSumAndY){

        // code here 
        double a1 = 0.1, b1 = 1, h1 = 0.1;
        int size1 = (int)((b1 - a1) / h1) + 1;
        firstSumAndY = new double[size1, 2];
        GetSumAndY(Factor, YFunc1, a1, b1, h1, firstSumAndY, 0);
        //
        double a2 = Math.PI/5, b2 = Math.PI, h2 = Math.PI/25;
        int size2 = (int)((b2 - a2) / h2) + 1;
        secondSumAndY = new double[size2, 2];
        GetSumAndY(Siign, YFunc2, a2, b2, h2, secondSumAndY, 1);

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
    public delegate void SwapDirection(double[] array);

    void SwapRight(double[] array) {
        for (int i = 0; i < array.Length - 1; i += 2)
            (array[i], array[i + 1]) = (array[i + 1], array[i]);
    }

    void SwapLeft(double[] array) {
        for (int i = array.Length - 1; i > 0; i -= 2)
            (array[i], array[i - 1]) = (array[i - 1], array[i]);
    }

    double GetSum(double[] array) {
        double s = 0;
        int n = array.Length;
        for (int i = 1; i < n; i += 2)
            s += array[i];
        return s;
    }
    public double Task_3_3(double[] array)
    {
        double answer = 0;
        SwapDirection swapper = default(SwapDirection); 

        // code here
        double average = 0;
        foreach (double i in array) 
            average += i;
        average /= array.Length;
        SwapDirection swap = default(SwapDirection);
        swap = array[0] > average ? SwapRight : SwapLeft;
        swap(array);
        answer = GetSum(array);
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
    public double YFirst(double x, double a, double b, double h) {
        return x * x - Math.Sin(x);
    }
    public double YSecond(double x, double a, double b, double h) {
        return Math.Exp(x) - 1;
    }
    public int CountSignFlips(YFunction Yf, double a, double b, double h) {
        int c = 0;
        bool flag = Yf(a, a, b, h) < 0 ? true : false;
        int n = (int)((b - a) / h);
        for (int i = 1; i <= n; i++) {
            double tmp = Yf(a + i * h, a, b, h);
            if ((tmp < 0 && !flag) || (tmp > 0 && flag)) { flag = !flag; c++; }
        }
        return c;
    }
    public void Task_3_5(out int func1, out int func2)
    {
        func1 = 0;
        func2 = 0;

        // code here
        func1 = CountSignFlips(YFirst, 0, 2, 0.1);
        func2 = CountSignFlips(YSecond, -1, 1, 0.2);


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
    public delegate int CountPositive(int[,] matrix, int index);

    public void Task_3_7(ref int[,] B, int[,] C)
    {
        // code here
        CountPositive count = CountRowPositive;
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
    public delegate void FindElementDelegate(int[,] matrix, out int x, out int y);
    public int[,] RemoveRows(int[,] matrix, FindElementDelegate findMax, FindElementDelegate findMin) {
        FindMaxIndex(matrix, out int rowMax, out int columnMax);
        FindMinIndex(matrix, out int rowMin, out int columnMin);
        switch (rowMax-rowMin) {
            case 0: RemoveRow(ref matrix, rowMax); break;
            case > 0: RemoveRow(ref matrix, rowMax); RemoveRow(ref matrix, rowMin); break;
            case < 0: RemoveRow(ref matrix, rowMin); RemoveRow(ref matrix, rowMax); break;
        }
        return matrix;
    }
    public void Task_3_13(ref int[,] matrix)
    {
        // code here

        // use public delegate FindElementDelegate(matrix) from Task_3_6
        // create and use method RemoveRows(matrix, findMax, findMin)
        matrix = RemoveRows(matrix, FindMaxIndex, FindMinIndex);

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
    public delegate void ReplaceMaxElement(int[,] matrix, int row, int column);
    public void EvenOddRowsTransform(int[,] matrix, ReplaceMaxElement replaceMaxElementOdd, ReplaceMaxElement replaceMaxElementEven) {
        int n = matrix.GetLength(0), m = matrix.GetLength(1);
        for (int i = 0; i < n; i++) {
            int columnIndex;
            FindRowMaxIndex(matrix, i, out columnIndex);
            if (i % 2 == 0) ReplaceMaxElementOdd(matrix, i, columnIndex);
            else ReplaceMaxElementEven(matrix, i, columnIndex);
        }
    }
    public void Task_3_27(int[,] A, int[,] B)
    {
        // code here
        EvenOddRowsTransform(A, ReplaceMaxElementOdd, ReplaceMaxElementEven);
        EvenOddRowsTransform(B, ReplaceMaxElementOdd, ReplaceMaxElementEven);
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
