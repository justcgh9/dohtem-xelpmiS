//import java.util.Arrays;
//import java.util.Scanner;
//
//class Matrix {
//    private int nCol;
//    private int nRow;
//    private double[][] matrix;
//
//    public Matrix(int nRow, int nCol) {
//        this.nCol = nCol;
//        this.nRow = nRow;
//        this.matrix = new double[nRow][nCol];
//    }
//
//    public Matrix(double[][] matrix) {
//        this.nCol = matrix[0].length;
//        this.nRow = matrix.length;
//        this.matrix = matrix;
//    }
//
//    public Matrix(double[] vector) {
//        this.nCol = 1;
//        this.nRow = vector.length;
//        this.matrix = new double[nRow][nCol];
//        for(int i = 0; i < nRow; i++) {
//            this.matrix[i][0] = vector[i];
//        }
//    }
//
//    public int getnCol() {
//        return nCol;
//    }
//
//    public int getnRow() {
//        return nRow;
//    }
//
//    public void multiplyByScalar(double scalar) {
//        for(int i = 0; i < nRow; i++) {
//            for(int j = 0; j < nCol; j++) {
//                this.matrix[i][j] *= scalar;
//            }
//        }
//    }
//
//    public void add(Matrix other) {
//        if(this.nCol != other.getnCol() || this.nRow != other.getnRow()) {
//            throw new IllegalArgumentException("Matrices must have the same dimensions");
//        }
//        for(int i = 0; i < nRow; i++) {
//            for(int j = 0; j < nCol; j++) {
//                this.matrix[i][j] += other.matrix[i][j];
//            }
//        }
//    }
//
//
//    public Matrix multiply(Matrix other) {
//        if (this.nCol != other.getnRow()) {
//            throw new IllegalArgumentException("Matrices must have the same dimensions");
//        }
//
//        Matrix product = new Matrix(this.nRow, other.getnCol());
//
//        for(int i = 0; i < this.nRow; i++) {
//            for(int j = 0; j < other.getnCol(); j++) {
//                for(int k = 0; k < this.nCol; k++) {
//                    product.matrix[i][j] += this.matrix[i][k]*other.matrix[k][j];
//                }
//            }
//        }
//        return product;
//    }
//
//    public Matrix inverse() {
//        if (this.nCol != this.nRow) {
//            throw new IllegalArgumentException("Matrix must be square");
//        }
//        double determinant = this.determinant();
//        if (determinant == 0) {
//            throw new IllegalArgumentException("Matrix must be invertible");
//        }
//
//
//        Matrix inverse = new Matrix(this.nRow, this.nCol);
//        for(int i = 0; i < this.nRow; i++) {
//            for(int j = 0; j < this.nCol; j++) {
//                inverse.matrix[i][j] = Math.pow(-1, i+j)*this.subMatrix(i, j).determinant()/determinant;
//            }
//        }
//
//        return inverse;
//
//    }
//
//    public double determinant() {
//        if (this.nCol != this.nRow) {
//            throw new IllegalArgumentException("Matrix must be square");
//        }
//        double determinant = 0;
//        if (this.nCol == 1) {
//            determinant = this.matrix[0][0];
//        } else if (this.nCol == 2) {
//            determinant = this.matrix[0][0]*this.matrix[1][1] - this.matrix[0][1]*this.matrix[1][0];
//        } else {
//            for (int i = 0; i < this.nCol; i++) {
//                determinant += Math.pow(-1, i)*this.matrix[0][i]*this.subMatrix(0, i).determinant();
//            }
//        }
//        return determinant;
//    }
//
//    public Matrix subMatrix(int i, int j) {
//        if (i > nRow || j > nCol) {
//            throw new IllegalArgumentException("Index out of bounds");
//        }
//        Matrix subMatrix = new Matrix(this.nRow - 1, this.nCol - 1);
//        int row = 0;
//        int col = 0;
//        for (int k = 0; k < this.nRow; k++) {
//            col = 0;
//            if (k != i) {
//                for (int l = 0; l < this.nCol; l++) {
//                    if (l != j) {
//                        subMatrix.matrix[row][col] = this.matrix[k][l];
//                        col++;
//                    }
//                }
//                row++;
//            }
//        }
//        return subMatrix;
//    }
//
//    public String toString() {
//        String s = "";
//        for (int i = 0; i < this.nRow; i++) {
//            s += "[";
//            for (int j = 0; j < this.nCol; j++) {
//                s += this.matrix[i][j];
//                if (j != this.nCol - 1) {
//                    s += ", ";
//                }
//            }
//            s += "]\n";
//        }
//        return s;
//    }
//
//}
//
//
//
///**
// * A class for solving linear programming problems using the Revised Simplex Algorithm.
// */
//public class RevisedSimplexSolver {
//    private double[][] tableau; // The tableau representing the linear programming problem
//    private int[] basis; // Array to keep track of the basic variables
//
//    /**
//     * Constructor to initialize the solver with the linear programming problem.
//     *
//     * @param A         The matrix of coefficients of constraint functions.
//     * @param b         The vector of right-hand side numbers.
//     * @param C         The vector of coefficients of the objective function.
//     * @param accuracy  The approximation accuracy.
//     */
//    public RevisedSimplexSolver(double[][] A, double[] b, double[] C, double accuracy) {
//        // Initialize the tableau with the given data
//        int numRows = A.length;
//        int numCols = A[0].length + numRows; // Include slack variables
//        tableau = new double[numRows][numCols];
//        basis = new int[numRows];
//        double[] b_new = new double[numCols];
//        for(int i = 0; i < numCols; i++) {
//            b_new[i] = i < numRows ? b[i] : 0;
//        }
//
//        // Initialize the tableau with A and slack variables
//        for (int i = 0; i < numRows; i++) {
//            for (int j = 0; j < A[i].length; j++) {
//                tableau[i][j] = A[i][j];
//            }
//            tableau[i][A[i].length + i] = 1.0; // Slack variables
//            basis[i] = A[i].length + i;
//        }
//        int iteration = 0;
//        // Perform the Revised Simplex Algorithm
//        while (true) {
//            int enteringCol = selectEnteringVariable(C);
//            if (enteringCol == -1) {
//                if(iteration == 0) {
//                    System.out.println("Method is not applicable");
//
//                }
//                return;
//            }
//            int leavingRow = selectLeavingVariable(enteringCol, b_new);
//            if (leavingRow == -1) {
//                if(iteration == 0)
//                System.out.println("Method is not applicable");
//                return;
//            }
//            iteration++;
//            pivot(enteringCol, leavingRow);
//        }
//    }
//
//    /**
//     * Check if the current solution is optimal.
//     *
//     * @param C         The vector of coefficients of the objective function.
//     * @param accuracy  The approximation accuracy.
//     * @return True if the solution is optimal; otherwise, false.
//     */
//    private boolean isOptimal(double[] C, double accuracy) {
//        for (double coefficient : C) {
//            if (coefficient < -accuracy) {
//                return false;
//            }
//        }
//        return true;
//    }
//
//    /**
//     * Select the entering variable with the most negative coefficient in the objective function.
//     *
//     * @param C The vector of coefficients of the objective function.
//     * @return The index of the entering variable.
//     */
//    private int selectEnteringVariable(double[] C) {
//        // Select the entering variable with the most negative coefficient in the objective function
//        int enteringCol = -1;
//        double max_neg = 0;
//        for(int i = 0; i < tableau[0].length; i++) {
//            double intermediateProduct = 1;
//            for(int j = 0; j < tableau.length; j++) {
//                intermediateProduct += tableau[j][i]*C[j];
//            }
//            enteringCol = intermediateProduct < max_neg ? i : enteringCol;
//            max_neg = intermediateProduct < max_neg ? intermediateProduct : max_neg;
//        }
//        return enteringCol;
//    }
//
//    /**
//     * Select the leaving variable using the ratio test.
//     *
//     * @param enteringCol The index of the entering variable.
//     * @return The index of the leaving variable.
//     */
//    private int selectLeavingVariable(int enteringCol, double[] b) {
//        // Select the leaving variable using the ratio test
//        int leavingRow = -1;
//        double minRatio = Double.POSITIVE_INFINITY;
//        for (int i = 0; i < tableau.length; i++) {
//            if (tableau[i][enteringCol] > 0) {
//                double ratio =  b[i]/ tableau[i][enteringCol];
//                if (ratio < minRatio) {
//                    minRatio = ratio;
//                    leavingRow = i;
//                }
//            }
//        }
//        return leavingRow;
//    }
//
//    /**
//     * Perform the pivot operation to update the tableau.
//     *
//     * @param enteringCol The index of the entering variable.
//     * @param leavingRow  The index of the leaving variable.
//     */
//    private void pivot(int enteringCol, int leavingRow) {
//        // Perform pivot operation to update the tableau
//        double pivotElement = tableau[leavingRow][enteringCol];
//        for (int j = 0; j < tableau[0].length; j++) {
//            tableau[leavingRow][j] /= pivotElement;
//        }
//        for (int i = 0; i < tableau.length; i++) {
//            if (i != leavingRow) {
//                double multiplier = tableau[i][enteringCol];
//                for (int j = 0; j < tableau[0].length; j++) {
//                    tableau[i][j] -= multiplier * tableau[leavingRow][j];
//                }
//            }
//        }
//        basis[leavingRow] = enteringCol;
//    }
//
//    /**
//     * Get the solution vector.
//     *
//     * @return The solution vector.
//     */
//    public double[] getSolution() {
//        // Extract the solution vector from the tableau
//        int numCols = tableau[0].length;
//        double[] solution = new double[numCols - tableau.length];
//        Arrays.fill(solution, 0.0);
//
//        for (int i = 0; i < tableau.length; i++) {
//            solution[basis[i]] = tableau[i][numCols - 1];
//        }
//        return solution;
//    }
//
//    /**
//     * Calculate the objective function value.
//     *
//     * @param C The vector of coefficients of the objective function.
//     * @return The objective function value.
//     */
//    public double getObjectiveValue(double[] C) {
//        // Calculate the objective function value
//        double[] solution = getSolution();
//        double objectiveValue = 0.0;
//        for (int j = 0; j < C.length; j++) {
//            objectiveValue += C[j] * solution[j];
//        }
//        return objectiveValue;
//    }
//
//
//    public static void main(String[] args) {
//        Scanner scanner = new Scanner(System.in);
//
//
//        System.out.print("Enter the number of constraints: ");
//        int numRows = scanner.nextInt();
//        System.out.print("Enter the number of variables: ");
//        int numCols = scanner.nextInt();
//
//
//        double[][] A = new double[numRows][numCols];
//        System.out.println("Enter the coefficients of the constraint matrix A:");
//        for (int i = 0; i < numRows; i++) {
//            for (int j = 0; j < numCols; j++) {
//                A[i][j] = scanner.nextDouble();
//            }
//        }
//
//
//        double[] b = new double[numRows];
//        System.out.println("Enter the right-hand side vector b:");
//        for (int i = 0; i < numRows; i++) {
//            b[i] = scanner.nextDouble();
//        }
//
//
//        double[] C = new double[numCols];
//        System.out.println("Enter the coefficients of the objective function vector C:");
//        for (int j = 0; j < numCols; j++) {
//            C[j] = scanner.nextDouble();
//        }
//
//
//        System.out.print("Enter the approximation accuracy: ");
//        double accuracy = Double.parseDouble(scanner.next());
//
//        System.out.println("c");
//
//        scanner.close();
//
//        System.out.println("a");
//
//        RevisedSimplexSolver solver = new RevisedSimplexSolver(A, b, C, accuracy);
//        double[] solution = solver.getSolution();
//        System.out.println("*");
//        double objectiveValue = solver.getObjectiveValue(C);
//
//        System.out.println("b");
//
//        System.out.println("Optimal Solution:");
//        System.out.println(Arrays.toString(solution));
//        System.out.println("Objective Function Value: " + objectiveValue);
//    }
//
//}
