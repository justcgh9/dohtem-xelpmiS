import java.util.Arrays;
import java.util.Scanner;

class Matrix {
    private int nCol;
    private int nRow;
    protected double[][] matrix;

    public Matrix(int nRow, int nCol) {
        this.nCol = nCol;
        this.nRow = nRow;
        this.matrix = new double[nRow][nCol];
    }

    public Matrix(double[][] matrix) {
        this.nCol = matrix[0].length;
        this.nRow = matrix.length;
        this.matrix = matrix;
    }

    public Matrix(double[] vector) {
        this.nCol = 1;
        this.nRow = vector.length;
        this.matrix = new double[nRow][nCol];
        for(int i = 0; i < nRow; i++) {
            this.matrix[i][0] = vector[i];
        }
    }

    public int getnCol() {
        return nCol;
    }

    public int getnRow() {
        return nRow;
    }

    public void multiplyByScalar(double scalar) {
        for(int i = 0; i < nRow; i++) {
            for(int j = 0; j < nCol; j++) {
                this.matrix[i][j] *= scalar;
            }
        }
    }

    public void add(Matrix other) {
        if(this.nCol != other.getnCol() || this.nRow != other.getnRow()) {
            throw new IllegalArgumentException("Matrices must have the same dimensions");
        }
        for(int i = 0; i < nRow; i++) {
            for(int j = 0; j < nCol; j++) {
                this.matrix[i][j] += other.matrix[i][j];
            }
        }
    }


    public Matrix multiply(Matrix other) {
        if (this.nCol != other.getnRow()) {
            throw new IllegalArgumentException("Matrices must have the same dimensions");
        }

        Matrix product = new Matrix(this.nRow, other.getnCol());

        for(int i = 0; i < this.nRow; i++) {
            for(int j = 0; j < other.getnCol(); j++) {
                for(int k = 0; k < this.nCol; k++) {
                    product.matrix[i][j] += this.matrix[i][k]*other.matrix[k][j];
                }
            }
        }
        return product;
    }

    public Matrix inverse() {
        if (this.nCol != this.nRow) {
            throw new IllegalArgumentException("Matrix must be square");
        }
        double determinant = this.determinant();
        if (determinant == 0) {
            throw new IllegalArgumentException("Matrix must be invertible");
        }


        Matrix inverse = new Matrix(this.nRow, this.nCol);
        for(int i = 0; i < this.nRow; i++) {
            for(int j = 0; j < this.nCol; j++) {
                inverse.matrix[i][j] = Math.pow(-1, i+j)*this.subMatrix(i, j).determinant()/determinant;
            }
        }

        return inverse.transpose();

    }

    public double determinant() {
        if (this.nCol != this.nRow) {
            throw new IllegalArgumentException("Matrix must be square");
        }
        double determinant = 0;
        if (this.nCol == 1) {
            determinant = this.matrix[0][0];
        } else if (this.nCol == 2) {
            determinant = this.matrix[0][0]*this.matrix[1][1] - this.matrix[0][1]*this.matrix[1][0];
        } else {
            for (int i = 0; i < this.nCol; i++) {
                determinant += Math.pow(-1, i)*this.matrix[0][i]*this.subMatrix(0, i).determinant();
            }
        }
        return determinant;
    }

    public Matrix subMatrix(int i, int j) {
        if (i > nRow || j > nCol) {
            throw new IllegalArgumentException("Index out of bounds");
        }
        Matrix subMatrix = new Matrix(this.nRow - 1, this.nCol - 1);
        int row = 0;
        int col = 0;
        for (int k = 0; k < this.nRow; k++) {
            col = 0;
            if (k != i) {
                for (int l = 0; l < this.nCol; l++) {
                    if (l != j) {
                        subMatrix.matrix[row][col] = this.matrix[k][l];
                        col++;
                    }
                }
                row++;
            }
        }
        return subMatrix;
    }

    public String toString() {
        String s = "";
        for (int i = 0; i < this.nRow; i++) {
            s += "[";
            for (int j = 0; j < this.nCol; j++) {
                s += this.matrix[i][j];
                if (j != this.nCol - 1) {
                    s += ", ";
                }
            }
            s += "]\n";
        }
        return s;
    }

    public Matrix subtract(Matrix cN) {
        if(this.nCol != cN.getnCol() || this.nRow != cN.getnRow()) {
            throw new IllegalArgumentException("Matrices must have the same dimensions");
        }
        Matrix difference = new Matrix(this.nRow, this.nCol);
        for(int i = 0; i < nRow; i++) {
            for(int j = 0; j < nCol; j++) {
                difference.matrix[i][j] = this.matrix[i][j] - cN.matrix[i][j];
            }
        }
        return difference;
    }

    public int rank() {
        int rank = 0;
        Matrix reducedRowEchelonForm = this.reducedRowEchelonForm();
        for(int i = 0; i < this.nRow; i++) {
            if(reducedRowEchelonForm.matrix[i][i] != 0) {
                rank++;
                break;
            }
        }
        return rank;
    }

    public Matrix transpose() {
        double[][] transpose = new double[this.nCol][this.nRow];
        for(int i = 0; i < this.nRow; i++) {
            for(int j = 0; j < this.nCol; j++) {
                transpose[j][i] = this.matrix[i][j];
            }
        }
        return new Matrix(transpose);
    }
    public Matrix reducedRowEchelonForm() {

        Matrix reducedRowEchelonForm = new Matrix(this.matrix);
        int pivot = 0;
        for(int i = 0; i < this.nRow; i++) {
            if(pivot >= this.nCol) {
                break;
            }
            if(reducedRowEchelonForm.matrix[i][pivot] == 0) {
                for(int j = i + 1; j < this.nRow; j++) {
                    if(reducedRowEchelonForm.matrix[j][pivot] != 0) {
                        reducedRowEchelonForm.swapRows(i, j);
                        break;
                    }
                }
            }
            if(reducedRowEchelonForm.matrix[i][pivot] != 0) {
                reducedRowEchelonForm.multiplyRow(i, 1/reducedRowEchelonForm.matrix[i][pivot]);
                for(int j = 0; j < this.nRow; j++) {
                    if(i != j && reducedRowEchelonForm.matrix[j][pivot] != 0) {
                        reducedRowEchelonForm.addRow(i, j, -reducedRowEchelonForm.matrix[j][pivot]);
                    }
                }
                pivot++;
            }
        }
        return reducedRowEchelonForm;
    }

    public void swapRows(int i, int j) {
        double[] temp = this.matrix[i];
        this.matrix[i] = this.matrix[j];
        this.matrix[j] = temp;
    }

    public void multiplyRow(int i, double scalar) {
        for(int j = 0; j < this.nCol; j++) {
            this.matrix[i][j] *= scalar;
        }
    }

    public static Matrix identityMatrix(int n) {
        Matrix identityMatrix = new Matrix(n, n);
        for(int i = 0; i < n; i++) {
            identityMatrix.matrix[i][i] = 1;
        }
        return identityMatrix;
    }

    public void addRow(int i, int j, double scalar) {
        for(int k = 0; k < this.nCol; k++) {
            this.matrix[j][k] += scalar*this.matrix[i][k];
        }
    }

    public double[] getColumn(int i) {
        double[] column = new double[this.nRow];
        for(int j = 0; j < this.nRow; j++) {
            column[j] = this.matrix[j][i];
        }
        return column;
    }

    public Matrix getRow(int i) {
        Matrix row = new Matrix(1, this.nCol);
        for(int j = 0; j < this.nCol; j++) {
            row.matrix[0][j] = this.matrix[i][j];
        }
        return row;
    }

    public void setRow(double[] row, int i) {
        this.matrix[i] = row;
    }

    public void setCol(double[] col, int i) {
        for(int j = 0; j < this.nRow; j++) {
            this.matrix[j][i] = col[j];
        }
    }
}

public class Main {


        public static void revisedSimplexAlgorithm(Matrix a, Matrix b, Matrix c,int numRows, int numCols, double accuracy, int [] pivots) {
            Matrix xB0 = new Matrix(numRows, numRows);
            int pivotIndexes[] = new int[numRows];
            for(int i = 0; i < numRows; i++) {
                for(int j = 0; j < numCols; j++) {
                    if(pivots[j] == i + 1) {
//                        xB0.setCol(a.getColumn(j), i);
                        pivotIndexes[i] = j;
                    }
                }
            }
            Matrix cB0 = new Matrix(numRows, 1);
            for(int i = 0; i < numRows; i++) {
                cB0.matrix[i][0] = c.matrix[pivotIndexes[i]][0];
            }


            Matrix B0 = Matrix.identityMatrix(numRows);

            Matrix B0Inverse;

            while (true) {
            B0Inverse = B0.inverse();

            xB0 = B0Inverse.multiply(b);

                Matrix deltas = cB0.transpose().multiply(B0Inverse).multiply(a).subtract(c.transpose());
                int enteringIndex = 0;

                for(int i = 0; i < deltas.getnCol(); i++) {
                    if(deltas.matrix[0][i] < deltas.matrix[0][enteringIndex]) {
                        enteringIndex = i;
                    }
                }
                if(deltas.matrix[0][enteringIndex] >= 0) {
                    double answer = cB0.transpose().multiply(B0Inverse).multiply(b).matrix[0][0];
                    System.out.println("Optimal solution found: " + answer);
                }

                int leavingIndex = 0;
                int trueLeavingIndex = 0;
                double minRatio = Double.MAX_VALUE;
                for(int i = 0; i < numRows; i++) {
                    if(xB0.matrix[i][0] > 0) {
                        double ratio1 = xB0.matrix[i][0];
                        double ratio2=a.matrix[i][enteringIndex];
                        double ratio = ratio1/ratio2;
                        System.out.println("Ratio1: " + ratio1 + " Ratio2: " + ratio2 + " Ratio: " + ratio);
                        if(ratio < minRatio) {
                            minRatio = ratio;
                            leavingIndex = i;
                        }
                    }
                }

                trueLeavingIndex = pivotIndexes[leavingIndex];
                pivotIndexes[leavingIndex] = enteringIndex;

                B0.setCol(a.getColumn(enteringIndex), leavingIndex);
                cB0.matrix[leavingIndex][0] = c.matrix[enteringIndex][0];

          }


        }

        public static int[] checkIdentity(Matrix a, int numRows, int numCols) {
            boolean identity = false;
            int[] pivots = new int[numCols];
            int numPivots = 0;
            for(int i = 0; i < numCols; i++) {
                int zeroCounter = 0;
                int oneCounter = 0;
                for(int j = 0 ; j < numRows; j++) {
                    if(a.matrix[j][i] ==  0) {
                        zeroCounter++;
                    } else if(a.matrix[j][i] == 1) {
                        oneCounter++;
                    }
                }
                if(oneCounter == 1 && zeroCounter == numRows - 1) {

                    pivots[i] = numPivots + 1;
                    numPivots++;
                } else {
                    pivots[i] = 0;
                }
            }
            if(numPivots == numRows) {
                return pivots;
            } else {
                return null;
            }
        }


        public static void main(String[] args) {
            Scanner scanner = new Scanner(System.in);

            System.out.print("Enter the number of constraints: ");
            int numRows = scanner.nextInt();
            System.out.print("Enter the number of variables: ");
            int numCols = scanner.nextInt();
            double[][] A = new double[numRows][numCols];
            System.out.println("Enter the coefficients of the constraint matrix A:");
            for (int i = 0; i < numRows; i++) {
                for (int j = 0; j < numCols; j++) {
                    A[i][j] = scanner.nextDouble();
                }
            }


            double[] b = new double[numRows];
            System.out.println("Enter the right-hand side vector b:");
            for (int i = 0; i < numRows; i++) {
                b[i] = scanner.nextDouble();
            }


            double[] C = new double[numCols];
            System.out.println("Enter the coefficients of the objective function vector C:");
            for (int j = 0; j < numCols; j++) {
                C[j] = scanner.nextDouble();
            }


            System.out.print("Enter the approximation accuracy: ");
            double accuracy = Double.parseDouble(scanner.next());
            int [] pivots = checkIdentity(new Matrix(A), numRows, numCols);
            if(pivots == null) {
                System.out.println("No valid basis in the given LPP");
                System.exit(0);
            }

            revisedSimplexAlgorithm(new Matrix(A), new Matrix(b), new Matrix(C), numRows, numCols, accuracy, pivots);



        }
}
