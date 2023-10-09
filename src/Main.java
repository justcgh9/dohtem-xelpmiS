import java.util.Arrays;
import java.util.Scanner;

class Matrix {
    private int nCol;
    private int nRow;
    private double[][] matrix;

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

        return inverse;

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

}

public class Main {

        public static void main(String[] args) {
            Scanner scanner = new Scanner(System.in);

            System.out.print("Enter the number of constraints: ");
            int nRows = scanner.nextInt();
            System.out.print("Enter the number of variables: ");
            int nCols = scanner.nextInt();


        }
}
