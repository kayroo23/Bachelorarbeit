import org.apache.commons.math3.linear.RealMatrix;

import java.util.ArrayList;
import java.util.List;

public class CalculationForAllPoints extends GeneralCalculation {

    protected double[][] calculateMatrix(List<List<Point>> newClusters, int numberOfPointsPerCluster){
        double[][] matrix = new double[numberOfPointsPerCluster * newClusters.size()][newClusters.get(0).get(0).getNumberOfAttributes()];

        int k = 0;
        for (List<Point> cluster : newClusters) {
            for (Point point: cluster) {
                matrix[k] = point.getAttributes();
                k++;
            }
        }
        return matrix;
    }

    protected List<Integer> calculateBestAttributes(int numberOfShownAttributes, RealMatrix spearmanMatrix){
        List<Integer> bestAttributes = new ArrayList<>();
        double[][] attributeValue = new double[numberOfShownAttributes][2];
        for(int i = 0; i < numberOfShownAttributes; i++){
            attributeValue[i][1] = 1;
        }
        for (int i = 0; i < spearmanMatrix.getColumnDimension(); i++) {
            double rowMax = 0;
            for(int j = 0; j < spearmanMatrix.getRowDimension(); j++){
                if(Math.abs(spearmanMatrix.getEntry(j , i)) > rowMax && Math.abs(spearmanMatrix.getEntry(j , i)) < 1){
                    rowMax = spearmanMatrix.getEntry(j , i);
                }
            }
            replaceWorst(Math.abs(rowMax), i, attributeValue);
        }
        for(int i = 0; i < attributeValue.length; i++){
            bestAttributes.add((int)attributeValue[i][0]);
        }
        return bestAttributes;
    }

    private static void replaceWorst(double rowMax, int rowMaxIndex, double[][] attributeValue){
        double max = 0;
        int arrayPosition = -1;
        for(int i = 0; i < attributeValue.length; i++){
            if(max < attributeValue[i][1]){
                max = attributeValue[i][1];
                arrayPosition = i;
            }
        }
        if(rowMax < max){
            attributeValue[arrayPosition][0] = rowMaxIndex;
            attributeValue[arrayPosition][1] = rowMax;
        }
    }

    protected void calculateResults(List<List<Point>> newClusters, List<Integer> bestAttributes){
        //Prints max and min value from the "best" attribute of each cluster
        double[][] result = new double[newClusters.size()][bestAttributes.size()*2];
        int k = 0;
        for (int x : bestAttributes) {
            double[][] tempMatrix = new double[newClusters.size()][2];
            for(int i = 0; i < newClusters.size(); i++){
                double min = Double.POSITIVE_INFINITY;
                double max = Double.NEGATIVE_INFINITY;
                for (Point p: newClusters.get(i)) {
                    if (p.getAttributes()[x] < min) {
                        min = p.getAttributes()[x];
                    }
                    if (p.getAttributes()[x] > max) {
                        max = p.getAttributes()[x];
                    }
                }
                tempMatrix[i][0] = min;
                tempMatrix[i][1] = max;
            }

            //Print the min/max Matrix
            System.out.println("Attribute: " + x);
            for(int i = 0; i < tempMatrix.length; i++){
                result[i][0 + 2 * k] = tempMatrix[i][0];
                result[i][1 + 2 * k] = tempMatrix[i][1];
            }
            k++;
        }
        for(int i = 0; i < result.length; i++){
            System.out.print("Cluster: "+ i + " : ");
            for(int j = 0; j < result[i].length; j++){
                System.out.print(result[i][j]);
                if(j % 2 == 1){
                    System.out.print(" ; ");
                }else{
                    System.out.print(" , ");
                }
            }
            System.out.println();
        }


    }

}
