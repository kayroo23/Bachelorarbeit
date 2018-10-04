import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.stat.descriptive.rank.Percentile;

import javax.swing.*;
import java.awt.*;
import java.util.ArrayList;
import java.util.List;

 class GeneralCalculation {

    private List<double[][]> listOfMatrices = new ArrayList<>();

    double[][] calculateMatrix(List<List<Point>> newClusters){
        double[][] matrix = new double[newClusters.get(0).size() * newClusters.size()][newClusters.get(0).get(0).getNumberOfAttributes()];

        int k = 0;
        for (List<Point> cluster : newClusters) {
            for (Point point: cluster) {
                matrix[k] = point.getAttributes();
                k++;
            }
        }
        return matrix;
    }

    List<double[][]> calculateMatrixList(List<List<Point>> newClusters) {
        for (List<Point> cluster : newClusters) {
            int k = 0;
            double[][] matrix = new double[newClusters.get(0).size()][newClusters.get(0).get(0).getNumberOfAttributes()];
            for (Point point: cluster) {
                matrix[k] = point.getAttributes();
                k++;
            }
            listOfMatrices.add(matrix);
        }
        return listOfMatrices;
    }

    List<Integer> calculateBestAttributes(int numberOfShownAttributes, RealMatrix spearmanMatrix){
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

     Object[][] calculateMinMaxResults(List<List<Point>> newClusters, List<Integer> bestAttributes){
        //Prints max and min value from the "best" attribute of each cluster
        System.out.print("             ");
        Object[][] result = new Object[newClusters.size()][bestAttributes.size()*2];
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
            System.out.print("Attribute: " + x + "                             ");
            for(int i = 0; i < tempMatrix.length; i++){
                result[i][0 + 2 * k] = tempMatrix[i][0];
                result[i][1 + 2 * k] = tempMatrix[i][1];
            }
            k++;
        }
        System.out.println();
        for(int i = 0; i < result.length; i++){
            System.out.print("Cluster: "+ i + " : ");
            for(int j = 0; j < result[i].length; j++){
                System.out.print(result[i][j]);
                if(j % 2 == 1){
                    System.out.print(" | ");
                }else{
                    System.out.print(" , ");
                }
            }
            System.out.println();
        }
        System.out.println();

        return result;

    }

    Object[][] calculateQuartileResult(List<List<Point>> newClusters, List<Integer> bestAttributes){
        Percentile p = new Percentile(0.25);
        double[] attribute1 = new double[newClusters.get(0).size()];
        Object[][] result = new Object[newClusters.size()][bestAttributes.size()*2];
        for(int j = 0; j < newClusters.size(); j++){
            System.out.print("             ");
            for(int k = 0; k < bestAttributes.size(); k++){
                for(int i = 0; i < newClusters.get(j).size(); i++){
                    attribute1[i] = newClusters.get(j).get(i).getAttributes()[bestAttributes.get(k)];
                }
                p.setData(attribute1);
                result[j][0 + 2*k] = p.evaluate(25);
                result[j][1 + 2*k] = p.evaluate(75);

                System.out.print(p.evaluate(25)  + " , ");
                System.out.print(p.evaluate(75) + " | ");
            }
            System.out.println();
        }
        System.out.println();

        return result;
    }

    void printTable(Object[][] objectToPrint, List<Integer> bestAttributes, String title, String first, String second){
        String[] test = new String[objectToPrint[0].length];
        test[0] = "";
        for(int loopVariable = 0; loopVariable < test.length-1; loopVariable++){
            if(loopVariable % 2 == 0){
                test[loopVariable+1] = "Attribute " + bestAttributes.get(loopVariable/2) + ":  " + first;
            }else {
                test[loopVariable+1] = "Attribute " + bestAttributes.get(loopVariable/2) + ":  " + second;
            }
        }
        JTable table = new JTable(objectToPrint, test);
        JScrollPane scrollPane = new JScrollPane(table);
        JFrame frame = new JFrame(title);
        frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
        Container c = frame.getContentPane();
        c.add(scrollPane, BorderLayout.CENTER);
        frame.setVisible(true);
        frame.setSize(900,200);

    }

     /**Fuegt eine Anfangsbezeichnung hinzu und fügt die Tabellen zusammen
      *
      * @param obj1 Entspricht den gerade Zeilennummern
      * @param obj2 Entspricht den ungerade Zeilennummern
      * @return Ergebnis aus obj1 und obj2
      */
    Object[][] calculateTable(Object[][] obj1, Object[][] obj2){
        Object[][] objResult = new Object[obj1.length + obj2.length][obj1[0].length+1];
        for(int loop = 0; loop < objResult.length; loop++){
            if(loop % 2 == 0){
                objResult[loop][0] = "Cluster " + (loop/2) + " MinMax";
                for(int l = 0; l < obj1[0].length; l++){
                    objResult[loop][l+1] =  obj1[loop/2][l];
                }
            }else{
                objResult[loop][0] = "Cluster " + (loop/2) + " Quartile";
                for(int l = 0; l < obj2[0].length; l++){
                    objResult[loop][l+1] =  obj2[loop/2][l];
                }
            }

        }
        return objResult;
    }

}
