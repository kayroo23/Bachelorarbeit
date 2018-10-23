import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.stat.correlation.SpearmansCorrelation;
import org.apache.commons.math3.stat.descriptive.moment.GeometricMean;
import org.apache.commons.math3.stat.descriptive.moment.StandardDeviation;
import org.apache.commons.math3.stat.descriptive.moment.Variance;
import org.apache.commons.math3.stat.descriptive.rank.Percentile;

import javax.swing.*;
import java.awt.*;
import java.util.ArrayList;
import java.util.List;

 class GeneralCalculation {

    private List<double[][]> listOfMatrices = new ArrayList<>();

     RealMatrix calculateStandardDeviation(double[][] matrix){
         double[] input = new double[matrix[0].length];
         RealMatrix result = new Array2DRowRealMatrix(input);
         double[] temp = new double[matrix[0].length];
         StandardDeviation stdabw = new StandardDeviation();

         for(int i = 0; i < matrix[0].length; i++){
             for(int k = 0; k < temp.length; k++){
                 temp[k] = matrix[k][i];
             }
             result.setEntry(i,0, stdabw.evaluate(temp));
         }
        return result;
     }

     RealMatrix calculateGeomMean(double[][] matrix){
         double[] input = new double[matrix[0].length];
         RealMatrix result = new Array2DRowRealMatrix(input);
         double[] temp = new double[matrix[0].length];
         GeometricMean geomMean = new GeometricMean();

         for(int i = 0; i < matrix[0].length; i++){
             for(int k = 0; k < temp.length; k++){
                 temp[k] = matrix[k][i];
             }
             result.setEntry(i,0, geomMean.evaluate(temp));
         }
         return result;
     }

     RealMatrix calculateVariance(double[][] matrix){
         double[] input = new double[matrix[0].length];
         RealMatrix result = new Array2DRowRealMatrix(input);
         double[] temp = new double[matrix[0].length];
         Variance variance = new Variance();

         for(int i = 0; i < matrix[0].length; i++){
             for(int k = 0; k < temp.length; k++){
                 temp[k] = matrix[k][i];
             }
             result.setEntry(i,0, variance.evaluate(temp));
         }
         return result;
     }

     List<Integer> calculateBestAttributesForVectors(int numberOfShownAttributes, RealMatrix vector){
         List<Integer> bestAttributes = new ArrayList<>();
         RealMatrix attributeValue = new Array2DRowRealMatrix(new double[vector.getRowDimension()]);
         attributeValue = vector;
         for(int l = 0; l < numberOfShownAttributes; l++){
             double min = Double.MAX_VALUE;
             int position = -1;
             for(int k = 0; k < attributeValue.getRowDimension(); k++){
                 if(attributeValue.getEntry(k,0) < min){
                     min = attributeValue.getEntry(k,0);
                     position = k;
                 }
             }
             bestAttributes.add(position);
             attributeValue.setEntry(position,0, Double.MAX_VALUE);
         }
         return bestAttributes;
     }


     /**Gibt die Anzahl an Punkten im kleinsten Cluster zurueck
      *
      * @param newClusters Liste an Clustern
      * @return Anzahl an Punkten im kleinsten Cluster zurueck
      */
    private int getMinimumClusterSize(List<List<Point>> newClusters){
        int min = Integer.MAX_VALUE;
        for(int i = 0; i < newClusters.size(); i++){
            if(newClusters.get(i).size() < min){
                min = newClusters.get(i).size();
            }
        }
        return min;
    }

    double[][] calculateMatrix(List<List<Point>> newClusters){
        int min = getMinimumClusterSize(newClusters);
        double[][] matrix = new double[min * newClusters.size()][newClusters.get(0).get(0).getNumberOfAttributes()];
        for (List<Point> cluster : newClusters) {
            for(int l = 0; l < min; l++){
                matrix[l] = cluster.get(l).getAttributes();
            }
        }
        return matrix;
    }

    private List<double[][]> calculateMatrixList(List<List<Point>> newClusters) {
        int min = getMinimumClusterSize(newClusters);
        for (List<Point> cluster : newClusters) {
            double[][] matrix = new double[min][newClusters.get(0).get(0).getNumberOfAttributes()];
            for (int l = 0; l < min; l++) {
                matrix[l] = cluster.get(l).getAttributes();
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
        //System.out.print("             ");
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
            //System.out.print("Attribute: " + x + "                             ");
            for(int i = 0; i < tempMatrix.length; i++){
                result[i][0 + 2 * k] = tempMatrix[i][0];
                result[i][1 + 2 * k] = tempMatrix[i][1];
            }
            k++;
        }
        /*System.out.println();
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
        System.out.println();*/

        return result;

    }

    Object[][] calculateQuartileResult(List<List<Point>> newClusters, List<Integer> bestAttributes){
        Percentile p = new Percentile();
        int min = getMinimumClusterSize(newClusters);
        double[] attribute1 = new double[min];
        Object[][] result = new Object[newClusters.size()][bestAttributes.size()*2];
        for(int j = 0; j < newClusters.size(); j++){
            //System.out.print("             ");
            for(int k = 0; k < bestAttributes.size(); k++){
                for(int i = 0; i < min; i++){
                    attribute1[i] = newClusters.get(j).get(i).getAttributes()[bestAttributes.get(k)];
                }

                p.setData(attribute1);
                result[j][0 + 2*k] = p.evaluate(25);
                result[j][1 + 2*k] = p.evaluate(75);

                //System.out.print(p.evaluate(25)  + " , ");
                //System.out.print(p.evaluate(75) + " | ");
            }
            //System.out.println();
        }
        //System.out.println();

        return result;
    }

     /**Calculates the overlap of all pairs
      *
      * @param newClusters List of clusters
      * @param bestAttributes List of best Attributes
      */
    void calculateOverlap(List<List<Point>> newClusters, List<Integer> bestAttributes){
        Object[][] calc = this.calculateMinMaxResults(newClusters,bestAttributes);
        double meanOfMeans = 0;
        for (int best = 0; best < bestAttributes.size(); best++) {
            //System.out.println(bestAttributes.get(best) + "   :");

            float mean = 0;
            int count = 0;
            for(int l = 0; l < calc.length-1; l++){
                for(int k = l+1; k < calc.length; k++){
                    mean += this.overlap((double)calc[l][(2*best)], (double)calc[l][(2*best)+1],
                            (double)calc[k][(2*best)], (double)calc[k][(2*best)+1]);
                    count++;
                    //System.out.println(this.overlap((double)calc[l][(2*best)], (double)calc[l][(2*best)+1],
                            //(double)calc[k][(2*best)], (double)calc[k][(2*best)+1]));
                }
            }
            System.out.println();
            mean = mean/count;
            /*System.out.println("Mean: " + mean);
            System.out.println();
            System.out.println();*/
            meanOfMeans += mean;
        }
        meanOfMeans = meanOfMeans/bestAttributes.size();
        System.out.println("Mean of Means: " + meanOfMeans);
        System.out.println();
        System.out.println();
    }

    private double overlap(double xMin, double xMax, double yMin, double yMax){
        return (Math.max(0, Math.max(xMax-yMin, yMax-xMin))/Math.min(xMax-xMin,yMax-yMin));
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

     /**errechnet die Tabellen der einzelnen Cluster und fügt sie zusammen
      *
      * @param newClusters Liste an Clustern
      * @param bestAttributes Liste an besten Attributen
      * @param numberOfShownAttributes Anzahl Attribute die angezeigt werden sollen
      * @return
      */
    Object[][] calculateTablePerClusterWithSpearman(List<List<Point>> newClusters, List<Integer> bestAttributes, int numberOfShownAttributes){
        List<double[][]> matrixList;
        matrixList = this.calculateMatrixList(newClusters);
        Object[][] objResult2 = new Object[newClusters.size()*4 - 2][];
        Object[][] test = new Object[newClusters.size()][(bestAttributes.size()*2)+1];
        Object[] empty = new Object[(bestAttributes.size()*2)+1];

        for(int l = 0; l < empty.length; l++){
            empty[l] = "";
        }

        for(int i = 0; i < newClusters.size(); i++) {
            //pearsonMatrix = new PearsonsCorrelation().computeCorrelationMatrix(matrixList.get(i));
            RealMatrix spearmanMatrix = new SpearmansCorrelation().computeCorrelationMatrix(matrixList.get(i));
            bestAttributes = this.calculateBestAttributes(numberOfShownAttributes, spearmanMatrix);
            List<List<Point>> oneCluster = new ArrayList<>();
            oneCluster.add(newClusters.get(i));
            Object[][] objResult1 = this.calculateTable(this.calculateMinMaxResults(oneCluster, bestAttributes), this.calculateQuartileResult(oneCluster, bestAttributes));

            test[i][0] = "";
            for (int loopVariable = 0; loopVariable < test[0].length - 1; loopVariable++) {
                if (loopVariable % 2 == 0) {
                    test[i][loopVariable + 1] = "Attribute " + bestAttributes.get(loopVariable / 2) + ":  " + "First";
                } else {
                    test[i][loopVariable + 1] = "Attribute " + bestAttributes.get(loopVariable / 2) + ":  " + "Second";
                }
                objResult2[4 * i] = objResult1[0];
                objResult2[(4 * i) + 1] = objResult1[1];
                if (i + 1 < newClusters.size()) {
                    objResult2[(4 * i) + 2] = empty;
                    objResult2[(4 * i) + 3] = test[i + 1];
                }
            }
        }

        return objResult2;
    }

     /** Errechnet die besten Attribute für das erste Cluster
      *
      * @param newClusters Liste an Clustern
      * @param numberOfShownAttributes Anzahl Attribute die angezeigt werden sollen
      * @return bestAttribute List
      */
    List<Integer> bestAttributesForFirstClusterWithSpearman(List<List<Point>> newClusters, int numberOfShownAttributes){
        List<double[][]> matrixList;
        matrixList = this.calculateMatrixList(newClusters);
        RealMatrix spearmanMatrix = new SpearmansCorrelation().computeCorrelationMatrix(matrixList.get(0));
        return this.calculateBestAttributes(numberOfShownAttributes, spearmanMatrix);
    }

     /** Errechnet die besten Attribute wenn sie zuerst für alle Cluster einzeln berechnet werden und davon die häufigsten ausgewählt werden
      *
      * @param newClusters Liste an Clustern
      * @param numberOfShownAttributes Anzahl Attribute die angezeigt werden sollen
      * @return bestAttribute List
      */
    List<Integer> bestAttributesFirstForClusterThenForAll(List<List<Point>> newClusters, int numberOfShownAttributes){
        List<double[][]> matrixList;
        matrixList = this.calculateMatrixList(newClusters);
        int[] findBestAttributes = new int[newClusters.get(0).get(0).getNumberOfAttributes()];
        for(int j = 0; j < findBestAttributes.length; j++){
            findBestAttributes[j] = 0;
        }
        List<Integer> bestAttributes = new ArrayList<>();
        for(int i = 0; i < newClusters.size(); i++){
            RealMatrix spearmanMatrix = new SpearmansCorrelation().computeCorrelationMatrix(matrixList.get(i));
            bestAttributes = this.calculateBestAttributes(numberOfShownAttributes, spearmanMatrix);
            for (int x : bestAttributes) {
                findBestAttributes[x]++;
            }
        }
        bestAttributes.clear();
        for(int k = 0; k < numberOfShownAttributes; k++){
            int highest = 0;
            for (int l = 0; l < findBestAttributes.length; l++) {
                if((findBestAttributes[l] > highest) && !bestAttributes.contains(l)){
                    highest = l;
                }
            }
            bestAttributes.add(highest);
        }
        return bestAttributes;
    }

}
