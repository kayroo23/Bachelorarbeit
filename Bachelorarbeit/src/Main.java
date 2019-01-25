import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.ml.clustering.CentroidCluster;
import org.apache.commons.math3.ml.clustering.KMeansPlusPlusClusterer;
import org.apache.commons.math3.stat.correlation.Covariance;
import org.apache.commons.math3.stat.correlation.KendallsCorrelation;
import org.apache.commons.math3.stat.correlation.PearsonsCorrelation;
import org.apache.commons.math3.stat.correlation.SpearmansCorrelation;
import org.jfree.ui.RefineryUtilities;

import java.io.*;
import java.util.ArrayList;
import java.util.List;
import java.util.Random;
import java.util.Scanner;

import java.io.File;
import java.io.FileNotFoundException;

public class Main {

    /**
     *
     * @param positionPairs eine Liste an Punkten, mit denen die uebergebene x und y Koordinate verglichen wird
     * @param positionX x-Koordinate of the point
     * @param positionY y-Koordinate of the point
     * @param border grenzt ein, wie nah zwei verschiedene Cluster sein duerfen
     * @return true, falls true zurueck, falls der Abstand des PositionPairs zu dem Punkt mit den Koordinaten positionX
     * und positionY
     * kleiner als border ist.
     * Ansonsten true
     */
    private static boolean tooClose(ArrayList<NewPair> positionPairs, int positionX, int positionY, int border){
        for (NewPair p : positionPairs) {
            double pythagoras = Math.sqrt((Math.pow(positionX - p.getFirst(), 2) + (Math.pow(positionY - p.getSecond(), 2))));
            if(pythagoras < border){
                return true;
            }
        }
        return false;
    }

    private static boolean isNotValidCenter(double[] center, ArrayList<Point> listOfCenters, double border){
        for (Point p: listOfCenters) {
            double sum = 0;
            for (int i = 0; i < center.length; i++) {
                sum += Math.pow(center[i] - p.getAttributes()[i],2);
            }
            if(Math.sqrt(sum) < border){
                return true;
            }
        }
        return false;
    }

    /**
     *
     * @param streuungsreichweite gibt die Verschiebung an
     * @param bereich gibt den moeglichen Ergebnisbereich an
     * @return Zufallszahl
     */
    private static int getRandomPosition(int streuungsreichweite, int bereich){
        return (int)((streuungsreichweite/2) + Math.random() * bereich);
    }

    /**
     * Fuellt die uebergebene Liste an Clustern mit Fakedaten
     * @param clusters Liste an Clustern
     * @param numberOfPointsPerCluster Anzahl der Datenpunkte innerhalb eines Clusters
     */
    private static void generateFakeData(List<List<NewPair>> clusters, int numberOfPointsPerCluster){
        int streuungsreichweite = 15;
        int bereich = (int)(clusters.size() * streuungsreichweite);
        int xK = getRandomPosition(streuungsreichweite, bereich);
        int yK = getRandomPosition(streuungsreichweite, bereich);
        ArrayList<NewPair> clusterPositions = new ArrayList<>();
        int first;
        int second;
        for (List<NewPair> x : clusters) {

            while(tooClose(clusterPositions, xK, yK, 2 * streuungsreichweite)){
                xK = getRandomPosition(streuungsreichweite, bereich);
                yK = getRandomPosition(streuungsreichweite, bereich);
            }

            for(int i = 0; i < numberOfPointsPerCluster; i++){
                first = (int)((xK - (streuungsreichweite/2)) + (Math.random() * streuungsreichweite));
                second = (int)((yK - (streuungsreichweite/2)) + (Math.random() * streuungsreichweite));
                NewPair p = new NewPair(first, second);
                x.add(p);
            }
            clusterPositions.add(new NewPair(xK,yK));
        }
    }

    private static void getRandomCenter(double[] center, int bereich){
        for(int i = 0; i < center.length; i++){
            center[i] = Math.random() * bereich;
        }
    }

    /**
     * Fuellt die uebergebene Liste an Clustern mit Fakedaten
     * @param clusters Liste an Clustern
     * @param numberOfPointsPerCluster Anzahl der Datenpunkte innerhalb eines Clusters
     * @param numberOfAtributes Anzahl an Attributen die jeder Punkt besitzt
     */
    private static void generateFakeData(List<List<Point>> clusters, int numberOfPointsPerCluster, int numberOfAtributes){
        double streuungsreichweite = 15;
        int bereich = (int)(clusters.size() * streuungsreichweite * 100);
        ArrayList<Point> clusterPositions = new ArrayList<>();

        for (List<Point> cluster: clusters) {
            double[] center = new double[numberOfAtributes];
            do{
                getRandomCenter(center, bereich);
            }while(isNotValidCenter(center, clusterPositions, 3 * streuungsreichweite));
            for(int i = 0; i < numberOfPointsPerCluster; i++){
                double[] pointArray = new double[numberOfAtributes];
                for(int k = 0; k < numberOfAtributes; k++){
                    pointArray[k] = (center[k] - (streuungsreichweite/2))  + (Math.random() * streuungsreichweite);
                }
                Point p = new Point(numberOfAtributes, pointArray);
                cluster.add(p);
                clusterPositions.add(new Point(numberOfAtributes, center));
            }
        }
    }

    private static void calculationFor2Attributes(List<List<NewPair>> clusters, int numberOfPointsPerCluster){
        new DrawPoints(clusters);

        double[] firstCoord = new double[numberOfPointsPerCluster * clusters.size()];
        double[] secondCoord = new double[numberOfPointsPerCluster * clusters.size()];
        double[][] matrix = new double[numberOfPointsPerCluster * clusters.size()][2];

        int k = 0;
        for (List<NewPair> x : clusters) {
            for (NewPair pair: x) {
                firstCoord[k] = pair.getFirst();
                secondCoord[k] = pair.getSecond();
                k++;
            }
        }

        //Fuellt die matrix mit den Punkten der Cluster
        for(int i = 0; i < numberOfPointsPerCluster * clusters.size(); i++){
            matrix[i][0] = firstCoord[i];
            matrix[i][1] = secondCoord[i];
        }

        double covariance = new Covariance().covariance(firstCoord,secondCoord);
        double pearson = new PearsonsCorrelation().correlation(firstCoord, secondCoord);
        double spearman = new SpearmansCorrelation().correlation(firstCoord, secondCoord);
        System.out.println("Covariance: " + covariance);
        System.out.println("Pearsons correlation coefficient: " + pearson);
        System.out.println("Spearmans rank correlation coefficient: " + spearman);

        //double covarianceMatrix = new Covariance().computeCovarianceMatrix(all);
        RealMatrix pearsonMatrix = new PearsonsCorrelation().computeCorrelationMatrix(matrix);
        RealMatrix spearmanMatrix = new SpearmansCorrelation().computeCorrelationMatrix(matrix);

        System.out.println("Pearson Matrix: ");
        for (int i = 0; i < pearsonMatrix.getColumnDimension(); i++) {
            for(int j = 0; j < pearsonMatrix.getRowDimension(); j++){
                System.out.print(pearsonMatrix.getEntry(j , i) + "   ");
            }
            System.out.println(" ");
        }
        System.out.println("Spearman Matrix: ");
        for (int i = 0; i < spearmanMatrix.getColumnDimension(); i++) {
            for(int j = 0; j < spearmanMatrix.getRowDimension(); j++){
                System.out.print(spearmanMatrix.getEntry(j , i) + "   ");
            }
            System.out.println(" ");
        }
    }

    private static void calculationForMoreAttributes(List<List<Point>> newClusters, int numberOfPointsPerCluster){
        double[][] matrix = new double[numberOfPointsPerCluster * newClusters.size()][newClusters.get(0).get(0).getNumberOfAttributes()];

        int k = 0;
        for (List<Point> cluster : newClusters) {
            for (Point point: cluster) {
                matrix[k] = point.getAttributes();
                k++;
            }
        }

        //Prints all Points sorted by the clusters
        for(int i = 0; i < matrix.length; i++){
            if(i % numberOfPointsPerCluster == 0){
                System.out.println("Next Cluster");
            }
            for(int j = 0; j < matrix[0].length; j++){
                System.out.print(matrix[i][j] + "   ");
            }
            System.out.println(" ");

        }
        System.out.println(" ");


        RealMatrix pearsonMatrix = new PearsonsCorrelation().computeCorrelationMatrix(matrix);
        RealMatrix spearmanMatrix = new SpearmansCorrelation().computeCorrelationMatrix(matrix);


        //Prints the Pearson Matrix
        System.out.println("Pearson Matrix: ");
        for (int i = 0; i < pearsonMatrix.getColumnDimension(); i++) {
            for(int j = 0; j < pearsonMatrix.getRowDimension(); j++){
                System.out.print(pearsonMatrix.getEntry(j , i) + "   ");
                if(i == j){
                    System.out.print("                ");
                }
            }
            System.out.println();
        }
        System.out.println();

        //Prints the Spearman Matrix
        System.out.println("Spearman Matrix: ");
        for (int i = 0; i < spearmanMatrix.getColumnDimension(); i++) {
            for(int j = 0; j < spearmanMatrix.getRowDimension(); j++){
                System.out.print(spearmanMatrix.getEntry(j , i) + "   ");
                if(i == j){
                    System.out.print("                ");
                }
            }
            System.out.println();
        }
        System.out.println();


        //Calculates the most significant attribute of pearson
        double pearsonColumnMin = 1;
        int pearsonAttribute = -1;
        for (int i = 0; i < pearsonMatrix.getColumnDimension(); i++) {
            double rowMax = 0;
            for(int j = 0; j < pearsonMatrix.getRowDimension(); j++){
                if(Math.abs(pearsonMatrix.getEntry(j , i)) > rowMax && Math.abs(pearsonMatrix.getEntry(j , i)) < 1){
                    rowMax = pearsonMatrix.getEntry(j , i);
                }
            }
            if(rowMax < pearsonColumnMin){
                pearsonColumnMin = rowMax;
                pearsonAttribute = i;
            }
        }
        System.out.println();
        System.out.println("Attribut: " + pearsonAttribute + " Pearson: " + pearsonColumnMin);

        //Calculates the most significant attribute of spearman
        double spearmanColumnMin = 1;
        int spearmanAttribute = -1;
        for (int i = 0; i < spearmanMatrix.getColumnDimension(); i++) {
            double rowMax = 0;
            for(int j = 0; j < spearmanMatrix.getRowDimension(); j++){
                if(Math.abs(spearmanMatrix.getEntry(j , i)) > rowMax && Math.abs(spearmanMatrix.getEntry(j , i)) < 1){
                    rowMax = spearmanMatrix.getEntry(j , i);
                }
            }
            if(rowMax < spearmanColumnMin){
                spearmanColumnMin = rowMax;
                spearmanAttribute = i;
            }
        }
        System.out.println();
        System.out.println("Attribut: " + spearmanAttribute + " Spearman: " + spearmanColumnMin);
        System.out.println();


        //Prints max and min value from the "best" attribute of each cluster
        double[][] tempMatrix = new double[newClusters.size()][2];
        for(int i = 0; i < newClusters.size(); i++){
            double min = Double.POSITIVE_INFINITY;
            double max = Double.NEGATIVE_INFINITY;
            for (Point p: newClusters.get(i)) {
                if (p.getAttributes()[spearmanAttribute] < min) {
                    min = p.getAttributes()[spearmanAttribute];
                }
                if (p.getAttributes()[spearmanAttribute] > max) {
                    max = p.getAttributes()[spearmanAttribute];
                }
            }
            tempMatrix[i][0] = min;
            tempMatrix[i][1] = max;
        }


        //Sort the min/max Matrix with respect to min
        for(int i = tempMatrix.length; i > 1; --i){
            double temp;
            for(int j = 0; j < i-1; ++j){
                if(tempMatrix[j][0] > tempMatrix[j+1][0]){
                    temp = tempMatrix[j][0];
                    tempMatrix[j][0] = tempMatrix[j+1][0];
                    tempMatrix[j+1][0] = temp;

                    temp = tempMatrix[j][1];
                    tempMatrix[j][1] = tempMatrix[j+1][1];
                    tempMatrix[j+1][1] = temp;
                }
            }
        }

        //Print the min/max Matrix
        for(int i = 0; i < tempMatrix.length; i++){
            System.out.println("Cluster: " + i + ": [" + tempMatrix[i][0] + "," + tempMatrix[i][1] + "]");
        }



    }

    private static double[][] calculateMatrix(List<List<Point>> newClusters, int numberOfPointsPerCluster){
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

    private static List<Integer> calculateBestAttributes(int numberOfShownAttributes, RealMatrix spearmanMatrix){
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

    private static void calculateResults(List<List<Point>> newClusters, List<Integer> bestAttributes){
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


    /**reads the txt file and stores the points in seperated clusters
     *
     * @param location location of the txt file
     * @param newClusters List of clusters
     * @param firstRelevantValue all Values from the input file before this value will be ignored
     */
    private static void readClassificationDataSet(String location, List<List<Point>> newClusters, int firstRelevantValue){
        File file = new File(location);
        newClusters.clear();
        List<String> knownClasses = new ArrayList<>();
        try{
            String s;
            Scanner sc = new Scanner(file);
            int positionInArrayList = 0;
            while(sc.hasNextLine()){
                while(sc.hasNext()){
                    s = sc.next();
                    String[] values = s.split(",");
                    //Add new cluster in Arraylist if necessary
                    if(!knownClasses.contains(values[values.length-1])){
                        knownClasses.add(values[values.length-1]);
                        newClusters.add(new ArrayList<Point>());
                    }
                    //Find right cluster for this point
                    for (int z = 0; z < knownClasses.size(); z++) {
                        if(knownClasses.get(z).equals(values[values.length-1])){
                            positionInArrayList = z;
                            break;
                        }
                    }
                    double[] attributes = new double[values.length - 1 - firstRelevantValue];
                    for(int l = firstRelevantValue; l < values.length-1; l++){
                        attributes[l - firstRelevantValue] = Double.parseDouble(values[l]);
                    }
                    newClusters.get(positionInArrayList).add(new Point(attributes.length, attributes));
                }
            }
        }catch(FileNotFoundException e){
            System.out.println("Einlesen der Datei fehlgeschlagen!");
        }
    }

    /**
     * reads CSV in location and converts every line in one Point without the last attribute (clusterlable)
     * @param location to the CSV-File
     * @param clusters to store the points
     */
    private static void readCSV(String location, List<Point> clusters){
        File file = new File(location);
        try{
            int i = 0;
            String s;
            Scanner sc = new Scanner(file);
            while(sc.hasNextLine()){
                //Exit when line is empty
                if(!sc.hasNext()){
                    break;
                }
                while(sc.hasNext()){
                    s = sc.next();
                    String[] values = s.split(",");
                    double[] attributes = new double[values.length - 1];
                    for(int l = 0; l < values.length - 1; l++){
                        attributes[l] = Double.parseDouble(values[l]);
                    }
                    clusters.add(new Point(attributes.length, attributes));
                }
            }
        }catch(FileNotFoundException e){
            System.out.println("Einlesen der Datei fehlgeschlagen!");
        }
    }

    public static void main(String[] args) {
        List<List<NewPair>> clusters = new ArrayList<>();
        List<List<Point>> newClusters = new ArrayList<>();
        int numberOfClusters = 5;
        int numberOfPointsPerCluster = 1000;
        int numberOfAttributes = 10;
        for(int i = 0; i < numberOfClusters; i++){
            clusters.add(new ArrayList<NewPair>());
            newClusters.add(new ArrayList<Point>());
        }
        //generateFakeData(clusters, numberOfPointsPerCluster);
        //calculationFor2Attributes(clusters, numberOfPointsPerCluster);

        //generateFakeData(newClusters, numberOfPointsPerCluster, numberOfAttributes);
        //newClusters = generateRandomGoldStandardDataset();
        //newClusters = generateGaussGoldStandardDataset(0.5, numberOfPointsPerCluster);
        newClusters = new generateSimpleData('n', 0.5, numberOfClusters, numberOfPointsPerCluster).getClusters();
        //calculationForMoreAttributes(newClusters, numberOfPointsPerCluster);

        //readClassificationDataSet("dataset_32_pendigits_changed.txt", newClusters, 0);
        //readClassificationDataSet("irisDataset.txt", newClusters, 0);
        //readClassificationDataSet("fertilityDataSet.txt", newClusters, 0);
        //readClassificationDataSet("ecoliDataSet.txt", newClusters, 1);

        int numberOfShownAttributes = 4;
        List<Integer> bestAttributes = new ArrayList<>();

        GeneralCalculation calculator = new GeneralCalculation();
        double[][] matrix = calculator.calculateMatrix(newClusters);


        RealMatrix pearsonMatrix = new PearsonsCorrelation().computeCorrelationMatrix(matrix);
        RealMatrix spearmanMatrix = new SpearmansCorrelation().computeCorrelationMatrix(matrix);
        RealMatrix kendallsTauMatrix = new KendallsCorrelation().computeCorrelationMatrix(matrix);

        RealMatrix standardDeviation = calculator.calculateStandardDeviation(matrix);
        RealMatrix geomMean = calculator.calculateGeomMean(matrix);
        RealMatrix variance = calculator.calculateVariance(matrix);
        RealMatrix medianDeviation = calculator.calculateMedianDeviation(matrix);
        RealMatrix variationsCoefficient = calculator.calculateVariationsCoefficient(matrix);
        RealMatrix quartilsDispersion = calculator.calculateQuartilsDispersion(matrix);

        //printMatrix(standardDeviation, "stdabw");
        //printMatrix(variance, "variance");
        //printMatrix(pearsonMatrix, "Pearson Matrix: ");
        //printMatrix(spearmanMatrix, "Spearman Matrix: ");
        //printMatrix(kendallsTauMatrix, "Kendalls Tau Matrix: ");


        //Calculates the best attributes in general with matrices
        /*
        bestAttributes = calculator.calculateBestAttributesForMatrix(numberOfShownAttributes, spearmanMatrix);
        Object[][] objResult = calculator.calculateTable(calculator.calculateMinMaxResults(newClusters, bestAttributes),
                calculator.calculateQuartileResult(newClusters, bestAttributes));
        calculator.printTable(objResult, bestAttributes, "Spearman for all Points: MinMax + Quartile",
                "First", "Second");
        System.out.println("For all attributes with Spearman");
        System.out.println(calculator.calculateOverlapInterval(newClusters, bestAttributes));

        bestAttributes = calculator.calculateBestAttributesForMatrix(numberOfShownAttributes, pearsonMatrix);
        Object[][] objResult1 = calculator.calculateTable(calculator.calculateMinMaxResults(newClusters, bestAttributes),
                calculator.calculateQuartileResult(newClusters, bestAttributes));
        calculator.printTable(objResult1, bestAttributes, "Pearson for all Points: MinMax + Quartile",
                "First", "Second");
        System.out.println("For all attributes with Pearson");
        System.out.println(calculator.calculateOverlapInterval(newClusters, bestAttributes));

        bestAttributes = calculator.calculateBestAttributesForMatrix(numberOfShownAttributes, kendallsTauMatrix);
        Object[][] objResult2 = calculator.calculateTable(calculator.calculateMinMaxResults(newClusters, bestAttributes),
                calculator.calculateQuartileResult(newClusters, bestAttributes));
        calculator.printTable(objResult2, bestAttributes, "Kendall for all Points: MinMax + Quartile",
                "First", "Second");
        System.out.println("For all attributes with Kendall");
        System.out.println(calculator.calculateOverlapInterval(newClusters, bestAttributes));



        //Calculates the best attributes in general with vector
        bestAttributes = calculator.calculateMinAttributesForVectors(numberOfShownAttributes, standardDeviation);
        Object[][] objResultVector1 = calculator.calculateTable(calculator.calculateMinMaxResults(newClusters, bestAttributes),
                calculator.calculateQuartileResult(newClusters, bestAttributes));
        calculator.printTable(objResultVector1, bestAttributes, "Stdabw for all Points: MinMax + Quartile",
                "First", "Second");
        System.out.println("For all attributes with Stdabw");
        System.out.println(calculator.calculateOverlapInterval(newClusters, bestAttributes));

        bestAttributes = calculator.calculateMinAttributesForVectors(numberOfShownAttributes, geomMean);
        Object[][] objResultVector2 = calculator.calculateTable(calculator.calculateMinMaxResults(newClusters, bestAttributes),
                calculator.calculateQuartileResult(newClusters, bestAttributes));
        calculator.printTable(objResultVector2, bestAttributes, "Geometric Mean for all Points: MinMax + Quartile",
                "First", "Second");
        System.out.println("For all attributes with Geometric Mean");
        System.out.println(calculator.calculateOverlapInterval(newClusters, bestAttributes));

        bestAttributes = calculator.calculateMinAttributesForVectors(numberOfShownAttributes, medianDeviation);
        Object[][] objResultVector4 = calculator.calculateTable(calculator.calculateMinMaxResults(newClusters, bestAttributes),
                calculator.calculateQuartileResult(newClusters, bestAttributes));
        calculator.printTable(objResultVector4, bestAttributes, "MedianDeviation for all Points: MinMax + Quartile",
                "First", "Second");
        System.out.println("For all attributes with MedianDeviation");
        System.out.println(calculator.calculateOverlapInterval(newClusters, bestAttributes));

        bestAttributes = calculator.calculateMinAttributesForVectors(numberOfShownAttributes, variance);
        Object[][] objResultVector3 = calculator.calculateTable(calculator.calculateMinMaxResults(newClusters, bestAttributes),
                calculator.calculateQuartileResult(newClusters, bestAttributes));
        calculator.printTable(objResultVector3, bestAttributes, "Variance for all Points: MinMax + Quartile",
                "First", "Second");
        System.out.println("For all attributes with Variance");
        System.out.println(calculator.calculateOverlapInterval(newClusters, bestAttributes));

        bestAttributes = calculator.calculateMinAttributesForVectors(numberOfShownAttributes, variationsCoefficient);
        Object[][] objResultVector5 = calculator.calculateTable(calculator.calculateMinMaxResults(newClusters, bestAttributes),
                calculator.calculateQuartileResult(newClusters, bestAttributes));
        calculator.printTable(objResultVector5, bestAttributes, "variationsCoefficient for all Points: MinMax + Quartile",
                "First", "Second");
        System.out.println("For all attributes with variationsCoefficient");
        System.out.println(calculator.calculateOverlapInterval(newClusters, bestAttributes));

        bestAttributes = calculator.calculateMinAttributesForVectors(numberOfShownAttributes, quartilsDispersion);
        Object[][] objResultVector6 = calculator.calculateTable(calculator.calculateMinMaxResults(newClusters, bestAttributes),
                calculator.calculateQuartileResult(newClusters, bestAttributes));
        calculator.printTable(objResultVector6, bestAttributes, "quartilsDispersion for all Points: MinMax + Quartile",
                "First", "Second");
        System.out.println("For all attributes with quartilsDispersion");
        System.out.println(calculator.calculateOverlapInterval(newClusters, bestAttributes));
        */





        //Calculates the best attributes per cluster
        /*
        List<double[][]> matrixList1 = calculator.calculateMatrixList(newClusters);
        List<RealMatrix> spearmanList = new ArrayList<>();
        for(int l = 0; l < matrixList1.size(); l++){
            spearmanList.add(new SpearmansCorrelation().computeCorrelationMatrix(matrixList1.get(l)));
        }
        for(int a = 0; a < spearmanList.size(); a++){
            for(int u = 0; u < spearmanList.get(a).getRowDimension(); u++){
                for(int v = 0; v < spearmanList.get(a).getColumnDimension(); v++){
                    if(Double.isNaN(spearmanList.get(a).getEntry(u,v))){
                        spearmanList.get(a).setEntry(u,v,0);
                    }
                }
            }
        }

        Object[][] objResult3 = calculator.calculateTablePerClusterWithMatrix(newClusters, numberOfShownAttributes, spearmanList);
        bestAttributes = calculator.calculateBestAttributesForMatrix(numberOfShownAttributes, spearmanList.get(0));
        calculator.printTable(objResult3, bestAttributes, "Spearman per Cluster: MinMax + Quartile ", "First",
                "Second");



        List<double[][]> matrixList2 = calculator.calculateMatrixList(newClusters);
        List<RealMatrix> varianceList = new ArrayList<>();
        for(int l = 0; l < matrixList2.size(); l++){
            varianceList.add(calculator.calculateVariance(matrixList2.get(l)));
        }
        Object[][] objResult31 = calculator.calculateMinimizingTablePerClusterWithVector(newClusters, numberOfShownAttributes, varianceList);
        bestAttributes = calculator.calculateMinAttributesForVectors(numberOfShownAttributes, varianceList.get(0));
        calculator.printTable(objResult31, bestAttributes, "Variance per Cluster: MinMax + Quartile ", "First",
                "Second");


        //Variance of all points - variance of cluster
        List<double[][]> matrixList21 = calculator.calculateMatrixList(newClusters);
        List<RealMatrix> varianceList2 = new ArrayList<>();
        for(int l = 0; l < matrixList21.size(); l++){
            RealMatrix tempVariance = calculator.calculateVariance(matrixList21.get(l));
            for(int k = 0; k < tempVariance.getRowDimension(); k++){
                tempVariance.setEntry(k,0, -(Math.abs(variance.getEntry(k,0)) - Math.abs(tempVariance.getEntry(k,0))));
            }
            varianceList2.add(tempVariance);

        }
        Object[][] objResult311 = calculator.calculateMinimizingTablePerClusterWithVector(newClusters, numberOfShownAttributes, varianceList2);
        bestAttributes = calculator.calculateMinAttributesForVectors(numberOfShownAttributes, varianceList2.get(0));
        calculator.printTable(objResult311, bestAttributes, "Variance difference per Cluster: MinMax + Quartile ", "First",
                "Second");


        //Stdabw of all points - stdabw of cluster
        List<double[][]> matrixList41 = calculator.calculateMatrixList(newClusters);
        List<RealMatrix> stdabwDif = new ArrayList<>();
        for(int l = 0; l < matrixList41.size(); l++){
            RealMatrix tempStdabw = calculator.calculateStandardDeviation(matrixList41.get(l));
            for(int k = 0; k < tempStdabw.getRowDimension(); k++){
                tempStdabw.setEntry(k,0, -(Math.abs(standardDeviation.getEntry(k,0)) - Math.abs(tempStdabw.getEntry(k,0))));
            }
            stdabwDif.add(tempStdabw);

        }
        Object[][] objResult41 = calculator.calculateMinimizingTablePerClusterWithVector(newClusters, numberOfShownAttributes, stdabwDif);
        bestAttributes = calculator.calculateMinAttributesForVectors(numberOfShownAttributes, stdabwDif.get(0));
        calculator.printTable(objResult41, bestAttributes, "Stdabw difference per Cluster: MinMax + Quartile ", "First",
                "Second");



        List<double[][]> matrixList3 = calculator.calculateMatrixList(newClusters);
        List<RealMatrix> stdabwList = new ArrayList<>();
        for(int l = 0; l < matrixList3.size(); l++){
            stdabwList.add(calculator.calculateStandardDeviation(matrixList3.get(l)));
        }
        Object[][] objResult32 = calculator.calculateMinimizingTablePerClusterWithVector(newClusters, numberOfShownAttributes, stdabwList);
        bestAttributes = calculator.calculateMinAttributesForVectors(numberOfShownAttributes, stdabwList.get(0));
        calculator.printTable(objResult32, bestAttributes, "Stdabw per Cluster: MinMax + Quartile ", "First",
                "Second");



        List<double[][]> matrixList4 = calculator.calculateMatrixList(newClusters);
        List<RealMatrix> variationskoeffizientList = new ArrayList<>();
        for(int l = 0; l < matrixList4.size(); l++){
            variationskoeffizientList.add(calculator.calculateVariationsCoefficient(matrixList4.get(l)));
        }
        Object[][] objResult33 = calculator.calculateMinimizingTablePerClusterWithVector(newClusters, numberOfShownAttributes, variationskoeffizientList);
        bestAttributes = calculator.calculateMinAttributesForVectors(numberOfShownAttributes, variationskoeffizientList.get(0));
        calculator.printTable(objResult33, bestAttributes, "Variationskoeffizient per Cluster: MinMax + Quartile ", "First",
                "Second");



        //VarCofDif of all points - VarCofDif of cluster
        List<double[][]> matrixList42 = calculator.calculateMatrixList(newClusters);
        List<RealMatrix> varCofDif = new ArrayList<>();
        for(int l = 0; l < matrixList41.size(); l++){
            RealMatrix tempVarCofDif = calculator.calculateVariationsCoefficient(matrixList42.get(l));
            for(int k = 0; k < tempVarCofDif.getRowDimension(); k++){
                tempVarCofDif.setEntry(k,0, -(Math.abs(variationsCoefficient.getEntry(k,0)) - Math.abs(tempVarCofDif.getEntry(k,0))));
            }
            varCofDif.add(tempVarCofDif);

        }
        Object[][] objResult42 = calculator.calculateMinimizingTablePerClusterWithVector(newClusters, numberOfShownAttributes, varCofDif);
        bestAttributes = calculator.calculateMinAttributesForVectors(numberOfShownAttributes, varCofDif.get(0));
        calculator.printTable(objResult42, bestAttributes, "VarCof Difference per Cluster: MinMax + Quartile ", "First",
                "Second");


        List<double[][]> matrixList5 = calculator.calculateMatrixList(newClusters);
        List<RealMatrix> geometricMeanList = new ArrayList<>();
        for(int l = 0; l < matrixList5.size(); l++){
            geometricMeanList.add(calculator.calculateGeomMean(matrixList5.get(l)));
        }
        Object[][] objResult34 = calculator.calculateMinimizingTablePerClusterWithVector(newClusters, numberOfShownAttributes, geometricMeanList);
        bestAttributes = calculator.calculateMinAttributesForVectors(numberOfShownAttributes, geometricMeanList.get(0));
        calculator.printTable(objResult34, bestAttributes, "GeometricMean per Cluster: MinMax + Quartile ", "First",
                "Second");



        List<double[][]> matrixList6 = calculator.calculateMatrixList(newClusters);
        List<RealMatrix> medianDeviationList = new ArrayList<>();
        for(int l = 0; l < matrixList6.size(); l++){
            medianDeviationList.add(calculator.calculateMedianDeviation(matrixList6.get(l)));
        }
        Object[][] objResult35 = calculator.calculateMinimizingTablePerClusterWithVector(newClusters, numberOfShownAttributes, medianDeviationList);
        bestAttributes = calculator.calculateMinAttributesForVectors(numberOfShownAttributes, medianDeviationList.get(0));
        calculator.printTable(objResult35, bestAttributes, "MedianDeviation per Cluster: MinMax + Quartile ", "First",
                "Second");


        //MedianDeviation of all points - MedianDeviation of cluster
        List<double[][]> matrixList61 = calculator.calculateMatrixList(newClusters);
        List<RealMatrix> medDefDif = new ArrayList<>();
        for(int l = 0; l < matrixList61.size(); l++){
            RealMatrix tempMedDefDif = calculator.calculateMedianDeviation(matrixList61.get(l));
            for(int k = 0; k < tempMedDefDif.getRowDimension(); k++){
                tempMedDefDif.setEntry(k,0, -(Math.abs(medianDeviation.getEntry(k,0)) - Math.abs(tempMedDefDif.getEntry(k,0))));
            }
            medDefDif.add(tempMedDefDif);

        }
        Object[][] objResult61 = calculator.calculateMinimizingTablePerClusterWithVector(newClusters, numberOfShownAttributes, medDefDif);
        bestAttributes = calculator.calculateMinAttributesForVectors(numberOfShownAttributes, medDefDif.get(0));
        calculator.printTable(objResult61, bestAttributes, "MedianDeviation difference per Cluster: MinMax + Quartile ", "First",
                "Second");






        List<double[][]> matrixList7 = calculator.calculateMatrixList(newClusters);
        List<RealMatrix> quartilsDispersionsList = new ArrayList<>();
        for(int l = 0; l < matrixList7.size(); l++){
            quartilsDispersionsList.add(calculator.calculateQuartilsDispersion(matrixList7.get(l)));
        }
        Object[][] objResult36 = calculator.calculateMinimizingTablePerClusterWithVector(newClusters, numberOfShownAttributes, quartilsDispersionsList);
        bestAttributes = calculator.calculateMinAttributesForVectors(numberOfShownAttributes, quartilsDispersionsList.get(0));
        calculator.printTable(objResult36, bestAttributes, "QuartilsDispersions per Cluster: MinMax + Quartile ", "First",
                "Second");




        //Quartilsdispersion of all points - Quartilsdispersion of cluster
        List<double[][]> matrixList71 = calculator.calculateMatrixList(newClusters);
        List<RealMatrix> quartDif = new ArrayList<>();
        for(int l = 0; l < matrixList71.size(); l++){
            RealMatrix tempQuartDif = calculator.calculateQuartilsDispersion(matrixList71.get(l));
            for(int k = 0; k < tempQuartDif.getRowDimension(); k++){
                tempQuartDif.setEntry(k,0, -(Math.abs(quartilsDispersion.getEntry(k,0)) - Math.abs(tempQuartDif.getEntry(k,0))));
            }
            quartDif.add(tempQuartDif);

        }
        Object[][] objResult71 = calculator.calculateMinimizingTablePerClusterWithVector(newClusters, numberOfShownAttributes, quartDif);
        bestAttributes = calculator.calculateMinAttributesForVectors(numberOfShownAttributes, quartDif.get(0));
        calculator.printTable(objResult71, bestAttributes, "Quartilsdispersion difference per Cluster: MinMax + Quartile ", "First",
                "Second");
        */



        //System.out.println("Gauss: ");
        //evaluateGaussDataset(numberOfPointsPerCluster,numberOfShownAttributes, 10);
        //System.out.println();
        //System.out.println("Equally: ");
        //evaluateEquallyDataset(numberOfPointsPerCluster, numberOfShownAttributes, 10);
        //System.out.println("Big Dataset Equally: ");
        /*try{
            if(args.length < 4){
                evaluateBigDataset('n',5, numberOfPointsPerCluster, numberOfShownAttributes, 100);
            }else {
                evaluateBigDataset(args[0].charAt(0), Integer.parseInt(args[1]), Integer.parseInt(args[2]), numberOfShownAttributes, Integer.parseInt(args[3]));
            }
        }catch(Exception e){
            e.printStackTrace();
        }*/










        /*
        try{
            if(args.length < 4){
                evaluation.evaluateClusteredDataset('n',5, numberOfPointsPerCluster, numberOfShownAttributes, 100);
            }else {
                evaluation.evaluateClusteredDataset(args[0].charAt(0), Integer.parseInt(args[1]), Integer.parseInt(args[2]), numberOfShownAttributes, Integer.parseInt(args[3]));
            }
        }catch(Exception e){
            e.printStackTrace();
        }
        */


        /*List<List<Point>> newClusters2 = new generateKMeansData('n', 0.3, numberOfClusters, numberOfPointsPerCluster).getClusters();
        for (List<Point> p : newClusters2) {
            for (Point point : p) {
                for(int k = 0; k < point.getNumberOfAttributes(); k++){
                    System.out.print(point.getAttributes()[k] + ", ");
                }
                System.out.println();
            }
        }*/



        List<Point> tempClusters = new ArrayList<>();
        System.out.println("Start reading");
        readCSV("data_gaussian_n100000_features10_k5_noise33.csv", tempClusters);
        System.out.println("Reading done");
        List<CentroidCluster<Point>> clusteringResults = new KMeansPlusPlusClusterer<Point>( 5, 100 ).cluster( tempClusters );
        System.out.println("KMeans done");
        List<List<Point>> clusters2 = new ArrayList<>();
        for(int i = 0; i < 5; i++){
            clusters2.add(new ArrayList<Point>());
        }
        for (int k = 0; k < clusteringResults.size(); k++){
            for(int l = 0; l < clusteringResults.get(k).getPoints().size(); l++){
                clusters2.get(k).add(clusteringResults.get(k).getPoints().get(l));
            }
        }

        /*for (List<Point> p : clusters2) {
            for (Point point : p) {
                for(int k = 0; k < point.getNumberOfAttributes(); k++){
                    System.out.print(point.getAttributes()[k] + ", ");
                }
                System.out.println();
            }
        }*/
        evaluation.evaluateDataset(5, clusters2, 1000000, 10, 33);




        //evaluation.evaluatePythonDataset('n',5, numberOfPointsPerCluster, numberOfShownAttributes, 100);


        /*
        //Calculates first the best attributes per cluster and then takes the most frequent of this list
        bestAttributes = calculator.bestAttributesFirstForClusterThenForAll(newClusters, numberOfShownAttributes);
        Object[][] objResult4 = calculator.calculateTable(calculator.calculateMinMaxResults(newClusters, bestAttributes),
                calculator.calculateQuartileResult(newClusters, bestAttributes));
        calculator.printTable(objResult4, bestAttributes, "Spearman first per Cluster then for all: MinMax + Quartile ",
                "First", "Second");
        System.out.println("First per Cluster, then the most frequent: ");
        calculator.calculateOverlapInterval(newClusters, bestAttributes);
        */

    }

}
