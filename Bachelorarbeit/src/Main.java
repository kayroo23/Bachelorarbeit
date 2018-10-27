import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.stat.correlation.Covariance;
import org.apache.commons.math3.stat.correlation.KendallsCorrelation;
import org.apache.commons.math3.stat.correlation.PearsonsCorrelation;
import org.apache.commons.math3.stat.correlation.SpearmansCorrelation;
import org.apache.commons.math3.stat.descriptive.moment.Mean;

import java.io.File;
import java.io.FileNotFoundException;
import java.util.ArrayList;
import java.util.List;
import java.util.Scanner;

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

    private static void printMatrix(RealMatrix spearmanMatrix, String s){
        System.out.println(s);
        for(int u = 0; u < spearmanMatrix.getRowDimension(); u++){
            System.out.println(spearmanMatrix.getRowVector(u));
        }
        System.out.println();
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


    private static List<List<Point>> generateGoldStandardDataset(){
        List<List<Point>> clusters = new ArrayList<>();
        for(int i = 0; i < 5; i++){
            clusters.add(new ArrayList<Point>());
        }
        for(int k = 0; k < 10000; k++){
            //alter,gehalt,groesse,kontonr,blutgrp,geschlecht, clusterkonstante, konstante
            double[] attList = new double[8];
            attList[0] = (int)((Math.random()*19)+1);
            attList[1] = (int)((Math.random()*9)+1);
            attList[2] = (int)((Math.random()*39)+1);
            attList[3] = (int)((Math.random()*99)+1);
            attList[4] = ((int)(Math.random()*4))*25;
            attList[5] = Math.random() > 0.5 ? 100 : 0;
            attList[6] = 0;
            attList[7] = 1234567;
            Point p = new Point(attList.length, attList);
            clusters.get(0).add(p);
        }
        for(int k = 0; k < 10000; k++){
            //alter,gehalt,groesse,kontonr,blutgrp,geschlecht, clusterkonstante, konstante
            double[] attList = new double[8];
            attList[0] = (int)((Math.random()*19)+21);
            attList[1] = (int)((Math.random()*39)+11);
            attList[2] = (int)((Math.random()*19)+41);
            attList[3] = (int)((Math.random()*99)+1);
            attList[4] = ((int)(Math.random()*4))*25;
            attList[5] = Math.random() > 0.5 ? 100 : 0;
            attList[6] = 1;
            attList[7] = 1234567;
            Point p = new Point(attList.length, attList);
            clusters.get(1).add(p);
        }
        for(int k = 0; k < 10000; k++){
            //alter,gehalt,groesse,kontonr,blutgrp,geschlecht, clusterkonstante, konstante
            double[] attList = new double[8];
            attList[0] = (int)((Math.random()*19)+41);
            attList[1] = (int)((Math.random()*29)+51);
            attList[2] = (int)((Math.random()*19)+61);
            attList[3] = (int)((Math.random()*99)+1);
            attList[4] = ((int)(Math.random()*4))*25;
            attList[5] = Math.random() > 0.5 ? 100 : 0;
            attList[6] = 2;
            attList[7] = 1234567;
            Point p = new Point(attList.length, attList);
            clusters.get(2).add(p);
        }
        for(int k = 0; k < 10000; k++){
            //alter,gehalt,groesse,kontonr,blutgrp,geschlecht, clusterkonstante, konstante
            double[] attList = new double[8];
            attList[0] = (int)((Math.random()*19)+61);
            attList[1] = (int)((Math.random()*5)+81);
            attList[2] = (int)((Math.random()*12)+81);
            attList[3] = (int)((Math.random()*99)+1);
            attList[4] = ((int)(Math.random()*4))*25;
            attList[5] = Math.random() > 0.5 ? 100 : 0;
            attList[6] = 3;
            attList[7] = 1234567;
            Point p = new Point(attList.length, attList);
            clusters.get(3).add(p);
        }
        for(int k = 0; k < 10000; k++){
            //alter,gehalt,groesse,kontonr,blutgrp,geschlecht, clusterkonstante, konstante
            double[] attList = new double[8];
            attList[0] = (int)((Math.random()*19)+81);
            attList[1] = (int)((Math.random()*5)+96);
            attList[2] = (int)((Math.random()*6)+94);
            attList[3] = (int)((Math.random()*99)+1);
            attList[4] = ((int)(Math.random()*4))*25;
            attList[5] = Math.random() > 0.5 ? 100 : 0;
            attList[6] = 4;
            attList[7] = 1234567;
            Point p = new Point(attList.length, attList);
            clusters.get(4).add(p);
        }

        return clusters;
    }

    public static void main(String[] args) {
        List<List<NewPair>> clusters = new ArrayList<>();
        List<List<Point>> newClusters = new ArrayList<>();
        int numberOfClusters = 5;
        int numberOfPointsPerCluster = 500;
        int numberOfAttributes = 10;
        for(int i = 0; i < numberOfClusters; i++){
            clusters.add(new ArrayList<NewPair>());
            newClusters.add(new ArrayList<Point>());
        }
        //generateFakeData(clusters, numberOfPointsPerCluster);
        //calculationFor2Attributes(clusters, numberOfPointsPerCluster);

        //generateFakeData(newClusters, numberOfPointsPerCluster, numberOfAttributes);
        newClusters = generateGoldStandardDataset();
        //calculationForMoreAttributes(newClusters, numberOfPointsPerCluster);

        //readClassificationDataSet("dataset_32_pendigits_changed.txt", newClusters, 0);
        //readClassificationDataSet("irisDataset.txt", newClusters, 0);
        //readClassificationDataSet("fertilityDataSet.txt", newClusters, 0);
        //readClassificationDataSet("ecoliDataSet.txt", newClusters, 1);

        int numberOfShownAttributes = 3;
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
        List<double[][]> matrixList1 = calculator.calculateMatrixList(newClusters);
        List<RealMatrix> spearmanList = new ArrayList<>();
        for(int l = 0; l < matrixList1.size(); l++){
            spearmanList.add(new SpearmansCorrelation().computeCorrelationMatrix(matrixList1.get(l)));
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
        Object[][] objResult31 = calculator.calculateTablePerClusterWithVector(newClusters, numberOfShownAttributes, varianceList);
        bestAttributes = calculator.calculateMinAttributesForVectors(numberOfShownAttributes, varianceList.get(0));
        calculator.printTable(objResult31, bestAttributes, "Variance per Cluster: MinMax + Quartile ", "First",
                "Second");

        List<double[][]> matrixList3 = calculator.calculateMatrixList(newClusters);
        List<RealMatrix> stdabwList = new ArrayList<>();
        for(int l = 0; l < matrixList3.size(); l++){
            stdabwList.add(calculator.calculateStandardDeviation(matrixList3.get(l)));
        }
        Object[][] objResult32 = calculator.calculateTablePerClusterWithVector(newClusters, numberOfShownAttributes, stdabwList);
        bestAttributes = calculator.calculateMinAttributesForVectors(numberOfShownAttributes, stdabwList.get(0));
        calculator.printTable(objResult32, bestAttributes, "Stdabw per Cluster: MinMax + Quartile ", "First",
                "Second");


        List<double[][]> matrixList4 = calculator.calculateMatrixList(newClusters);
        List<RealMatrix> variationskoeffizientList = new ArrayList<>();
        for(int l = 0; l < matrixList4.size(); l++){
            variationskoeffizientList.add(calculator.calculateVariationsCoefficient(matrixList4.get(l)));
        }
        Object[][] objResult33 = calculator.calculateTablePerClusterWithVector(newClusters, numberOfShownAttributes, variationskoeffizientList);
        bestAttributes = calculator.calculateMinAttributesForVectors(numberOfShownAttributes, variationskoeffizientList.get(0));
        calculator.printTable(objResult33, bestAttributes, "Variationskoeffizient per Cluster: MinMax + Quartile ", "First",
                "Second");


        List<double[][]> matrixList5 = calculator.calculateMatrixList(newClusters);
        List<RealMatrix> geometricMeanList = new ArrayList<>();
        for(int l = 0; l < matrixList5.size(); l++){
            geometricMeanList.add(calculator.calculateGeomMean(matrixList5.get(l)));
        }
        Object[][] objResult34 = calculator.calculateTablePerClusterWithVector(newClusters, numberOfShownAttributes, geometricMeanList);
        bestAttributes = calculator.calculateMinAttributesForVectors(numberOfShownAttributes, geometricMeanList.get(0));
        calculator.printTable(objResult34, bestAttributes, "GeometricMean per Cluster: MinMax + Quartile ", "First",
                "Second");


        List<double[][]> matrixList6 = calculator.calculateMatrixList(newClusters);
        List<RealMatrix> medianDeviationList = new ArrayList<>();
        for(int l = 0; l < matrixList6.size(); l++){
            medianDeviationList.add(calculator.calculateMedianDeviation(matrixList6.get(l)));
        }
        Object[][] objResult35 = calculator.calculateTablePerClusterWithVector(newClusters, numberOfShownAttributes, medianDeviationList);
        bestAttributes = calculator.calculateMinAttributesForVectors(numberOfShownAttributes, medianDeviationList.get(0));
        calculator.printTable(objResult35, bestAttributes, "MedianDeviation per Cluster: MinMax + Quartile ", "First",
                "Second");


        List<double[][]> matrixList7 = calculator.calculateMatrixList(newClusters);
        List<RealMatrix> quartilsDispersionsList = new ArrayList<>();
        for(int l = 0; l < matrixList7.size(); l++){
            quartilsDispersionsList.add(calculator.calculateQuartilsDispersion(matrixList7.get(l)));
        }
        Object[][] objResult36 = calculator.calculateTablePerClusterWithVector(newClusters, numberOfShownAttributes, quartilsDispersionsList);
        bestAttributes = calculator.calculateMinAttributesForVectors(numberOfShownAttributes, quartilsDispersionsList.get(0));
        calculator.printTable(objResult36, bestAttributes, "QuartilsDispersions per Cluster: MinMax + Quartile ", "First",
                "Second");








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
