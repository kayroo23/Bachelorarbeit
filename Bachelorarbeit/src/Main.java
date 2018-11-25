import org.apache.commons.math3.distribution.AbstractRealDistribution;
import org.apache.commons.math3.distribution.CauchyDistribution;
import org.apache.commons.math3.distribution.NormalDistribution;
import org.apache.commons.math3.distribution.UniformRealDistribution;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.stat.correlation.Covariance;
import org.apache.commons.math3.stat.correlation.KendallsCorrelation;
import org.apache.commons.math3.stat.correlation.PearsonsCorrelation;
import org.apache.commons.math3.stat.correlation.SpearmansCorrelation;
import org.jfree.ui.RefineryUtilities;

import java.io.File;
import java.io.FileNotFoundException;
import java.util.ArrayList;
import java.util.List;
import java.util.Random;
import java.util.Scanner;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.PrintWriter;

public class Main {

    private static String[] nameOfStatistic = {"Varianz", "StdAbw", "GeomMean", "MedianDev", "Quartilsdisp", "VariationsCoeff",
            "VarianzDiff", "StdAbwDiff", "MedianDevDiff", "VariationsCoeffDiff", "QuartilsdispDiff", "GeomMeanDiff"};


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


    private static List<List<Point>> generateEquallyGoldStandardDataset(double overlap, int pointsPerCluster){
        List<List<Point>> clusters = new ArrayList<>();
        for(int i = 0; i < 5; i++){
            clusters.add(new ArrayList<Point>());
        }
        double overlap1 = overlap * Math.random();
        double overlap2 = overlap * Math.random();
        for(int k = 0; k < pointsPerCluster; k++){
            //alter,gehalt,groesse,kontonr,blutgrp,geschlecht, clusterkonstante, konstante
            double[] attList = new double[8];
            attList[0] = (int)((Math.random()*19)+1);
            attList[1] = (int)((Math.random()*9* (1 + overlap1))+1 * (1 - overlap2));
            attList[2] = (int)((Math.random()*39 * (1 + overlap1))+1 * (1 - overlap2));
            attList[3] = (int)((Math.random()*99)+1);
            attList[4] = ((int)(Math.random()*4))*25;
            attList[5] = Math.random() > 0.5 ? 100 : 0;
            attList[6] = 1;
            attList[7] = 1234567;
            Point p = new Point(attList.length, attList);
            clusters.get(0).add(p);
        }
        overlap1 = overlap * Math.random();
        overlap2 = overlap * Math.random();
        for(int k = 0; k < pointsPerCluster; k++){
            //alter,gehalt,groesse,kontonr,blutgrp,geschlecht, clusterkonstante, konstante
            double[] attList = new double[8];
            attList[0] = (int)((Math.random()*19)+21);
            attList[1] = (int)((Math.random()*39 * (1 + overlap1))+11 * (1 - overlap2));
            attList[2] = (int)((Math.random()*19 * (1 + overlap1))+41 * (1 - overlap2));
            attList[3] = (int)((Math.random()*99)+1);
            attList[4] = ((int)(Math.random()*4))*25;
            attList[5] = Math.random() > 0.5 ? 100 : 0;
            attList[6] = 25;
            attList[7] = 1234567;
            Point p = new Point(attList.length, attList);
            clusters.get(1).add(p);
        }
        overlap1 = overlap * Math.random();
        overlap2 = overlap * Math.random();
        for(int k = 0; k < pointsPerCluster; k++){
            //alter,gehalt,groesse,kontonr,blutgrp,geschlecht, clusterkonstante, konstante
            double[] attList = new double[8];
            attList[0] = (int)((Math.random()*19)+41);
            attList[1] = (int)((Math.random()*29 * (1 + overlap1))+51 * (1 - overlap2));
            attList[2] = (int)((Math.random()*19 * (1 + overlap1))+61 * (1 - overlap2));
            attList[3] = (int)((Math.random()*99)+1);
            attList[4] = ((int)(Math.random()*4))*25;
            attList[5] = Math.random() > 0.5 ? 100 : 0;
            attList[6] = 50;
            attList[7] = 1234567;
            Point p = new Point(attList.length, attList);
            clusters.get(2).add(p);
        }
        overlap1 = overlap * Math.random();
        overlap2 = overlap * Math.random();
        for(int k = 0; k < pointsPerCluster; k++){
            //alter,gehalt,groesse,kontonr,blutgrp,geschlecht, clusterkonstante, konstante
            double[] attList = new double[8];
            attList[0] = (int)((Math.random()*19)+61);
            attList[1] = (int)((Math.random()*15 * (1 + overlap1))+81 * (1 - overlap2));
            attList[2] = (int)((Math.random()*12 * (1 + overlap1))+81 * (1 - overlap2));
            attList[3] = (int)((Math.random()*99)+1);
            attList[4] = ((int)(Math.random()*4))*25;
            attList[5] = Math.random() > 0.5 ? 100 : 0;
            attList[6] = 75;
            attList[7] = 1234567;
            Point p = new Point(attList.length, attList);
            clusters.get(3).add(p);
        }
        overlap1 = overlap * Math.random();
        overlap2 = overlap * Math.random();
        for(int k = 0; k < pointsPerCluster; k++){
            //alter,gehalt,groesse,kontonr,blutgrp,geschlecht, clusterkonstante, konstante
            double[] attList = new double[8];
            attList[0] = (int)((Math.random()*19)+81);
            attList[1] = (int)((Math.random()*4 * (1 + overlap1))+96 * (1 - overlap2));
            attList[2] = (int)((Math.random()*6 * (1 + overlap1))+94 * (1 - overlap2));
            attList[3] = (int)((Math.random()*99)+1);
            attList[4] = ((int)(Math.random()*4))*25;
            attList[5] = Math.random() > 0.5 ? 100 : 0;
            attList[6] = 100;
            attList[7] = 1234567;
            Point p = new Point(attList.length, attList);
            clusters.get(4).add(p);
        }

        return clusters;
    }

    /**
     *
     * @param overlap percentage overlap
     * @param pointsPerCluster number of points per cluster
     * @return goldstandard dataset
     */
    private static List<List<Point>> generateGaussGoldStandardDataset(double overlap, int pointsPerCluster){
        List<List<Point>> clusters = new ArrayList<>();
        for(int i = 0; i < 5; i++){
            clusters.add(new ArrayList<Point>());
        }
        Random rnd = new Random();
        double overlap1 = overlap * Math.random();
        double overlap2 = overlap * Math.random();
        for(int k = 0; k < pointsPerCluster; k++){
            //alter,gehalt,groesse,kontonr,blutgrp,geschlecht, clusterkonstante, konstante
            double[] attList = new double[8];
            attList[0] = (int)((Math.random()*19)+1);
            attList[1] = calculateGaussian(1,10,overlap1,rnd);
            attList[2] = calculateGaussian(1,40,overlap2,rnd);
            attList[3] = (int)((Math.random()*99)+1);
            attList[4] = ((int)(Math.random()*4))*25;
            attList[5] = Math.random() > 0.5 ? 100 : 0;
            attList[6] = 1;
            attList[7] = 1234567;
            Point p = new Point(attList.length, attList);
            clusters.get(0).add(p);
        }
        rnd = new Random();
        overlap1 = overlap * Math.random();
        overlap2 = overlap * Math.random();
        for(int k = 0; k < pointsPerCluster; k++){
            //alter,gehalt,groesse,kontonr,blutgrp,geschlecht, clusterkonstante, konstante
            double[] attList = new double[8];
            attList[0] = (int)((Math.random()*19)+21);
            attList[1] = calculateGaussian(11,50,overlap1,rnd);
            attList[2] = calculateGaussian(41,60,overlap2,rnd);
            attList[3] = (int)((Math.random()*99)+1);
            attList[4] = ((int)(Math.random()*4))*25;
            attList[5] = Math.random() > 0.5 ? 100 : 0;
            attList[6] = 25;
            attList[7] = 1234567;
            Point p = new Point(attList.length, attList);
            clusters.get(1).add(p);
        }
        rnd = new Random();
        overlap1 = overlap * Math.random();
        overlap2 = overlap * Math.random();
        for(int k = 0; k < pointsPerCluster; k++){
            //alter,gehalt,groesse,kontonr,blutgrp,geschlecht, clusterkonstante, konstante
            double[] attList = new double[8];
            attList[0] = (int)((Math.random()*19)+41);
            attList[1] = calculateGaussian(51,80,overlap1,rnd);
            attList[2] = calculateGaussian(61,80,overlap2,rnd);
            attList[3] = (int)((Math.random()*99)+1);
            attList[4] = ((int)(Math.random()*4))*25;
            attList[5] = Math.random() > 0.5 ? 100 : 0;
            attList[6] = 50;
            attList[7] = 1234567;
            Point p = new Point(attList.length, attList);
            clusters.get(2).add(p);
        }
        rnd = new Random();
        overlap1 = overlap * Math.random();
        overlap2 = overlap * Math.random();
        for(int k = 0; k < pointsPerCluster; k++){
            //alter,gehalt,groesse,kontonr,blutgrp,geschlecht, clusterkonstante, konstante
            double[] attList = new double[8];
            attList[0] = (int)((Math.random()*19)+61);
            attList[1] = calculateGaussian(81,95,overlap1,rnd);
            attList[2] = calculateGaussian(81,93,overlap2,rnd);
            attList[3] = (int)((Math.random()*99)+1);
            attList[4] = ((int)(Math.random()*4))*25;
            attList[5] = Math.random() > 0.5 ? 100 : 0;
            attList[6] = 75;
            attList[7] = 1234567;
            Point p = new Point(attList.length, attList);
            clusters.get(3).add(p);
        }
        rnd = new Random();
        overlap1 = overlap * Math.random();
        overlap2 = overlap * Math.random();
        for(int k = 0; k < pointsPerCluster; k++){
            //alter,gehalt,groesse,kontonr,blutgrp,geschlecht, clusterkonstante, konstante
            double[] attList = new double[8];
            attList[0] = (int)((Math.random()*19)+81);
            attList[1] = calculateGaussian(95,100,overlap1,rnd);
            attList[2] = calculateGaussian(94,100,overlap2,rnd);
            attList[3] = (int)((Math.random()*99)+1);
            attList[4] = ((int)(Math.random()*4))*25;
            attList[5] = Math.random() > 0.5 ? 100 : 0;
            attList[6] = 100;
            attList[7] = 1234567;
            Point p = new Point(attList.length, attList);
            clusters.get(4).add(p);
        }

        return clusters;
    }

    /**Generates the y Value of an logistic function
     *
     * @param x Zähler
     * @param maxX Nenner
     * @return y Value
     */
    private static double getLogisticFunctionValue(double x, double maxX){
        double t= -2 + (5*(x/maxX));
        return (1/(1+Math.exp(-t)));
    }

    /**Generates a fakeData goldtandard
     *
     *
     * @param c     char to choose the distribution type
     * @param overlap   noisecoefficient
     * @param numberOfClusters number of clusters
     * @param pointsPerCluster points per cluster
     * @return a List of Clusters with fake Data
     */
    private static List<List<Point>> generateGoldStandardMoreClusters(char c, double overlap, int numberOfClusters, int pointsPerCluster){
        List<List<Point>> clusters = new ArrayList<>();
        AbstractRealDistribution randomGenAge;
        AbstractRealDistribution randomGen1;
        AbstractRealDistribution randomGen2;


        for(int i = 0; i < numberOfClusters; i++){
            clusters.add(new ArrayList<Point>());
        }
        int maxValue = numberOfClusters * 10;
        for(int clusterNr = 0; clusterNr < numberOfClusters; clusterNr++){
            Random randomNr = new Random();
            double overlap1 = 1 + (overlap * Math.random());
            double min = (10*clusterNr) - (overlap1*10);
            double max = 10 + (overlap1*10)+(10*clusterNr);

            if(c == 'u'){
                randomGen1 = new UniformRealDistribution(min, max);
                randomGen2 = new UniformRealDistribution(getLogisticFunctionValue(clusterNr, numberOfClusters) * (10*clusterNr),
                        getLogisticFunctionValue(clusterNr+1, numberOfClusters) * (10*(clusterNr+1)));
                randomGenAge = new UniformRealDistribution(10*clusterNr, (10*clusterNr)+10);
            }else if(c == 'c'){
                //je niedriger scale, desto unwarscheinlicher sind abweichungen vom median
                randomGen1 = new CauchyDistribution((min+max)/2,(max-min)/100);
                randomGen2 = new CauchyDistribution(getLogisticFunctionValue(clusterNr + 0.5, numberOfClusters) *  (10*clusterNr),
                        ((getLogisticFunctionValue(clusterNr +1 , numberOfClusters) -
                                getLogisticFunctionValue(clusterNr, numberOfClusters))/100) *  (10*(clusterNr + 1)));
                randomGenAge = new CauchyDistribution((min+max)/2,(max-min)/100);
            }else {
                randomGen1 = new NormalDistribution((min+max)/2, (max-min)/3);
                randomGen2 = new NormalDistribution(getLogisticFunctionValue(clusterNr + 0.5, numberOfClusters) *  (10*clusterNr),
                        ((getLogisticFunctionValue(clusterNr +1 , numberOfClusters) -
                                getLogisticFunctionValue(clusterNr, numberOfClusters))/3) *  (10*(clusterNr + 1)));
                randomGenAge = new NormalDistribution((min+max)/2, (max-min)/3);
            }

            for(int k = 0; k < pointsPerCluster; k++){
                double[] attList = new double[8];
                //good attributes
                attList[0] = randomGenAge.sample();
                attList[1] = randomGen1.sample();
                attList[2] = randomGen2.sample();
                attList[6] = 10 * clusterNr;

                //Random attributes
                attList[3] = randomNr.nextInt(maxValue);
                attList[4] = (randomNr.nextInt(4))*((float)numberOfClusters*2.5);
                attList[5] = Math.random() > 0.5 ? maxValue : 0;
                attList[7] = 1234567;

                Point p = new Point(attList.length, attList);
                clusters.get(clusterNr).add(p);
            }
        }
        return clusters;
    }

    /**Calculates gaussian distribution where the max/min values are reached at +distance/-distance
     * The higher the distance the lower the chance to get a value out of the [min,max] interval
     *
     * @param min value after -distance
     * @param max value after +distance
     * @param overlap increases the range
     * @return gaussian distribution
     */
    private static double calculateGaussian(double min, double max, double overlap, Random rnd){
        double localMin = min * (1 - overlap);
        double localMax = max * (1 + overlap);
        double mean = (localMin + localMax) / 2;
        return (mean + rnd.nextGaussian()*((localMax - mean)/3));
    }

    private static void evaluateEquallyDataset(int numberOfPointsPerCluster, int numberOfShownAttributes, int iterations){
        GeneralCalculation calculator = new GeneralCalculation();
        List<Integer> knownBestAttributes = new ArrayList<>();
        knownBestAttributes.add(0);
        knownBestAttributes.add(1);
        knownBestAttributes.add(2);
        knownBestAttributes.add(6);
        double overlap = 0;
        int steps = 10;
        double[][] count = new double[steps+1][12];
        List<List<Point>> newClusters2 = new ArrayList<>();
        List<List<Integer>> listOfBestAtt = new ArrayList<>();
        for(int i = 0; i <= steps; i++){
            overlap = ((double)i)/10;
            for(int k = 0; k < iterations; k++){
                newClusters2 = generateEquallyGoldStandardDataset(overlap, numberOfPointsPerCluster);
                List<double[][]> matrixList = calculator.calculateMatrixList(newClusters2);
                for(int l = 0; l < matrixList.size(); l++){

                    //standard metriken
                    listOfBestAtt.add(calculator.calculateMinAttributesForVectors(numberOfShownAttributes,
                            calculator.calculateVariance(matrixList.get(l))));
                    listOfBestAtt.add(calculator.calculateMinAttributesForVectors(numberOfShownAttributes,
                            calculator.calculateStandardDeviation(matrixList.get(l))));
                    listOfBestAtt.add(calculator.calculateMinAttributesForVectors(numberOfShownAttributes,
                            calculator.calculateGeomMean(matrixList.get(l))));
                    listOfBestAtt.add(calculator.calculateMinAttributesForVectors(numberOfShownAttributes,
                            calculator.calculateMedianDeviation(matrixList.get(l))));
                    listOfBestAtt.add(calculator.calculateMinAttributesForVectors(numberOfShownAttributes,
                            calculator.calculateQuartilsDispersion(matrixList.get(l))));
                    listOfBestAtt.add(calculator.calculateMinAttributesForVectors(numberOfShownAttributes,
                            calculator.calculateVariationsCoefficient(matrixList.get(l))));



                    double[][] matrix = calculator.calculateMatrix(newClusters2);
                    RealMatrix standardDeviation = calculator.calculateStandardDeviation(matrix);
                    RealMatrix geomMean = calculator.calculateGeomMean(matrix);
                    RealMatrix variance = calculator.calculateVariance(matrix);
                    RealMatrix medianDeviation = calculator.calculateMedianDeviation(matrix);
                    RealMatrix variationsCoefficient = calculator.calculateVariationsCoefficient(matrix);
                    RealMatrix quartilsDispersion = calculator.calculateQuartilsDispersion(matrix);

                    List<RealMatrix> varianceList2 = new ArrayList<>();
                    List<RealMatrix> stdabwList = new ArrayList<>();
                    List<RealMatrix> medDevList = new ArrayList<>();
                    List<RealMatrix> varCofList = new ArrayList<>();
                    List<RealMatrix> quartDispList = new ArrayList<>();
                    List<RealMatrix> geomMeanList = new ArrayList<>();

                    //Berechnung für metrikdifferenzen
                    for(int x = 0; x < matrixList.size(); x++){
                        RealMatrix tempVariance = calculator.calculateVariance(matrixList.get(x));
                        RealMatrix tempStdabw = calculator.calculateStandardDeviation(matrixList.get(x));
                        RealMatrix tempMedDev = calculator.calculateMedianDeviation(matrixList.get(x));
                        RealMatrix tempVarCof = calculator.calculateVariationsCoefficient(matrixList.get(x));
                        RealMatrix tempQuartDisp = calculator.calculateQuartilsDispersion(matrixList.get(x));
                        RealMatrix tempGeomMean = calculator.calculateGeomMean(matrixList.get(x));

                        for(int y = 0; y < tempVariance.getRowDimension(); y++){
                            tempVariance.setEntry(y,0, -(Math.abs(variance.getEntry(y,0)) - Math.abs(tempVariance.getEntry(y,0))));
                            tempStdabw.setEntry(y,0, -(Math.abs(standardDeviation.getEntry(y,0)) - Math.abs(tempStdabw.getEntry(y,0))));
                            tempMedDev.setEntry(y,0, -(Math.abs(medianDeviation.getEntry(y,0)) - Math.abs(tempMedDev.getEntry(y,0))));
                            tempVarCof.setEntry(y,0, -(Math.abs(variationsCoefficient.getEntry(y,0)) - Math.abs(tempVarCof.getEntry(y,0))));
                            tempQuartDisp.setEntry(y,0, -(Math.abs(quartilsDispersion.getEntry(y,0)) - Math.abs(tempQuartDisp.getEntry(y,0))));
                            tempGeomMean.setEntry(y,0, -(Math.abs(geomMean.getEntry(y,0)) - Math.abs(tempGeomMean.getEntry(y,0))));
                        }
                        varianceList2.add(tempVariance);
                        stdabwList.add(tempStdabw);
                        medDevList.add(tempMedDev);
                        varCofList.add(tempVarCof);
                        quartDispList.add(tempQuartDisp);
                        geomMeanList.add(tempGeomMean);
                    }
                    listOfBestAtt.add(calculator.calculateMinAttributesForVectors(numberOfShownAttributes, varianceList2.get(l)));
                    listOfBestAtt.add(calculator.calculateMinAttributesForVectors(numberOfShownAttributes, stdabwList.get(l)));
                    listOfBestAtt.add(calculator.calculateMinAttributesForVectors(numberOfShownAttributes, medDevList.get(l)));
                    listOfBestAtt.add(calculator.calculateMinAttributesForVectors(numberOfShownAttributes, varCofList.get(l)));
                    listOfBestAtt.add(calculator.calculateMinAttributesForVectors(numberOfShownAttributes, quartDispList.get(l)));
                    listOfBestAtt.add(calculator.calculateMinAttributesForVectors(numberOfShownAttributes, geomMeanList.get(l)));


                    for(int u = 0; u < listOfBestAtt.size(); u++){
                        for (int att : listOfBestAtt.get(u)) {
                            if(knownBestAttributes.contains(att)){
                                count[i][u]++;
                            }
                        }
                    }
                    listOfBestAtt.clear();
                }
                matrixList.clear();
            }
            //5 ist die Anzahl an Clustern
            double divisor = 5*numberOfShownAttributes*iterations;
            System.out.print("Rauschwert der Daten: ");
            System.out.println(overlap);
            System.out.print("Genauigkeit variance: ");
            System.out.println(count[i][0]/(divisor));
            System.out.print("Genauigkeit stdabw: ");
            System.out.println(count[i][1]/(divisor));
            System.out.print("Genauigkeit Geometric Mean: ");
            System.out.println(count[i][2]/(divisor));
            System.out.print("Genauigkeit Median Deviation: ");
            System.out.println(count[i][3]/(divisor));
            System.out.print("Genauigkeit Quartilsdispersionscoeff.: ");
            System.out.println(count[i][4]/(divisor));
            System.out.print("Genauigkeit Variationscoeff.: ");
            System.out.println(count[i][5]/(divisor));
            System.out.print("Genauigkeit Varianzdifferenz: ");
            System.out.println(count[i][6]/(divisor));
            System.out.print("Genauigkeit Stdabw.differenz: ");
            System.out.println(count[i][7]/(divisor));
            System.out.print("Genauigkeit Mediandev.differenz: ");
            System.out.println(count[i][8]/(divisor));
            System.out.print("Genauigkeit Variationscoeff.differenz: ");
            System.out.println(count[i][9]/(divisor));
            System.out.print("Genauigkeit Quartilsdisp.differenz: ");
            System.out.println(count[i][10]/(divisor));
            System.out.print("Genauigkeit GeometricMean differenz: ");
            System.out.println(count[i][11]/(divisor));


            System.out.println();
        }
        /*
        //5 ist die Anzahl an Clustern
        double divisor = 5*numberOfShownAttributes*steps*iterations;
        System.out.print("Genauigkeit variance: ");
        System.out.println(count[0]/(divisor));
        System.out.print("Genauigkeit stdabw: ");
        System.out.println(count[1]/(divisor));
        System.out.print("Genauigkeit Geometric Mean: ");
        System.out.println(count[2]/(divisor));
        System.out.print("Genauigkeit Median Deviation: ");
        System.out.println(count[3]/(divisor));
        System.out.print("Genauigkeit Quartilsdispersionscoeff.: ");
        System.out.println(count[4]/(divisor));
        System.out.print("Genauigkeit Variationscoeff.: ");
        System.out.println(count[5]/(divisor));
        System.out.print("Genauigkeit Varianzdifferenz: ");
        System.out.println(count[6]/(divisor));
        System.out.print("Genauigkeit Stdabw.differenz: ");
        System.out.println(count[7]/(divisor));
        System.out.print("Genauigkeit Mediandev.differenz: ");
        System.out.println(count[8]/(divisor));
        System.out.print("Genauigkeit Variationscoeff.differenz: ");
        System.out.println(count[9]/(divisor));
        System.out.print("Genauigkeit Quartilsdisp.differenz: ");
        System.out.println(count[10]/(divisor));
        System.out.print("Genauigkeit GeometricMean differenz: ");
        System.out.println(count[11]/(divisor));
        */

        PlotLineChart chart = new PlotLineChart(
                "Equallydistribution" ,
                "Equallydistribution");
        for(int i = 0; i < count[0].length; i++){
            for(int j = 0; j < count.length; j++){
                chart.addToDataset(count[j][i]/(5*numberOfShownAttributes*iterations),
                        nameOfStatistic[i], Integer.toString(j));
            }
        }

        chart.setFinalData("");
        chart.pack( );
        RefineryUtilities.centerFrameOnScreen( chart );
        chart.setVisible( true );

    }

    private static void evaluateGaussDataset(int numberOfPointsPerCluster, int numberOfShownAttributes, int iterations){
        GeneralCalculation calculator = new GeneralCalculation();
        List<Integer> knownBestAttributes = new ArrayList<>();
        knownBestAttributes.add(0);
        knownBestAttributes.add(1);
        knownBestAttributes.add(2);
        knownBestAttributes.add(6);

        double overlap = 0;
        int steps = 10;
        double[][] count = new double[11][12];
        List<List<Point>> newClusters2 = new ArrayList<>();
        List<List<Integer>> listOfBestAtt = new ArrayList<>();
        for(int i = 0; i <= steps; i++){
            overlap = ((double)i)/10;
            for(int k = 0; k < iterations; k++){
                newClusters2 = generateGaussGoldStandardDataset(overlap, numberOfPointsPerCluster);
                List<double[][]> matrixList = calculator.calculateMatrixList(newClusters2);
                for(int l = 0; l < matrixList.size(); l++){

                    //standard metriken
                    listOfBestAtt.add(calculator.calculateMinAttributesForVectors(numberOfShownAttributes,
                            calculator.calculateVariance(matrixList.get(l))));
                    listOfBestAtt.add(calculator.calculateMinAttributesForVectors(numberOfShownAttributes,
                            calculator.calculateStandardDeviation(matrixList.get(l))));
                    listOfBestAtt.add(calculator.calculateMinAttributesForVectors(numberOfShownAttributes,
                            calculator.calculateGeomMean(matrixList.get(l))));
                    listOfBestAtt.add(calculator.calculateMinAttributesForVectors(numberOfShownAttributes,
                            calculator.calculateMedianDeviation(matrixList.get(l))));
                    listOfBestAtt.add(calculator.calculateMinAttributesForVectors(numberOfShownAttributes,
                            calculator.calculateQuartilsDispersion(matrixList.get(l))));
                    listOfBestAtt.add(calculator.calculateMinAttributesForVectors(numberOfShownAttributes,
                            calculator.calculateVariationsCoefficient(matrixList.get(l))));



                    double[][] matrix = calculator.calculateMatrix(newClusters2);
                    RealMatrix standardDeviation = calculator.calculateStandardDeviation(matrix);
                    RealMatrix geomMean = calculator.calculateGeomMean(matrix);
                    RealMatrix variance = calculator.calculateVariance(matrix);
                    RealMatrix medianDeviation = calculator.calculateMedianDeviation(matrix);
                    RealMatrix variationsCoefficient = calculator.calculateVariationsCoefficient(matrix);
                    RealMatrix quartilsDispersion = calculator.calculateQuartilsDispersion(matrix);

                    List<RealMatrix> varianceList2 = new ArrayList<>();
                    List<RealMatrix> stdabwList = new ArrayList<>();
                    List<RealMatrix> medDevList = new ArrayList<>();
                    List<RealMatrix> varCofList = new ArrayList<>();
                    List<RealMatrix> quartDispList = new ArrayList<>();
                    List<RealMatrix> geomMeanList = new ArrayList<>();

                    //Berechnung für metrikdifferenzen
                    for(int x = 0; x < matrixList.size(); x++){
                        RealMatrix tempVariance = calculator.calculateVariance(matrixList.get(x));
                        RealMatrix tempStdabw = calculator.calculateStandardDeviation(matrixList.get(x));
                        RealMatrix tempMedDev = calculator.calculateMedianDeviation(matrixList.get(x));
                        RealMatrix tempVarCof = calculator.calculateVariationsCoefficient(matrixList.get(x));
                        RealMatrix tempQuartDisp = calculator.calculateQuartilsDispersion(matrixList.get(x));
                        RealMatrix tempGeomMean = calculator.calculateGeomMean(matrixList.get(x));

                        for(int y = 0; y < tempVariance.getRowDimension(); y++){
                            tempVariance.setEntry(y,0, -(Math.abs(variance.getEntry(y,0)) - Math.abs(tempVariance.getEntry(y,0))));
                            tempStdabw.setEntry(y,0, -(Math.abs(standardDeviation.getEntry(y,0)) - Math.abs(tempStdabw.getEntry(y,0))));
                            tempMedDev.setEntry(y,0, -(Math.abs(medianDeviation.getEntry(y,0)) - Math.abs(tempMedDev.getEntry(y,0))));
                            tempVarCof.setEntry(y,0, -(Math.abs(variationsCoefficient.getEntry(y,0)) - Math.abs(tempVarCof.getEntry(y,0))));
                            tempQuartDisp.setEntry(y,0, -(Math.abs(quartilsDispersion.getEntry(y,0)) - Math.abs(tempQuartDisp.getEntry(y,0))));
                            tempGeomMean.setEntry(y,0, -(Math.abs(geomMean.getEntry(y,0)) - Math.abs(tempGeomMean.getEntry(y,0))));
                        }
                        varianceList2.add(tempVariance);
                        stdabwList.add(tempStdabw);
                        medDevList.add(tempMedDev);
                        varCofList.add(tempVarCof);
                        quartDispList.add(tempQuartDisp);
                        geomMeanList.add(tempGeomMean);
                    }
                    listOfBestAtt.add(calculator.calculateMinAttributesForVectors(numberOfShownAttributes, varianceList2.get(l)));
                    listOfBestAtt.add(calculator.calculateMinAttributesForVectors(numberOfShownAttributes, stdabwList.get(l)));
                    listOfBestAtt.add(calculator.calculateMinAttributesForVectors(numberOfShownAttributes, medDevList.get(l)));
                    listOfBestAtt.add(calculator.calculateMinAttributesForVectors(numberOfShownAttributes, varCofList.get(l)));
                    listOfBestAtt.add(calculator.calculateMinAttributesForVectors(numberOfShownAttributes, quartDispList.get(l)));
                    listOfBestAtt.add(calculator.calculateMinAttributesForVectors(numberOfShownAttributes, geomMeanList.get(l)));


                    for(int u = 0; u < listOfBestAtt.size(); u++){
                        for (int att : listOfBestAtt.get(u)) {
                            if(knownBestAttributes.contains(att)){
                                count[i][u]++;
                            }
                        }
                    }
                    listOfBestAtt.clear();
                }
                matrixList.clear();
            }
            //5 ist die Anzahl an Clustern
            double divisor = 5*numberOfShownAttributes*iterations;
            System.out.print("Rauschwert der Daten: ");
            System.out.println(overlap);
            System.out.print("Genauigkeit variance: ");
            System.out.println(count[i][0]/(divisor));
            System.out.print("Genauigkeit stdabw: ");
            System.out.println(count[i][1]/(divisor));
            System.out.print("Genauigkeit Geometric Mean: ");
            System.out.println(count[i][2]/(divisor));
            System.out.print("Genauigkeit Median Deviation: ");
            System.out.println(count[i][3]/(divisor));
            System.out.print("Genauigkeit Quartilsdispersionscoeff.: ");
            System.out.println(count[i][4]/(divisor));
            System.out.print("Genauigkeit Variationscoeff.: ");
            System.out.println(count[i][5]/(divisor));
            System.out.print("Genauigkeit Varianzdifferenz: ");
            System.out.println(count[i][6]/(divisor));
            System.out.print("Genauigkeit Stdabw.differenz: ");
            System.out.println(count[i][7]/(divisor));
            System.out.print("Genauigkeit Mediandev.differenz: ");
            System.out.println(count[i][8]/(divisor));
            System.out.print("Genauigkeit Variationscoeff.differenz: ");
            System.out.println(count[i][9]/(divisor));
            System.out.print("Genauigkeit Quartilsdisp.differenz: ");
            System.out.println(count[i][10]/(divisor));
            System.out.print("Genauigkeit GeometricMean differenz: ");
            System.out.println(count[i][11]/(divisor));

            System.out.println();
        }
        /*
        //5 ist die Anzahl an Clustern
        double divisor = 5*numberOfShownAttributes*steps*iterations;
        System.out.print("Genauigkeit variance: ");
        System.out.println(count[0]/(divisor));
        System.out.print("Genauigkeit stdabw: ");
        System.out.println(count[1]/(divisor));
        System.out.print("Genauigkeit Geometric Mean: ");
        System.out.println(count[2]/(divisor));
        System.out.print("Genauigkeit Median Deviation: ");
        System.out.println(count[3]/(divisor));
        System.out.print("Genauigkeit Quartilsdispersionscoeff.: ");
        System.out.println(count[4]/(divisor));
        System.out.print("Genauigkeit Variationscoeff.: ");
        System.out.println(count[5]/(divisor));
        System.out.print("Genauigkeit Varianzdifferenz: ");
        System.out.println(count[6]/(divisor));
        System.out.print("Genauigkeit Stdabw.differenz: ");
        System.out.println(count[7]/(divisor));
        System.out.print("Genauigkeit Mediandev.differenz: ");
        System.out.println(count[8]/(divisor));
        System.out.print("Genauigkeit Variationscoeff.differenz: ");
        System.out.println(count[9]/(divisor));
        System.out.print("Genauigkeit Quartilsdisp.differenz: ");
        System.out.println(count[10]/(divisor));
        System.out.print("Genauigkeit GeometricMean differenz: ");
        System.out.println(count[11]/(divisor));
        */

        PlotLineChart chart = new PlotLineChart(
                "Gaussdistribution" ,
                "Gaussdistribution");
        for(int i = 0; i < count[0].length; i++){
            for(int j = 0; j < count.length; j++){
                chart.addToDataset(count[j][i]/(5*numberOfShownAttributes*iterations),
                        nameOfStatistic[i], Integer.toString(j));
            }
        }

        chart.setFinalData("");
        chart.pack( );
        RefineryUtilities.centerFrameOnScreen( chart );
        chart.setVisible( true );
    }

    private static void evaluateBigDataset(char distribution, int numberOfClusters, int numberOfPointsPerCluster, int numberOfShownAttributes, int iterations)
            throws FileNotFoundException{
        GeneralCalculation calculator = new GeneralCalculation();
        List<Integer> knownBestAttributes = new ArrayList<>();

        knownBestAttributes.add(0);
        knownBestAttributes.add(1);
        knownBestAttributes.add(2);
        knownBestAttributes.add(6);
        double overlap;
        int steps = 10;
        double[][] count = new double[steps+1][12];
        long[] time = new long[12];

        List<List<Point>> newClusters2;
        List<List<Integer>> listOfBestAtt = new ArrayList<>();
        for(int i = 0; i <= steps; i++){
            overlap = ((double)i)/10;
            for(int k = 0; k < iterations; k++){
                newClusters2 = generateGoldStandardMoreClusters(distribution, overlap, numberOfClusters, numberOfPointsPerCluster);
                List<double[][]> matrixList = calculator.calculateMatrixList(newClusters2);

                long tempTime;
                double[][] matrix = calculator.calculateMatrix(newClusters2);
                tempTime = System.currentTimeMillis();
                RealMatrix standardDeviation = calculator.calculateStandardDeviation(matrix);
                time[6] += System.currentTimeMillis() - tempTime;
                tempTime = System.currentTimeMillis();
                RealMatrix geomMean = calculator.calculateGeomMean(matrix);
                time[7] += System.currentTimeMillis() - tempTime;
                tempTime = System.currentTimeMillis();
                RealMatrix variance = calculator.calculateVariance(matrix);
                time[8] += System.currentTimeMillis() - tempTime;
                tempTime = System.currentTimeMillis();
                RealMatrix medianDeviation = calculator.calculateMedianDeviation(matrix);
                time[9] += System.currentTimeMillis() - tempTime;
                tempTime = System.currentTimeMillis();
                RealMatrix variationsCoefficient = calculator.calculateVariationsCoefficient(matrix);
                time[10] += System.currentTimeMillis() - tempTime;
                tempTime = System.currentTimeMillis();
                RealMatrix quartilsDispersion = calculator.calculateQuartilsDispersion(matrix);
                time[11] += System.currentTimeMillis() - tempTime;


                for(int l = 0; l < matrixList.size(); l++){
                    //standard metriken
                    tempTime = System.currentTimeMillis();
                    listOfBestAtt.add(calculator.calculateMinAttributesForVectors(numberOfShownAttributes,
                            calculator.calculateVariance(matrixList.get(l))));
                    time[0] += System.currentTimeMillis() - tempTime;
                    tempTime = System.currentTimeMillis();
                    listOfBestAtt.add(calculator.calculateMinAttributesForVectors(numberOfShownAttributes,
                            calculator.calculateStandardDeviation(matrixList.get(l))));
                    time[1] += System.currentTimeMillis() - tempTime;
                    tempTime = System.currentTimeMillis();
                    listOfBestAtt.add(calculator.calculateMinAttributesForVectors(numberOfShownAttributes,
                            calculator.calculateGeomMean(matrixList.get(l))));
                    time[2] += System.currentTimeMillis() - tempTime;
                    tempTime = System.currentTimeMillis();
                    listOfBestAtt.add(calculator.calculateMinAttributesForVectors(numberOfShownAttributes,
                            calculator.calculateMedianDeviation(matrixList.get(l))));
                    time[3] += System.currentTimeMillis() - tempTime;
                    tempTime = System.currentTimeMillis();
                    listOfBestAtt.add(calculator.calculateMinAttributesForVectors(numberOfShownAttributes,
                            calculator.calculateQuartilsDispersion(matrixList.get(l))));
                    time[4] += System.currentTimeMillis() - tempTime;
                    tempTime = System.currentTimeMillis();
                    listOfBestAtt.add(calculator.calculateMinAttributesForVectors(numberOfShownAttributes,
                            calculator.calculateVariationsCoefficient(matrixList.get(l))));
                    time[5] += System.currentTimeMillis() - tempTime;




                    List<RealMatrix> varianceList2 = new ArrayList<>();
                    List<RealMatrix> stdabwList = new ArrayList<>();
                    List<RealMatrix> medDevList = new ArrayList<>();
                    List<RealMatrix> varCofList = new ArrayList<>();
                    List<RealMatrix> quartDispList = new ArrayList<>();
                    List<RealMatrix> geomMeanList = new ArrayList<>();

                    //Berechnung für metrikdifferenzen
                    for(int x = 0; x < matrixList.size(); x++){

                        tempTime = System.currentTimeMillis();
                        RealMatrix tempVariance = calculator.calculateVariance(matrixList.get(x));
                        time[6] += System.currentTimeMillis() - tempTime;
                        tempTime = System.currentTimeMillis();
                        RealMatrix tempStdabw = calculator.calculateStandardDeviation(matrixList.get(x));
                        time[7] += System.currentTimeMillis() - tempTime;
                        tempTime = System.currentTimeMillis();
                        RealMatrix tempMedDev = calculator.calculateMedianDeviation(matrixList.get(x));
                        time[8] += System.currentTimeMillis() - tempTime;
                        tempTime = System.currentTimeMillis();
                        RealMatrix tempVarCof = calculator.calculateVariationsCoefficient(matrixList.get(x));
                        time[9] += System.currentTimeMillis() - tempTime;
                        tempTime = System.currentTimeMillis();
                        RealMatrix tempQuartDisp = calculator.calculateQuartilsDispersion(matrixList.get(x));
                        time[10] += System.currentTimeMillis() - tempTime;
                        tempTime = System.currentTimeMillis();
                        RealMatrix tempGeomMean = calculator.calculateGeomMean(matrixList.get(x));
                        time[11] += System.currentTimeMillis() - tempTime;

                        for(int y = 0; y < tempVariance.getRowDimension(); y++){
                            tempTime = System.currentTimeMillis();
                            tempVariance.setEntry(y,0, -(Math.abs(variance.getEntry(y,0)) - Math.abs(tempVariance.getEntry(y,0))));
                            time[6] += System.currentTimeMillis() - tempTime;
                            tempTime = System.currentTimeMillis();
                            tempStdabw.setEntry(y,0, -(Math.abs(standardDeviation.getEntry(y,0)) - Math.abs(tempStdabw.getEntry(y,0))));
                            time[7] += System.currentTimeMillis() - tempTime;
                            tempTime = System.currentTimeMillis();
                            tempMedDev.setEntry(y,0, -(Math.abs(medianDeviation.getEntry(y,0)) - Math.abs(tempMedDev.getEntry(y,0))));
                            time[8] += System.currentTimeMillis() - tempTime;
                            tempTime = System.currentTimeMillis();
                            tempVarCof.setEntry(y,0, -(Math.abs(variationsCoefficient.getEntry(y,0)) - Math.abs(tempVarCof.getEntry(y,0))));
                            time[9] += System.currentTimeMillis() - tempTime;
                            tempTime = System.currentTimeMillis();
                            tempQuartDisp.setEntry(y,0, -(Math.abs(quartilsDispersion.getEntry(y,0)) - Math.abs(tempQuartDisp.getEntry(y,0))));
                            time[10] += System.currentTimeMillis() - tempTime;
                            tempTime = System.currentTimeMillis();
                            tempGeomMean.setEntry(y,0, -(Math.abs(geomMean.getEntry(y,0)) - Math.abs(tempGeomMean.getEntry(y,0))));
                            time[11] += System.currentTimeMillis() - tempTime;
                        }
                        varianceList2.add(tempVariance);
                        stdabwList.add(tempStdabw);
                        medDevList.add(tempMedDev);
                        varCofList.add(tempVarCof);
                        quartDispList.add(tempQuartDisp);
                        geomMeanList.add(tempGeomMean);
                    }

                    tempTime = System.currentTimeMillis();
                    listOfBestAtt.add(calculator.calculateMinAttributesForVectors(numberOfShownAttributes, varianceList2.get(l)));
                    time[6] += System.currentTimeMillis() - tempTime;
                    tempTime = System.currentTimeMillis();
                    listOfBestAtt.add(calculator.calculateMinAttributesForVectors(numberOfShownAttributes, stdabwList.get(l)));
                    time[7] += System.currentTimeMillis() - tempTime;
                    tempTime = System.currentTimeMillis();
                    listOfBestAtt.add(calculator.calculateMinAttributesForVectors(numberOfShownAttributes, medDevList.get(l)));
                    time[8] += System.currentTimeMillis() - tempTime;
                    tempTime = System.currentTimeMillis();
                    listOfBestAtt.add(calculator.calculateMinAttributesForVectors(numberOfShownAttributes, varCofList.get(l)));
                    time[9] += System.currentTimeMillis() - tempTime;
                    tempTime = System.currentTimeMillis();
                    listOfBestAtt.add(calculator.calculateMinAttributesForVectors(numberOfShownAttributes, quartDispList.get(l)));
                    time[10] += System.currentTimeMillis() - tempTime;
                    tempTime = System.currentTimeMillis();
                    listOfBestAtt.add(calculator.calculateMinAttributesForVectors(numberOfShownAttributes, geomMeanList.get(l)));
                    time[11] += System.currentTimeMillis() - tempTime;

                    for(int u = 0; u < listOfBestAtt.size(); u++){
                        for (int att : listOfBestAtt.get(u)) {
                            if(knownBestAttributes.contains(att)){
                                count[i][u]++;
                            }
                        }
                    }
                    listOfBestAtt.clear();
                }
                matrixList.clear();
            }
        }

        double divisor = numberOfClusters*numberOfShownAttributes*iterations;
        StringBuilder sb = new StringBuilder();
        sb.append(distribution);
        sb.append("_");
        sb.append(numberOfClusters);
        sb.append("_");
        sb.append(numberOfPointsPerCluster);
        sb.append("_");
        sb.append(iterations);
        sb.append(".csv");
        PrintWriter csvWriter = new PrintWriter(new File(sb.toString()));
        sb = new StringBuilder();
        sb.append("Distribution");
        sb.append(',');
        sb.append("NumClust");
        sb.append(',');
        sb.append("PointsCluster");
        sb.append(',');
        sb.append("Iterations per Noise");
        sb.append('\n');

        sb.append(distribution);
        sb.append(',');
        sb.append(numberOfClusters);
        sb.append(',');
        sb.append(numberOfPointsPerCluster);
        sb.append(',');
        sb.append(iterations);
        sb.append('\n');

        sb.append("Rauschwert");
        sb.append(',');
        sb.append("var");
        sb.append(',');
        sb.append("stdabw");
        sb.append(',');
        sb.append("geoMean");
        sb.append(',');
        sb.append("MedDev");
        sb.append(',');
        sb.append("QuartDisp");
        sb.append(',');
        sb.append("VariCoeff");
        sb.append(',');
        sb.append("VarDiff");
        sb.append(',');
        sb.append("stdabwDiff");
        sb.append(',');
        sb.append("MedDevDiff");
        sb.append(',');
        sb.append("VarCoffDiff");
        sb.append(',');
        sb.append("QuartDispDiff");
        sb.append(',');
        sb.append("geoMeanDiff");
        sb.append('\n');

        for(int i = 0; i <= steps; i++){
            sb.append(((double)i)/10);
            sb.append(',');
            sb.append(count[i][0]/(divisor));
            sb.append(',');
            sb.append(count[i][1]/(divisor));
            sb.append(',');
            sb.append(count[i][2]/(divisor));
            sb.append(',');
            sb.append(count[i][3]/(divisor));
            sb.append(',');
            sb.append(count[i][4]/(divisor));
            sb.append(',');
            sb.append(count[i][5]/(divisor));
            sb.append(',');
            sb.append(count[i][6]/(divisor));
            sb.append(',');
            sb.append(count[i][7]/(divisor));
            sb.append(',');
            sb.append(count[i][8]/(divisor));
            sb.append(',');
            sb.append(count[i][9]/(divisor));
            sb.append(',');
            sb.append(count[i][10]/(divisor));
            sb.append(',');
            sb.append(count[i][11]/(divisor));
            sb.append('\n');
        }
        sb.append('\n');
        sb.append("Time");
        sb.append('\n');
        sb.append(0);
        sb.append(',');
        for(int i = 0; i < 12; i++){
            sb.append(time[i]);
            sb.append(',');
        }
        csvWriter.write(sb.toString());
        csvWriter.close();
        System.out.println("Done");
    }

    public static void main(String[] args) {
        List<List<NewPair>> clusters = new ArrayList<>();
        List<List<Point>> newClusters = new ArrayList<>();
        int numberOfClusters = 5;
        int numberOfPointsPerCluster = 100;
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
        newClusters = generateGoldStandardMoreClusters('n', 0.5, numberOfClusters, numberOfPointsPerCluster);
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
        try{
            if(args.length < 4){
                evaluateBigDataset('n',5, numberOfPointsPerCluster, numberOfShownAttributes, 100);
            }else {
                evaluateBigDataset(args[0].charAt(0), Integer.parseInt(args[1]), Integer.parseInt(args[2]), numberOfShownAttributes, Integer.parseInt(args[3]));
            }
        }catch(Exception e){
            System.out.println("FileNotFound");
        }




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
