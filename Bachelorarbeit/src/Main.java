import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.stat.correlation.Covariance;
import org.apache.commons.math3.stat.correlation.PearsonsCorrelation;
import org.apache.commons.math3.stat.correlation.SpearmansCorrelation;

import java.util.ArrayList;
import java.util.List;

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
        int bereich = (int)(clusters.size() * streuungsreichweite * 10);
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
            System.out.println("Cluster: " + i + ": [" + min + "," + max + "]");
        }


    }

    public static void main(String[] args) {
        List<List<NewPair>> clusters = new ArrayList<>();
        List<List<Point>> newClusters = new ArrayList<>();
        int numberOfClusters = 10;
        int numberOfPointsPerCluster = 100;
        int numberOfAttributes = 10;

        for(int i = 0; i < numberOfClusters; i++){
            clusters.add(new ArrayList<NewPair>());
            newClusters.add(new ArrayList<Point>());

        }
        //generateFakeData(clusters, numberOfPointsPerCluster);
        //calculationFor2Attributes(clusters, numberOfPointsPerCluster);
        generateFakeData(newClusters, numberOfPointsPerCluster, numberOfAttributes);
        /*for(int i = 0; i < clusters.size(); i++){
            for(int j = 0; j < clusters.get(i).size(); j++){
                newClusters.get(i).add(new Point(2, clusters.get(i).get(j).getAttribute()));
            }
        }*/
        calculationForMoreAttributes(newClusters, numberOfPointsPerCluster);
    }

}
