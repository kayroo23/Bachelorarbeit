import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.RealMatrix;
import org.jfree.ui.RefineryUtilities;

import java.io.*;
import java.util.ArrayList;
import java.util.List;
import java.util.Random;

public final class evaluation {

    private static String[] nameOfStatistic = {"Varianz", "StdAbw", "GeomMean", "MedianDev", "Quartilsdisp", "VariationsCoeff",
            "VarianzDiff", "StdAbwDiff", "MedianDevDiff", "VariationsCoeffDiff", "QuartilsdispDiff", "GeomMeanDiff"};

    static void evaluateEquallyDataset(int numberOfPointsPerCluster, int numberOfShownAttributes, int iterations){
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

    static void evaluateGaussDataset(int numberOfPointsPerCluster, int numberOfShownAttributes, int iterations){
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

    static void evaluateBigDataset(char distribution, int numberOfClusters, int numberOfPointsPerCluster, int numberOfShownAttributes, int iterations) {
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
                newClusters2 = new generateSimpleData(distribution, overlap, numberOfClusters, numberOfPointsPerCluster).getClusters();
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
        try{
            generateOutputFile(count, time, distribution, numberOfClusters, numberOfShownAttributes, numberOfPointsPerCluster, iterations, steps);
        }catch(Exception e){
            e.printStackTrace();
        }
    }

    static void evaluatePythonDataset(char distribution, int numberOfClusters, int numberOfPointsPerCluster, int numberOfShownAttributes, int iterations){
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
                String[][] input = new String[(numberOfClusters*numberOfPointsPerCluster) + 1][];
                try{
                    String args = distribution + " " + Double.toString(overlap) + " " + Integer.toString(numberOfClusters) + " " + Integer.toString(numberOfPointsPerCluster);
                    String command = "python ..\\PythonScript\\generateSimpleData.py " + args;
                    Process p = Runtime.getRuntime().exec(command);
                    p.waitFor();
                    BufferedReader br = new BufferedReader(new FileReader("newCsv.csv"));
                    String newLine;
                    int f = 0;
                    while ((newLine = br.readLine()) != null) {
                        input[f] = newLine.split(",");
                        f++;
                    }
                }catch(Exception e){
                    e.printStackTrace();
                }

                newClusters2 = stringsToCluster(input, numberOfClusters);

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
        try{
            generateOutputFile(count, time, distribution, numberOfClusters, numberOfShownAttributes, numberOfPointsPerCluster, iterations, steps);
        }catch(Exception e){
            e.printStackTrace();
        }
    }

    static void evaluateClusteredDataset(char distribution, int numberOfClusters, int numberOfPointsPerCluster, int numberOfShownAttributes, int iterations) {
        GeneralCalculation calculator = new GeneralCalculation();
        List<Integer> knownBestAttributes = new ArrayList<>();
        List<Integer> tempBestAttributes = new ArrayList<>();

        knownBestAttributes.add(0);
        knownBestAttributes.add(1);
        knownBestAttributes.add(2);
        knownBestAttributes.add(6);
        double overlap;
        int steps = 11;
        double[][] count = new double[steps+1][12];
        double[][] notFound = new double[4 * steps][12];
        long[] time = new long[12];

        List<List<Point>> newClusters2;
        List<List<Integer>> listOfBestAtt = new ArrayList<>();
        for(int i = 0; i < steps; i++){
            overlap = ((double)i)/10;
            for(int k = 0; k < iterations; k++){
                newClusters2 = new generateKMeansData(distribution, overlap, numberOfClusters, numberOfPointsPerCluster).getClusters();
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
                        tempBestAttributes = new ArrayList<>(knownBestAttributes);
                        for (int att : listOfBestAtt.get(u)) {
                            if(knownBestAttributes.contains(att)){
                                count[i][u]++;
                                for(int z = 0; z < tempBestAttributes.size(); z++){
                                    if(tempBestAttributes.get(z) == att){
                                        tempBestAttributes.remove(z);
                                        break;
                                    }
                                }
                            }
                        }
                        for(int f = 0; f < tempBestAttributes.size(); f++){
                            if(tempBestAttributes.get(f) == 6){
                                notFound[3 + (i * 4)][u]++;
                            }else{
                                notFound[tempBestAttributes.get(f) + (i * 4)][u]++;
                            }
                        }
                    }
                    listOfBestAtt.clear();
                }
                matrixList.clear();
            }
        }
        try{
            generateOutputFile1(count, notFound, time, distribution, numberOfClusters, numberOfShownAttributes, numberOfPointsPerCluster, iterations, steps);
        }catch(Exception e){
            e.printStackTrace();
        }
    }

    static void evaluateDataset(int numberOfShownAttributes, List<List<Point>> newClusters2, int points, int features, int clusterNr, double noise,
                                double jaccard, double mutual, double rand) {
        GeneralCalculation calculator = new GeneralCalculation();
        List<Integer> knownBestAttributes = new ArrayList<>();
        List<Integer> tempBestAttributes = new ArrayList<>();

        knownBestAttributes.add(0);
        knownBestAttributes.add(1);
        knownBestAttributes.add(2);
        knownBestAttributes.add(3);
        knownBestAttributes.add(4);

        int numberOfMetrics = 10;
        double[] count = new double[numberOfMetrics];
        double[][] notFound = new double[5][numberOfMetrics];
        long[] time = new long[numberOfMetrics];
        int[][] clusterGood = new int[4][newClusters2.size()];

        List<List<double[][]>> listOfBestAtt = new ArrayList<>();
        List<double[][]> matrixList = calculator.calculateMatrixList(newClusters2);

        for(int z = 0; z < numberOfMetrics; z++){
            listOfBestAtt.add(new ArrayList<>());
        }

        long tempTime;
        //double[][] matrix = calculator.calculateMatrix(newClusters2);
        double[][] matrix = calculator.matchClusterPoints(newClusters2);

        tempTime = System.nanoTime();
        RealMatrix variance = calculator.calculateVariance(matrix);
        //RealMatrix m = new Array2DRowRealMatrix(matrix);
        //GeneralCalculation.printMatrix(m, "");


        time[5] += System.nanoTime() - tempTime;
        tempTime = System.nanoTime();
        tempTime = System.nanoTime();
        RealMatrix standardDeviation = calculator.calculateStandardDeviation(matrix);
        time[6] += System.nanoTime() - tempTime;
        //tempTime = System.nanoTime();
        //RealMatrix geomMean = calculator.calculateGeomMean(matrix);
        //time[7] += System.nanoTime() - tempTime;
        RealMatrix medianDeviation = calculator.calculateMedianDeviation(matrix);
        time[7] += System.nanoTime() - tempTime;
        tempTime = System.nanoTime();
        RealMatrix variationsCoefficient = calculator.calculateVariationsCoefficient(matrix);
        time[8] += System.nanoTime() - tempTime;
        tempTime = System.nanoTime();
        RealMatrix quartilsDispersion = calculator.calculateQuartilsDispersion(matrix);

        time[9] += System.nanoTime() - tempTime;

        for(int l = 0; l < matrixList.size(); l++){
            //standard metriken
            tempTime = System.nanoTime();
            listOfBestAtt.get(0).add(calculator.calculateMinForVectors(numberOfShownAttributes,
                    calculator.calculateVariance(matrixList.get(l))));
            time[0] += System.nanoTime() - tempTime;
            tempTime = System.nanoTime();
            listOfBestAtt.get(1).add(calculator.calculateMinForVectors(numberOfShownAttributes,
                    calculator.calculateStandardDeviation(matrixList.get(l))));
            time[1] += System.nanoTime() - tempTime;
            /*tempTime = System.nanoTime();
            listOfBestAtt.get(2).add(calculator.calculateMinForVectors(numberOfShownAttributes,
                    calculator.calculateGeomMean(matrixList.get(l))));
            time[2] += System.nanoTime() - tempTime;*/
            tempTime = System.nanoTime();
            listOfBestAtt.get(2).add(calculator.calculateMinForVectors(numberOfShownAttributes,
                    calculator.calculateMedianDeviation(matrixList.get(l))));
            time[2] += System.nanoTime() - tempTime;
            tempTime = System.nanoTime();
            listOfBestAtt.get(3).add(calculator.calculateMinForVectors(numberOfShownAttributes,
                    calculator.calculateQuartilsDispersion(matrixList.get(l))));
            time[3] += System.nanoTime() - tempTime;
            tempTime = System.nanoTime();
            listOfBestAtt.get(4).add(calculator.calculateMinForVectors(numberOfShownAttributes,
                    calculator.calculateVariationsCoefficient(matrixList.get(l))));
            time[4] += System.nanoTime() - tempTime;




            List<RealMatrix> varianceList2 = new ArrayList<>();
            List<RealMatrix> stdabwList = new ArrayList<>();
            List<RealMatrix> medDevList = new ArrayList<>();
            List<RealMatrix> varCofList = new ArrayList<>();
            List<RealMatrix> quartDispList = new ArrayList<>();
            List<RealMatrix> geomMeanList = new ArrayList<>();

            //Berechnung für metrikdifferenzen
            for(int x = 0; x < matrixList.size(); x++){

                tempTime = System.nanoTime();
                RealMatrix tempVariance = calculator.calculateVariance(matrixList.get(x));
                //GeneralCalculation.printMatrix(variance, "");
                //GeneralCalculation.printMatrix(tempVariance, "");


                time[5] += System.nanoTime() - tempTime;
                tempTime = System.nanoTime();
                RealMatrix tempStdabw = calculator.calculateStandardDeviation(matrixList.get(x));
                time[6] += System.nanoTime() - tempTime;
                tempTime = System.nanoTime();
                RealMatrix tempMedDev = calculator.calculateMedianDeviation(matrixList.get(x));
                time[7] += System.nanoTime() - tempTime;
                tempTime = System.nanoTime();
                RealMatrix tempVarCof = calculator.calculateVariationsCoefficient(matrixList.get(x));
                time[8] += System.nanoTime() - tempTime;
                tempTime = System.nanoTime();
                RealMatrix tempQuartDisp = calculator.calculateQuartilsDispersion(matrixList.get(x));
                time[9] += System.nanoTime() - tempTime;
                //tempTime = System.nanoTime();
                //RealMatrix tempGeomMean = calculator.calculateGeomMean(matrixList.get(x));
                //time[11] += System.nanoTime() - tempTime;

                for(int y = 0; y < tempVariance.getRowDimension(); y++){
                    tempTime = System.nanoTime();
                    tempVariance.setEntry(y,0, - (Math.abs(variance.getEntry(y,0)) - Math.abs(tempVariance.getEntry(y,0))));
                    time[5] += System.nanoTime() - tempTime;
                    tempTime = System.nanoTime();
                    tempStdabw.setEntry(y,0, -(Math.abs(standardDeviation.getEntry(y,0)) - Math.abs(tempStdabw.getEntry(y,0))));
                    time[6] += System.nanoTime() - tempTime;
                    tempTime = System.nanoTime();
                    tempMedDev.setEntry(y,0, -(Math.abs(medianDeviation.getEntry(y,0)) - Math.abs(tempMedDev.getEntry(y,0))));
                    time[7] += System.nanoTime() - tempTime;
                    tempTime = System.nanoTime();
                    tempVarCof.setEntry(y,0, -(Math.abs(variationsCoefficient.getEntry(y,0)) - Math.abs(tempVarCof.getEntry(y,0))));
                    time[8] += System.nanoTime() - tempTime;
                    tempTime = System.nanoTime();
                    tempQuartDisp.setEntry(y,0, -(Math.abs(quartilsDispersion.getEntry(y,0)) - Math.abs(tempQuartDisp.getEntry(y,0))));
                    time[9] += System.nanoTime() - tempTime;
                    //tempTime = System.nanoTime();
                    //tempGeomMean.setEntry(y,0, -(Math.abs(geomMean.getEntry(y,0)) - Math.abs(tempGeomMean.getEntry(y,0))));
                    //time[11] += System.nanoTime() - tempTime;
                }
                varianceList2.add(tempVariance);
                //GeneralCalculation.printMatrix(tempVariance, "");
                stdabwList.add(tempStdabw);
                medDevList.add(tempMedDev);
                varCofList.add(tempVarCof);
                quartDispList.add(tempQuartDisp);
                //geomMeanList.add(tempGeomMean);
            }

            tempTime = System.nanoTime();
            listOfBestAtt.get(5).add(calculator.calculateMinForVectors(numberOfShownAttributes, varianceList2.get(l)));
            time[5] += System.nanoTime() - tempTime;
            tempTime = System.nanoTime();
            listOfBestAtt.get(6).add(calculator.calculateMinForVectors(numberOfShownAttributes, stdabwList.get(l)));
            time[6] += System.nanoTime() - tempTime;
            tempTime = System.nanoTime();
            listOfBestAtt.get(7).add(calculator.calculateMinForVectors(numberOfShownAttributes, medDevList.get(l)));
            time[7] += System.nanoTime() - tempTime;
            tempTime = System.nanoTime();
            listOfBestAtt.get(8).add(calculator.calculateMinForVectors(numberOfShownAttributes, varCofList.get(l)));
            time[8] += System.nanoTime() - tempTime;
            tempTime = System.nanoTime();
            listOfBestAtt.get(9).add(calculator.calculateMinForVectors(numberOfShownAttributes, quartDispList.get(l)));
            time[9] += System.nanoTime() - tempTime;
            //tempTime = System.nanoTime();
            //listOfBestAtt.get(11).add(calculator.calculateMinForVectors(numberOfShownAttributes, geomMeanList.get(l)));
            //time[11] += System.nanoTime() - tempTime;

            for(int u = 0; u < listOfBestAtt.size(); u++){
                tempBestAttributes = new ArrayList<>(knownBestAttributes);
                int good = 0;
                for(int s = 0; s < listOfBestAtt.get(u).get(l).length; s++){
                    int att = (int)listOfBestAtt.get(u).get(l)[s][0];
                    if(knownBestAttributes.contains(att)){
                        good++;
                        count[u]++;
                        for(int z = 0; z < tempBestAttributes.size(); z++){
                            if(tempBestAttributes.get(z) == att){
                                tempBestAttributes.remove(z);
                                break;
                            }
                        }
                    }
                }
                for(int f = 0; f < tempBestAttributes.size(); f++){
                    notFound[tempBestAttributes.get(f)][u]++;
                }
                //MedDevDifference
                if(u == 7){
                    clusterGood[0][l] = good;
                    //VarDiff
                } else if(u == 5){
                    clusterGood[1][l] = good;
                    //MedDev
                }else if(u == 2){
                    clusterGood[2][l] = good;
                    //Var
                }else if(u == 0){
                    clusterGood[3][l] = good;
                }
            }

        }
        matrixList.clear();
        try{
            generateOutputFile1(count, listOfBestAtt, time, clusterNr, numberOfShownAttributes, points, features, noise, jaccard, mutual, rand);
            generatePlotFile(clusterGood, clusterNr, points, features, noise);
        }catch(Exception e){
            e.printStackTrace();
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

    private static List<List<Point>> stringsToCluster(String[][] input, int numberOfClusters){
        List<List<Point>> clusters = new ArrayList<>();
        for(int i = 0; i < numberOfClusters; i++){
            clusters.add(new ArrayList<Point>());
        }
        for(int k = 1; k < input.length -1; k++){
            int len = input[k].length - 2;
            double[] attList = new double[len];
            for(int l = 1; l < input[k].length -1; l++){
                attList[l - 1] = Double.parseDouble(input[k][l]);
            }
            Point p = new Point(attList.length, attList);
            clusters.get(Integer.parseInt(input[k][input[k].length-1])).add(p);
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

    private static void generateOutputFile(double[][] count, long[] time, char distribution, int numberOfClusters,
                                           int numberOfShownAttributes, int numberOfPointsPerCluster, int iterations,
                                           int steps)throws FileNotFoundException {
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

    private static void generateOutputFile1(double[][] count, double[][] notFound, long[] time, char distribution, int numberOfClusters,
                                            int numberOfShownAttributes, int numberOfPointsPerCluster, int iterations,
                                            int steps)throws FileNotFoundException{
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
        sb.append('\n');
        sb.append('\n');
        sb.append("Not Found:");
        sb.append('\n');
        sb.append("Attribute");
        sb.append('\n');
        for(int k = 0; k < notFound.length; k++){
            if((k % 4) != 3){
                if(k % 4 == 0){
                    sb.append('\n');
                }
                sb.append(k % 4);
            }else{
                sb.append(6);
            }

            sb.append(',');
            for(int l = 0; l < notFound[k].length; l++){
                sb.append(notFound[k][l]);
                sb.append(',');
            }
            sb.append('\n');
        }

        csvWriter.write(sb.toString());
        csvWriter.close();
        System.out.println("Done");
    }


    private static void generateOutputFile1(double[] count, List<List<double[][]>> listOfBestAtt, long[] time, int numberOfClusters,
                                            int numberOfShownAttributes, int points, int features, double noise,
                                            double jaccard, double mutual, double rand)throws FileNotFoundException{
        double divisor = numberOfClusters*numberOfShownAttributes;
        StringBuilder sb = new StringBuilder();
        sb.append("new_eval");
        sb.append("_");
        sb.append("g");
        sb.append("_");
        sb.append(points);
        sb.append("_");
        sb.append(features);
        sb.append("_");
        sb.append(numberOfClusters);
        sb.append("_");
        sb.append(noise);
        sb.append(".csv");
        PrintWriter csvWriter = new PrintWriter(new File(sb.toString()));

        sb = new StringBuilder();

        sb.append("Metrik:");
        sb.append(',');
        sb.append("var");
        sb.append(',');
        sb.append("stdabw");
        sb.append(',');
        //sb.append("geoMean");
        //sb.append(',');
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
        //sb.append(',');
        //sb.append("geoMeanDiff");
        sb.append('\n');


        sb.append("Mean:");
        sb.append(',');
        sb.append(count[0]/(divisor));
        sb.append(',');
        sb.append(count[1]/(divisor));
        sb.append(',');
        sb.append(count[2]/(divisor));
        sb.append(',');
        sb.append(count[3]/(divisor));
        sb.append(',');
        sb.append(count[4]/(divisor));
        sb.append(',');
        sb.append(count[5]/(divisor));
        sb.append(',');
        sb.append(count[6]/(divisor));
        sb.append(',');
        sb.append(count[7]/(divisor));
        sb.append(',');
        sb.append(count[8]/(divisor));
        sb.append(',');
        sb.append(count[9]/(divisor));
        sb.append('\n');
        sb.append("cluster-statistic");
        sb.append(jaccard);
        sb.append(',');
        sb.append(mutual);
        sb.append(',');
        sb.append(rand);
        sb.append(',');
        sb.append('\n');
        sb.append("Time in ns:");
        sb.append(',');
        for(int i = 0; i < time.length; i++){
            sb.append(time[i]);
            sb.append(',');
        }
        sb.append('\n');
        sb.append('\n');
        sb.append("Found:");
        sb.append('\n');

        for(int k = 0; k < listOfBestAtt.get(0).size(); k++){
            sb.append("Attribute: ");
            sb.append(',');
            for(int l = 0; l < listOfBestAtt.size(); l++){
                for(int m = 0; m < listOfBestAtt.get(l).get(k).length; m++){
                    sb.append(listOfBestAtt.get(l).get(k)[m][0]);
                    sb.append('|');
                }
                sb.append(',');
            }
            sb.append('\n');
            sb.append("Metrik: ");
            sb.append(',');
            for(int l = 0; l < listOfBestAtt.size(); l++){
                for(int m = 0; m < listOfBestAtt.get(l).get(k).length; m++){
                    sb.append(listOfBestAtt.get(l).get(k)[m][1]);
                    sb.append('|');
                }
                sb.append(',');
            }
            sb.append('\n');
            sb.append('\n');
        }

        csvWriter.write(sb.toString());
        csvWriter.close();
        System.out.println("Done");
    }

    private static void generatePlotFile(int[][] clusterGood, int numberOfClusters, int points, int features, double noise)throws FileNotFoundException{
        StringBuilder sb = new StringBuilder();
        sb.append("./plotData/plot");
        sb.append("_");
        sb.append("g");
        sb.append("_");
        sb.append(points);
        sb.append("_");
        sb.append(features);
        sb.append("_");
        sb.append(numberOfClusters);
        sb.append("_");
        sb.append(noise);
        sb.append(".csv");
        PrintWriter csvWriter = new PrintWriter(new File(sb.toString()));

        sb = new StringBuilder();
        for(int k = 0; k < clusterGood.length; k++){
            for(int i : clusterGood[k]){
                sb.append(',');
                sb.append(i);
            }
            sb.append('\n');
        }


        csvWriter.write(sb.toString());
        csvWriter.close();
        System.out.println("Done");
    }

}
