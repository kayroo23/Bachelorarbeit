import org.apache.commons.math3.distribution.*;

import java.util.ArrayList;
import java.util.Random;

class multivariateGaussian extends generateData{

    multivariateGaussian(final double overlap, final int numberOfClusters, final int pointsPerCluster){

        //initialisierung der Liste
        for(int i = 0; i < numberOfClusters; i++){
            clusters.add(new ArrayList<Point>());
        }
        int maxValue = numberOfClusters * 10;
        double[] sample;

        for(int clusterNr = 0; clusterNr < numberOfClusters; clusterNr++){
            Random randomNr = new Random();
            double overlap1 = 1 + (overlap * Math.random());
            double min = (10*clusterNr) - (overlap1*10);
            double max = 10 + (overlap1*10)+(10*clusterNr);

            double[] means = {(min+max)/2, (min+max)/2, getLogisticFunctionValue(clusterNr + 0.5, numberOfClusters) *  (10*clusterNr)};
            double[][] covariances = new double[3][3];
            for(int i = 0; i < covariances.length; i++){
                for(int j = i; j < covariances[i].length; j++){
                    //maximale Abhängigkeit der relevanten Attribute
                    covariances[i][j] = 1;
                }
            }

            AbstractMultivariateRealDistribution randomGen = new MultivariateNormalDistribution(means, covariances);

            //Hinzufuegen der eigentlichen Daten mit der obigen Verteilung
            for(int k = 0; k < pointsPerCluster; k++){
                double[] attList = new double[8];
                //good attributes
                sample = randomGen.sample();
                attList[0] = sample[0];
                attList[1] = sample[1];
                attList[2] = sample[2];
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
    }
}
