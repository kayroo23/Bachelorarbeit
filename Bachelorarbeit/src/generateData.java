import org.apache.commons.math3.distribution.AbstractRealDistribution;
import org.apache.commons.math3.distribution.CauchyDistribution;
import org.apache.commons.math3.distribution.NormalDistribution;
import org.apache.commons.math3.distribution.UniformRealDistribution;

import java.util.ArrayList;
import java.util.List;
import java.util.Random;

class generateData {
    private List<List<Point>> clusters = new ArrayList<>();

    /**Generates a fakeData goldtandard
     *
     *
     * @param c     char to choose the distribution type
     * @param overlap   noisecoefficient
     * @param numberOfClusters number of clusters
     * @param pointsPerCluster points per cluster
     */
    generateData(char c, double overlap, int numberOfClusters, int pointsPerCluster){
        AbstractRealDistribution randomGenAge;
        AbstractRealDistribution randomGen1;
        AbstractRealDistribution randomGen2;

        //initialisierung der Liste
        for(int i = 0; i < numberOfClusters; i++){
            clusters.add(new ArrayList<Point>());
        }
        int maxValue = numberOfClusters * 10;

        for(int clusterNr = 0; clusterNr < numberOfClusters; clusterNr++){
            Random randomNr = new Random();
            double overlap1 = 1 + (overlap * Math.random());
            double min = (10*clusterNr) - (overlap1*10);
            double max = 10 + (overlap1*10)+(10*clusterNr);

            //Auswahl der Verteilung
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

            //Hinzufuegen der eigentlichen Daten mit der obigen Verteilung
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

    List<List<Point>> getClusters(){
        return clusters;
    }
}
