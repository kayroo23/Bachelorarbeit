import org.apache.commons.math3.distribution.AbstractRealDistribution;
import org.apache.commons.math3.distribution.CauchyDistribution;
import org.apache.commons.math3.distribution.NormalDistribution;
import org.apache.commons.math3.distribution.UniformRealDistribution;

import java.util.ArrayList;
import java.util.Random;

class generateSimpleData extends generateData{

    /**Generates a fakeData goldtandard
     *
     *
     * @param c     char to choose the distribution type
     * @param overlap   noisecoefficient
     * @param numberOfClusters number of clusters
     * @param pointsPerCluster points per cluster
     */
    generateSimpleData(char c, double overlap, int numberOfClusters, int pointsPerCluster){
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
                Point p = getRandomPoint(randomGen1, randomGen2, randomGenAge, clusterNr, numberOfClusters);
                clusters.get(clusterNr).add(p);
            }
        }
    }
}
