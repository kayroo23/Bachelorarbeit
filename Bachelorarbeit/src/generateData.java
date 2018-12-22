import org.apache.commons.math3.distribution.AbstractRealDistribution;

import java.util.ArrayList;
import java.util.List;
import java.util.Random;

abstract class generateData {
    List<List<Point>> clusters = new ArrayList<>();

    Point getRandomPoint(AbstractRealDistribution randomGen1, AbstractRealDistribution randomGen2,
                         AbstractRealDistribution randomGenAge, int clusterNr, int numberOfClusters){
        Random randomNr = new Random();
        double[] attList = new double[8];
        int maxValue = numberOfClusters * 10;
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

        return new Point(attList.length, attList);
    }

    /**Generates the y Value of an logistic function
     *
     * @param x Zähler
     * @param maxX Nenner
     * @return y Value
     */
    static double getLogisticFunctionValue(double x, double maxX){
        double t= -2 + (5*(x/maxX));
        return (1/(1+Math.exp(-t)));
    }

    List<List<Point>> getClusters(){
        return clusters;
    }
}
