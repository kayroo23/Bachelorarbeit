import java.util.ArrayList;
import java.util.List;

abstract class generateData {
    List<List<Point>> clusters = new ArrayList<>();


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
